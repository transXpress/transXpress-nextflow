#!/usr/bin/env nextflow

/*
 * Notes: 
 * - must avoid spaces in filenames in samples.txt
 * - samples.txt must contain complete (absolute) paths to read files
 */

//This ensures that if the workflow is rerun on a different day from the same input files, then 
//it doesn't recalculate the assembly_prefix
def outFile = new File('transXpress_results/assembly_prefix.txt')
if ( ! outFile.exists() )
{

def resultsDir = new File('transXpress_results') 
if ( ! resultsDir.exists())
{
resultsDir.mkdirs() //Make the transXpress_results directory
}


theDate = ""
if ( params.prefix_add_date == true) {
theDate = new java.util.Date().format( params.prefix_add_date_formatting ) //yyMMdd by default
}

metadata = ""
if ( params.prefix_add_metadata_file != "" ) {
theText = file(params.prefix_add_metadata_file).text
metadata = theText.replace(" ","_").trim()
}

dateMetadataToJoin = []
if (theDate != "") {
dateMetadataToJoin.add(theDate)
}
if (metadata != "") {
dateMetadataToJoin.add(metadata)
}

dateMetadataPrefix = dateMetadataToJoin.join("_")
if (dateMetadataPrefix != "") {
  dateMetadataPrefix += "_"
}


outFile.withWriter('UTF-8') { writer ->
    writer.write(dateMetadataPrefix)
}
}

//Convert into a global string variable
dateMetadataPrefix = outFile.text

log.info """
 transXpress
 ===================================
 """+dateMetadataPrefix+" assembling with "+params.assembler

/*
 * Step 0 Download reference databases
 */

process downloadEggNOG {
  executor 'local'
  storeDir params.storeDB
  output:
    file "NOG.annotations.tsv" into eggNOGDb
  script:
    """
    wget "http://eggnogdb.embl.de/download/latest/data/NOG/NOG.annotations.tsv.gz"
    gunzip NOG.annotations.tsv.gz
    """
}

process downloadVirusesUniref50 {
  executor 'local'
  storeDir params.storeDB
  errorStrategy 'ignore'
  output:
    set file("virusesUniref50.pep.fasta"), file("virusesUniref50.pep.fasta.p??") into virusDb
  script:
    """
    wget -t 3 -O virusesUniref50.pep.fasta.gz "https://www.uniprot.org/uniref/?query=uniprot%3A%28taxonomy%3A%22Viruses+%5B10239%5D%22%29+AND+identity%3A0.5&format=fasta&compress=yes"
    gunzip virusesUniref50.pep.fasta.gz
    makeblastdb -in virusesUniref50.pep.fasta -dbtype prot
    """
}

process downloadRfam {
  executor 'local'
  storeDir params.storeDB
  output:
    set file("Rfam.cm"), file("Rfam.cm.???") into rfamDb
  script:
    """
    wget "ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz"
    gunzip Rfam.cm.gz
    cmpress Rfam.cm
    """
}

process downloadSprot {
  executor 'local'
  storeDir params.storeDB
  output:
    set file("uniprot_sprot.fasta"), file("uniprot_sprot.fasta.p??") into sprotDb
  script:
    """
    wget "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
    gunzip uniprot_sprot.fasta.gz
    makeblastdb -in uniprot_sprot.fasta -dbtype prot
    """
}

process downloadPfam {
  executor 'local'
  storeDir params.storeDB
  output:
    set file("Pfam-A.hmm"), file("Pfam-A.hmm.h??") into pfamDb
  script:
    """
    wget "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
    gunzip Pfam-A.hmm.gz
    hmmpress Pfam-A.hmm
    """
}

/*
 * Step 1. De novo assembly
 */


///
/// Load read files into channels
///

process loadSamples {
executor 'local'
input:
  file "samples.txt" from file(params.samples)
output:
  file "samples.txt" into toParse, toRelative
script:
"""
##do nothing. This process is just so the dag looks nicer
"""
}
//Load samples into filepair tuples.
toParse.splitCsv(sep:'\t',header:false)
     .map{ row ->
     println row
     return tuple(file(row[2]), file(row[3])) }
     .set{ readPairs_ch }

readPairs_ch.into{ trimReadPairs_ch ; fastqcReadPairs_ch }

//Once the sample files are loaded into Nextflow channels, everything should be specified relatively
process convertSamplesToRelative {
executor 'local'
input:
    file toRelative
output:
    file "relative_samples.txt" into relativeSamples_ch
script:
"""
#!/usr/bin/env python
import re
import os
import os.path

output_handle = open("relative_samples.txt", "w")

with open("${toRelative}", "r") as input_handle:
      for line in input_handle:
        row = re.split("[\t ]+", line)
        newRow = [row[0],row[1]]
        forwardReads = os.path.basename(row[2]).strip()+".R1-P.qtrim.fastq.gz"        
        reverseReads = os.path.basename(row[3]).strip()+".R2-P.qtrim.fastq.gz"        
        newRow.append(forwardReads)
        newRow.append(reverseReads)
        newRowString = "\t".join(newRow)+os.linesep
        output_handle.write(newRowString)
output_handle.close()        
"""
}
relativeSamples_ch.into{ samples_file_toTrinity; relativeSamples_toTrinityFinish; relativeSamples_toKallisto; samples_file_toYAMLConvert}


process trimmomatic {
cpus 4
input:
 set file(R1_reads),file(R2_reads) from trimReadPairs_ch
 file "adapters.fasta" from file(params.trimmomatic_adapter_file)
tag {"$R1_reads"+" and " +"$R2_reads"}
output:
  set file("${R1_reads}.R1-P.qtrim.fastq.gz"), file("${R2_reads}.R2-P.qtrim.fastq.gz") into filteredPairedReads_toChoice,filteredPairedReads_toKallisto
  //file "*U.qtrim.fastq.gz" into filteredSingleReads
script:
"""
##TODO:Adjust the params.TRIMMOMATIC_PARAMS in the nextflow.config file to change the parameters  
trimmomatic PE -threads ${task.cpus} ${R1_reads} ${R2_reads} ${R1_reads}.R1-P.qtrim.fastq.gz ${R1_reads}.R1-U.qtrim.fastq.gz ${R2_reads}.R2-P.qtrim.fastq.gz ${R2_reads}.R2-U.qtrim.fastq.gz ${params.TRIMMOMATIC_PARAMS} 
"""
}

///This switch controls which assembler is executed.
filteredPairedReads_toTrinity = Channel.create()
filteredPairedReads_toRnaspades = Channel.create()
filteredPairedReads_toChoice.choice(filteredPairedReads_toTrinity,filteredPairedReads_toRnaspades) { assemblerString = params.assembler.toLowerCase().trim()
                                                                                                     if (assemblerString == "trinity") return 0
                                                                                                     else if (assemblerString == "rnaspades") return 1
                                                                                                     else {
                                                                                                      println "Error:"+params.assembler+" doesn't match an implemented assembler"
                                                                                                      return -1 //Error status
                                                                                                      }
                                                                                                      }

filteredPairedReads_toTrinity.collect().into{ trinityInchwormPairedReads ; trinityFinishPairedReads }

process fastqc {
publishDir "transXpress_results/fastqc_results/", mode: "copy"
cpus 2
input:
 set file(R1_reads),file(R2_reads) from fastqcReadPairs_ch
tag {"$R1_reads"+" and " +"$R2_reads"}
output:
 set file("${R1_reads}.fastqc.ok"), file("${R2_reads}.fastqc.ok") into fastqcPassResults
 file "*.html" into fastqcHtmlResults
script:
"""
fastqc ${R1_reads} &
fastqc ${R2_reads}

##TODO: fix this, so the pipeline doesn't execute if bad data is detected.
##Check for bad run
#cat *.html | grep -oe "\\[FAIL\\].{1,30}Per sequence quality scores" > fastqc.fail
##If no bad run, produce the .fastqc.ok files

#if [[ ! -s fastqc.fail ]]
#then
# touch ${R1_reads}.fastqc.ok
# touch ${R2_reads}.fastqc.ok
#else
# echo "Error. FastQC detected a failure in the per read sequencing quality. This possibly indicates a bad sequencing run, which would result in a poor transcriptome assembly. Please solve the quality issue in your raw read files provided in samples.txt"
# exit 1
#fi

##Bash conditional wasn't working. Just doing a dummy:
touch ${R1_reads}.fastqc.ok
touch ${R2_reads}.fastqc.ok

sleep 15 ##Not a super important process, so might as well put in a delay to help with filesystem latency issues.
"""
}

//This is the samples file for rnaspades
process relativeSamplesToYAML {
executor 'local'
input:
    file samples_file_toYAMLConvert
    //file filteredPairedReads from filteredPairedReads_toYAML.collect() //collect flattens the tuple. This input ensures the process waits until trimmomatic is all done, and also allows for assertions as a sanity check
output:
    file "samples_trimmed.yaml" into yaml_rnaSPAdes_ch

script:
    """
    #!/usr/bin/env python
    import re
    import os
    import os.path
    import pprint

    sample_list = []
    with open("${samples_file_toYAMLConvert}", "r") as input_handle:
      for line in input_handle:
        row = re.split("[\t ]+", line)
        if (len(row) > 3): # paired reads
          paired_dict = {}
          paired_dict['orientation'] = 'fr'
          paired_dict['type'] = 'paired-end'
          f_reads = row[2].strip()
          r_reads = row[3].strip()
          #assert os.path.isfile(f_reads)
          #assert os.path.isfile(r_reads)
          paired_dict['left reads'] = [f_reads]
          paired_dict['right reads'] = [r_reads]
          ##Depends on where the unpaired reads are written to
          ##Presumably could get from sample txt as well
          ##Note: Commented out, to make it consistent with Trinity. May change in future
          ##unpaired_reads_f = f[:-16]+"U.qtrim.fastq.gz"
          ##unpaired_reads_r = r[:-16]+"U.qtrim.fastq.gz"
          ##assert os.path.isfile(unpaired_reads_f)
          ##assert os.path.isfile(unpaired_reads_r)
          ##unpaired_reads = [unpaired_reads_f] + [unpaired_reads_r]
          ##if len(unpaired_reads) > 0:
          ##    sample_dict['single reads'] = unpaired_reads
          sample_list.append(paired_dict)

        if (len(row) == 3): # unpaired reads
          unpaired_dict = {}
          unpaired_dict['type'] = 'single'
          u_reads = row[2].strip()
          assert os.path.isfile(u_reads)
          unpaired_dict['single reads'] = [u_reads]
          sample_list.append(unpaired_dict)

    with open("samples_trimmed.yaml", "w") as output_handle:
      output_handle.write(pprint.pformat(sample_list))
    """
}

process trinityInchwormChrysalis {
  cache 'lenient'

  cpus params.assembly_CPUs
  memory params.assembly_MEM+" GB"

  tag { dateMetadataPrefix+"Trinity" }

  afterScript 'echo \"(Above completion message is from Trinity. transXpress will continue the pipeline execution.)\"'
  //afterScript 'exit(1)'
  input:
     file trinityInchwormPairedReads //This flattens the tuple
     file samples_file from samples_file_toTrinity
  output:
    file "trinity_out_dir/[!Tcr]*" into trinityWorkDirRootFiles_ch1, trinityWorkDirRootFiles_ch2 //Not files starting with c or r, so not chrysalis, read_partitions, recursive trinity cmds, Trinity.timing
    file "trinity_out_dir/chrysalis/*" into trinityWorkDirChrysalisFiles_ch1, trinityWorkDirChrysalisFiles_ch2
    file "trinity_out_dir/read_partitions/Fb*/C*/*.trinity.reads.fa" into trinityPhase1ReadPartitionsFiles_ch
    //

  script:
    """
    Trinity --no_distributed_trinity_exec --max_memory ${task.memory.toGiga()}G --CPU ${task.cpus} --samples_file ${samples_file} ${params.TRINITY_PARAMS}
    sleep 5 ##Try to prevent filesystem latency errors
    #chmod -R a-w ./trinity_out_dir ##<- make the results read only, to troubleshoot Trinity processes writing into directories they shouldn't.
    """
}

//This is a workaround, of sorts, to parallelize the Trinity butterfly steps
//As Trinity usually expects to be writing into a single directory over time
//which doesn't play nice with Nextflow, have to do some tricks
//to recreate the structure of the Trinity work directory
//so each process ends up with a symbolic linked copy of the Trinity work directory
//and can write real files in the proper directory structure.
//But that directory structure can't be remade by Nextflow using Queue Channels, that I am aware of
//So the below mapping trick makes the directory structure into variables
//that the "trinityFinish" process can use to recreate the structure
trinityPhase1ReadPartitionsFiles_ch.flatten().map{ file ->
                                        def filePath = file.toAbsolutePath().toString()
                                        dir1Match = filePath =~ /Fb_[0-9]+/ ///Regex matching
                                        dir2Match = filePath =~ /CBin_[0-9]+/ ///Regex matching
					dir1String = dir1Match[0]
					dir2String = dir2Match[0]
 
                                        return tuple( [dir1String,dir2String] ,file)
                                        }
                                        .groupTuple()
                                        .set{ partitionedReadGroups_ch }
	//Add .groupTuple() to execute in groups by the directories

process trinityButterflyParallelVersion2 {
  //TODO This process has hardcoded parameters.  It should really be getting them from the TRINITY params...
  cache 'lenient'
  cpus 10
  input:
    file "trinity_out_dir/*" from trinityWorkDirRootFiles_ch1 //An attempt to relativize the butterfly processes
    file "trinity_out_dir/chrysalis/*" from trinityWorkDirChrysalisFiles_ch1 //An attempt to relativize the butterfly processes
    set dir,file(fastaFiles) from partitionedReadGroups_ch
    
  output:
   file "commands.completed" into trinityFinishedCmds
   file "commands.txt" into trinityCmds
   file "trinity_out_dir/read_partitions/*/*/*out.Trinity.fasta" into butterflyTrinityFiles
  tag { dateMetadataPrefix+"Trinity-"+dir[0]+"/"+dir[1]}
  script:
    """
    ##Have to recreate the directory structure for the read_partition files
    mkdir "trinity_out_dir/read_partitions"
    mkdir "trinity_out_dir/read_partitions/${dir[0]}"
    mkdir "trinity_out_dir/read_partitions/${dir[0]}/${dir[1]}"
    for f in ./*.fa
    do
     ##Once the directory structure is made, link the FASTA file back into it.
     ##Note the dollar sign escaping for nextflow
     ln -s "../../../../\$f" "trinity_out_dir/read_partitions/${dir[0]}/${dir[1]}/\$f"
    done
    for f in ./trinity_out_dir/read_partitions/${dir[0]}/${dir[1]}/*.fa
    do
     ##Make our own Trinity commands, using relative paths
     echo "Trinity --single '\$f' --output '\$f.out' --CPU 1 --max_memory 1G --run_as_paired --seqType fa --trinity_complete --full_cleanup --no_distributed_trinity_exec --min_glue 2 --min_kmer_cov 2" >> commands.txt
    done
    ##Execute in parallel
    parallel --jobs ${task.cpus} < commands.txt
    cp commands.txt commands.completed
    #chmod -R a-w ./trinity_out_dir
    sleep 15 ##Try to prevent filesystem latency errors
    """ 
}


process trinityFinish {
   publishDir "transXpress_results", mode: "copy", saveAs: { filename -> filename.replaceAll("trinity_out_dir/Trinity", "transcriptome") }
  input:
    file "trinity_out_dir/*" from trinityWorkDirRootFiles_ch2 //An attempt to relativize the butterfly processes
    file "trinity_out_dir/chrysalis/*" from trinityWorkDirChrysalisFiles_ch2 //An attempt to relativize the butterfly processes
    file butterflyTrinityFilesCollected from butterflyTrinityFiles.collect() 
    file relativeSamples from relativeSamples_toTrinityFinish
    file finishedCommands from trinityFinishedCmds.collectFile(name: "recursive_trinity.cmds.completed")
    file trinityCmdsCollected from trinityCmds.collectFile(name: "recursive_trinity.cmds")
    file trinityFinishPairedReads
  output:
    set val("Trinity"), file("./trinity_out_dir/Trinity.fasta.gene_trans_map"),file("./trinity_out_dir/Trinity.fasta") into trinityFinalOutput
    file "./trinity_out_dir/recursive_trinity.cmds.completed"
  memory "1 GB"
  tag { dateMetadataPrefix+"Trinity" }
  script:
    """
    mkdir trinity_out_dir/read_partitions/
    mkdir trinity_out_dir/read_partitions/Fb_0/
    mkdir trinity_out_dir/read_partitions/Fb_0/CBin_0 ##This is just a dummy directory to fool Trinity
    for f in ./*.fasta
    do
     ##Link the files into a directory structure that Trinity can deal with correctly, even if it isn't 100% right
     ln -s ../../../../\$f trinity_out_dir/read_partitions/Fb_0/CBin_0/
    done
    ##Have to produce these files to trick Trinity into thinking things are done
    cp ${trinityCmdsCollected} trinity_out_dir/recursive_trinity.cmds
    cp ${finishedCommands} trinity_out_dir/recursive_trinity.cmds.completed
    touch trinity_out_dir/recursive_trinity.cmds.ok
    Trinity --samples_file ${relativeSamples} --max_memory ${task.memory.toGiga()}G ${params.TRINITY_PARAMS}
    """ 
}


process runSPAdes {
cpus params.assembly_CPUs
memory params.assembly_MEM+" GB"

input:
   file filteredPairedReads from filteredPairedReads_toRnaspades.collect()
   file datasets_YAML from yaml_rnaSPAdes_ch
   //file filteredSingleReads from filteredSingleReads_ch2.collect()
output:
   set val("rnaSPAdes"), file("rnaSPAdes.gene_trans_map"),file(dateMetadataPrefix+"rnaSPAdes/transcripts.fasta") into rnaSPAdesFinalOutput
script:
"""
rnaspades.py --dataset ${datasets_YAML} ${params.STRAND_SPECIFIC_RNASPADES} -t ${task.cpus} -m ${task.memory.toGiga()} -o ${dateMetadataPrefix}rnaSPAdes --only-assembler -k 47
##Make a fake gene to transcript file:
cat "${dateMetadataPrefix}rnaSPAdes/transcripts.fasta" | grep ">" | tr -d ">" | cut -f 1 -d " " > tmp.txt
paste tmp.txt tmp.txt > rnaSPAdes.gene_trans_map
"""
}

trinityFinalOutput.mix(rnaSPAdesFinalOutput).set{ finishedAssemblies }

process renameAssembly {
   executor 'local'
   publishDir "transXpress_results", mode: "copy"
   input:
    set val(assembler), file(geneTransMap), file(transcriptome_fasta) from finishedAssemblies
    file "species.txt" from file(params.prefix_add_metadata_file) //Just a dummy input to include the file on the DAG
   output:
    file "${dateMetadataPrefix}${assembler}.transcripts.fasta" into transcriptomeToSplit 
    set assembler, file("${dateMetadataPrefix}${assembler}.transcripts.fasta") into transcriptomeToTransdecoder, transcriptomeToStats, transcriptomeToAnnotation //Also pass the assembler type along
    set assembler, file("${dateMetadataPrefix}${assembler}.transcripts.fasta"), file("${assembler}_renamed.fasta.gene_trans_map") into transcriptomeGeneTransMapKallisto
   tag { dateMetadataPrefix+"${assembler}" }
   script:
   """
   seqkit replace -p 'TRINITY_' -r '' ${transcriptome_fasta} | seqkit replace -p '^' -r '${dateMetadataPrefix}${assembler}_' > ${dateMetadataPrefix}${assembler}.transcripts.fasta 
   cat ${geneTransMap} | sed 's/TRINITY_//g' | sed 's/^/${dateMetadataPrefix}${assembler}_/g' | sed 's/\t/\t${dateMetadataPrefix}${assembler}_/g' > ${assembler}_renamed.fasta.gene_trans_map
   """
}

/*
 * Step 3. Annotate the assembly
 */

process transdecoderLongOrfs {
  publishDir "transXpress_results", mode: "copy"
  input:
    set val(assemblerName),file(transcriptomeTransdecoder) from transcriptomeToTransdecoder
  output:
    file "${transcriptomeTransdecoder}.transdecoder_dir/*.pep" into longOrfsProteomeSplit
    file "${transcriptomeTransdecoder}.transdecoder_dir/*" into transdecoderLongOrfsDirFiles
    file "${transcriptomeTransdecoder}.transdecoder_dir.__checkpoints_longorfs/*" into longOrfsCheckpointsFiles
    file "*.cmds" into longOrfsRootCmds
    set val("${assemblerName}"), file(transcriptomeTransdecoder) into transcriptomeTransdecoderPredict
  tag { dateMetadataPrefix+"${assemblerName}" }
  script:
    """
    TransDecoder.LongOrfs -t ${transcriptomeTransdecoder}
    #chmod -R a-w ${transcriptomeTransdecoder}.transdecoder_dir/ ##write protect the output to troubleshoot downstream accessing.
    """
}

process kallisto {
  publishDir "transXpress_results", mode: "copy"
  cpus params.assembly_CPUs
  input:
    file filteredReadsFromPairs from filteredPairedReads_toKallisto.collect() //This flattens the tuples
    set val(assemblerKallisto), file(transcriptomeKallisto), file(geneTransMap) from transcriptomeGeneTransMapKallisto
    file relativeSamples from relativeSamples_toKallisto
  output:
    file "kallisto.isoform.TPM.not_cross_norm" into rawKallistoTable
    file "kallisto.isoform.TMM.EXPR.matrix" optional true into normalizedKallistoTable
    file "kallisto.gene.TPM.not_cross_norm" 
    file "kallisto.gene.TMM.EXPR.matrix" optional true
  tag { dateMetadataPrefix+"${assemblerKallisto}" }
  script:
    """
    export TRINITY_HOME=\$(dirname \$(readlink -f \$(which Trinity)))
    echo TRINITY_HOME set to \${TRINITY_HOME}
    \${TRINITY_HOME}/util/align_and_estimate_abundance.pl --transcripts ${transcriptomeKallisto} ${params.STRAND_SPECIFIC_TRINITY} --seqType fq --samples_file ${relativeSamples} --prep_reference --thread_count ${task.cpus} --est_method kallisto --gene_trans_map ${geneTransMap}
    \${TRINITY_HOME}/util/abundance_estimates_to_matrix.pl --est_method kallisto --name_sample_by_basedir --gene_trans_map $geneTransMap */abundance.tsv
    """
}

// Produce a new channel that will emit the expression data file, preferrably TMM normalized.
// If normalized expression file is not available (e.g., a single sample was used), the channel
// will emit the raw TPM values.
normalizedKallistoTable
    .concat( rawKallistoTable )
    .first()
    .into { transcriptExpression; expressionStats }

process trinityStats {
  publishDir "transXpress_results", mode: "copy"
  cpus 1
  input:
    set val(assemblerStats), file(transcriptomeStats) from transcriptomeToStats
    file expressionStats
  output:
    file "transcriptome_stats.txt"
    file "transcriptome_exN50.plot.pdf"
  tag { dateMetadataPrefix+"${assemblerStats}" }
  script:
    """
    export TRINITY_HOME=\$(dirname \$(readlink -f \$(which Trinity)))
    echo TRINITY_HOME set to \${TRINITY_HOME}
    \${TRINITY_HOME}/util/TrinityStats.pl ${transcriptomeStats} > transcriptome_stats.txt
    \${TRINITY_HOME}/util/misc/contig_ExN50_statistic.pl ${expressionStats} ${transcriptomeStats} > transcriptome_exN50
    \${TRINITY_HOME}/util/misc/plot_ExN50_statistic.Rscript transcriptome_exN50
    """
}


transcriptomeToSplit
  .splitFasta(by: 100, file: true)
  .into { sprotBlastxChunks; rfamChunks }

longOrfsProteomeSplit
  .splitFasta(by: 100, file: true)
  .into { sprotBlastpChunks; transdecoderPfamChunks }


process sprotBlastxParallel {
  cpus 2
  input:
    file chunk from sprotBlastxChunks
    set sprotDb, sprotDbIndex from sprotDb
  output:
    file "blastx_out" into sprotBlastxResults
  tag { chunk }
  script:
    """
    echo blastx ${chunk} using database ${sprotDb}
    blastx -query ${chunk} -db ${sprotDb} -num_threads ${task.cpus} -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt "6 std stitle" -out blastx_out
    """
}


process sprotBlastpParallel {
  cpus 2
  input:
    file chunk from sprotBlastpChunks
    set sprotDb, sprotDbIndex from sprotDb
  output:
    file "blastp_out" into sprotBlastpResults
  tag { dateMetadataPrefix+"-"+chunk }
  script:
    """
    echo blastp ${chunk} using database ${sprotDb}
    blastp -query ${chunk} -db ${sprotDb} -num_threads ${task.cpus} -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt "6 std stitle" -out blastp_out
    """
}
sprotBlastxResults.collectFile(name: 'blastx_annotations.tsv').set { blastxResult }
sprotBlastpResults.collectFile(name: 'blastp_annotations.tsv').into { blastpForTransdecoder; blastpResult }


process pfamParallel {
  cpus 2
  input:
    file chunk from transdecoderPfamChunks
    set pfamDb, pfamDbIndex from pfamDb
  output:
    file "pfam_out" into pfamResults
    file "pfam_dom_out" into pfamDomResults
    set file("${chunk}"),file("${"pfam_dom_out"}") into revisePfamChunks
  tag { dateMetadataPrefix+chunk }
  script:
    """
    echo pfam ${chunk} using database ${pfamDb}
    hmmscan --cpu ${task.cpus} --domtblout pfam_dom_out --tblout pfam_out ${pfamDb} ${chunk}
    """
}
pfamResults.collectFile(name: 'pfam_annotations.txt').set { pfamResult }
pfamDomResults.collectFile(name: 'pfam_dom_annotations.txt').into { pfamDomResult ; pfamForTransdecoder ; pfamToGff3Doms}


process rfamParallel {
  cpus 2
  input:
    file chunk from rfamChunks
    set rfamDb, rfamDbIndex from rfamDb
  output:
    file "rfam_out" into rfamResults
    //file "rfam_dom_out" into rfamDomResults
  tag { chunk }
  script:
    """
    echo rfam ${chunk} using database ${rfamDb}
    cmscan -E 0.00001 --incE 0.00001 --rfam --cpu ${task.cpus} --tblout rfam_out ${rfamDb} ${chunk}
    """
}

process publishRfamResults {
  executor 'local'
  publishDir "transXpress_results", mode: "copy"
  input:
    file rfamResult from rfamResults.collectFile(name: 'rfam_annotations_unsorted.txt')
  output:
    file "rfam_annotations.txt" into rfamResultPub
  tag { dateMetadataPrefix }
  script:
  """
  cat ${rfamResult} | head -n 2 > header.txt
  cat ${rfamResult} | tail -n 9 > footer.txt
  cat header.txt ${rfamResult} footer.txt | grep -v "#" | sort -k3,3 -k15nr,15 | sort -u -k3,3 --merge | sort -k15nr,15 > rfam_annotations.txt
  """
}

process transdecoderPredict {
  publishDir "transXpress_results", mode: "copy" // , saveAs: { filename -> "transcriptome_after_predict.pep" }
  input:
    set val(assembler),file(transcriptomeFile) from transcriptomeTransdecoderPredict
    file "${transcriptomeFile}.transdecoder_dir/*" from transdecoderLongOrfsDirFiles
    file "${transcriptomeFile}.transdecoder_dir.__checkpoints_longorfs/*" from longOrfsCheckpointsFiles
    file longOrfsRootCmds 
    file blastpForTransdecoder
    file pfamForTransdecoder
  output:
    file "${transcriptomeFile}.transdecoder.pep" into predictProteome, predictProteomeSplitBy100
    file "${transcriptomeFile}.transdecoder.bed"
    file "${transcriptomeFile}.transdecoder.cds"
    file "${transcriptomeFile}.transdecoder.gff3" into transdecoderGFF3ToPfam, transdecoderGFF3ToLiftover
  tag { dateMetadataPrefix+"${assembler}" }
  script:
    """
    TransDecoder.Predict -t ${transcriptomeFile} --retain_pfam_hits ${pfamForTransdecoder} --retain_blastp_hits ${blastpForTransdecoder}
    """
}

predictProteome.into{ proteomeToAnnotation ; proteomeToPfamRevise }

predictProteomeSplitBy100
  .splitFasta(by: 100, file: true)
  .into{ tmhmmChunks ; deeplocChunks ; signalpChunks }

signalpChunks.into{ signalp4Chunks ; signalp5Chunks}

proteomeToPfamRevise.combine(revisePfamChunks).set{ combinedToPfamRevise }

process revisePfamResults {
input:
 set file(proteome),file(longOrfChunk),file(hmmerscanDomFile) from combinedToPfamRevise
output:
 file "test.txt"
tag { dateMetadataPrefix+"${longOrfChunk}" }
script:
"""
seqkit fx2tab -n -i ${longOrfChunk} | cut -f 1 > ids.txt
seqkit grep -f ids.txt ${proteome} > predict_subset.fasta
seqkit grep -s -f <(seqkit seq -s ${longOrfChunk}) predict_subset.fasta > common.fasta

predictLen=`cat predict_subset.fasta | grep ">" | wc -l`
commonLen=`cat common.fasta | grep ">" | wc -l`
if [ \$predictLen -eq \$commonLen ]
then
      echo "no changes between the longorfs and transdecoder predict versions"
      touch test.txt
else
      ##TODO, maybe delete this whole process. I was under the impression that 
      ##the transdecoder peptides would change after transdecoder predict
      ##But seems this isn't true?
      ##Or at least I haven't seen this get triggered w/ the test datasets
      echo "changes detected"
fi


"""
}

process pfamToGff3 {
publishDir "transXpress_results", mode: "copy"

input:
file pfamDomResult from pfamToGff3Doms
file refGFF3 from transdecoderGFF3ToPfam

output:
 file "pfam_domains.gff3" //TODO hook this up to 
tag { "$refGFF3" }
script:
"""
##Using wrf's pfam2gff script.
##Disabled for now until the dependency situation is figured out.
##python ../../../pfam2gff.py -g ${refGFF3} -i ${pfamDomResult} -T > pfam_domains.gff3
touch pfam_domains.gff3
"""
}


process signalp4Parallel {
  cpus 1
  input:
    file chunk from signalp4Chunks
  output:
    file "signalp_out" into signalp4Results
  tag { chunk }
  script:
    """
    if [ -f "/usr/local/bin/signalp" ]; then
     signalp -t ${params.SIGNALP_ORGANISMS} -f short ${chunk} > signalp_out
    else
    echo "Unable to find signalP 4, so making dummy file instead"
     touch signalp_out
    fi
    """
}

process signalp5Parallel {
  cpus 1
  input:
    file chunk from signalp5Chunks
  output:
    file "*signalp5" into signalp5Results
    file "*.gff3" into signalp5ResultsGff3
  tag { chunk }
  script:
    """
    
    ##Weng lab specific stuff for testing
    if [ -f "/lab/solexa_weng/testtube/signalp-5.0/bin/signalp" ]; then
    /lab/solexa_weng/testtube/signalp-5.0/bin/signalp -prefix "signalp5_"${chunk} -org ${params.SIGNALP_ORGANISMS} -format short -fasta ${chunk} -gff3
    else
    echo "Unable to find signalP 5, so making dummy files instead"
    touch blank.gff3
    touch blank.signalp5
    fi
    
    """
}
signalp5ResultsGff3.collectFile(name: 'signalp5_annotations.gff3').into{ signalp5ResultGff3ToAnnotate ; signalp5ResultGff3ToLiftover}


process deeplocParallel {
  input:
    file chunk from deeplocChunks
  output:
    file "${chunk}.out.txt" into deeplocResults
  tag { chunk }
  script:
    """
    if hash deeploc 2>/dev/null;
    then
    export MKL_THREADING_LAYER=GNU
    eexport PATH="/lab/solexa_weng/testtube/miniconda3/bin:$PATH"
    deeploc -f ${chunk} -o ${chunk}.out
    else
    echo "Unable to find deeploc, so making dummy files instead"
    touch ${chunk}.out.txt
    fi 
    """
}

//TODO: Update this to tmhmm.py, which is on conda?
process tmhmmParallel {
  cpus 1
  input:
    file chunk from tmhmmChunks
  output:
    file "tmhmm_out" into tmhmmResults
  tag { chunk }
  script:
    """
    if hash deeploc 2>/dev/null;
    then
    echo tmhmm ${chunk}
    tmhmm --short < ${chunk} > tmhmm_out
    else
    echo "Unable to find tmhmm, so making dummy files instead"
    touch tmhmm_out
    fi
    """
}
tmhmmResults.collectFile(name: 'tmhmm_annotations.tsv').into{ tmhmmResultToAnnotate ; tmhmmResultToGff3 }

process tmhmmMakeGff3 {
executor 'local'
input:
 file tmhmmResultToGff3
tag{ dateMetadataPrefix }
output:
 file "*.gff3" into tmhmmGff3ToLiftover, tmhmmGff3
script:
"""
#!/usr/bin/env python

import os
import re

write_handle = open("tmhmm.pep.gff3","w")
read_handle = open("${tmhmmResultToGff3}","r")
for line in read_handle.readlines():
    if line[0] == "#" or len(line.strip()) == 0:
        continue
    splitrow = re.split("\\t+",line)
    protein_id = splitrow[0]
    topology = splitrow[5]
    matches = re.findall("([oi][0-9]+-[0-9]+)",topology)
    for m in matches:
        type = m[0]
        start,end = m[1:].split("-")
        outstring = "\\t".join([protein_id,"tmhmm","transmembrane_region",start,end,".",".",".","."])+os.linesep
        write_handle.write(outstring)
read_handle.close()
write_handle.close()
"""

}


signalp5ResultGff3ToLiftover.mix(tmhmmGff3ToLiftover).set{ combinedGff3s }
transdecoderGFF3ToLiftover.combine(combinedGff3s.collectFile(name: 'collected.gff3')).set{ liftOverGffsTuple }

process liftoverPeptideGff3ToTranscript {
executor 'local'
input:
 set file(transdecoderGFF3ToLiftover), file(gff3FilesToLiftover) from liftOverGffsTuple

tag { "${gff3FilesToLiftover}"  }
output:
 file "liftover.gff3" into liftoverResults
script:
"""
#!/usr/bin/env python

import re

read_handle = open("${transdecoderGFF3ToLiftover}","r")
transcript_len_dict = dict()
cds_dict = dict()
protein_to_transcript = dict()
transcript_to_protein = dict()
for line in read_handle.readlines():
    if line[0] == "#" or len(line.strip()) == 0:
        continue
    splitrow = re.split("\t",line)
    if "gene" == splitrow[2].strip().lower():
        transcript_len_dict[splitrow[0]] = int(splitrow[4].strip()) ##Length of transcipt.
    if "CDS" == splitrow[2].strip().upper():
        transcript_id = splitrow[0]
        cds_start = int(splitrow[3].strip())
        cds_end = int(splitrow[4].strip())
        cds_strand = splitrow[6].strip()
        ##Note the extra backslash for nextflow
        protein_id = re.search("ID=cds\\.(.+);",splitrow[8]).group(1)
        print(protein_id)
        cds_dict[protein_id] = {"cds_start":cds_start,"cds_end":cds_end,"cds_strand":cds_strand}
        protein_to_transcript[protein_id] = transcript_id
        transcript_to_protein[transcript_id] = protein_id
read_handle.close()

write_handle = open("liftover.gff3","w")
read_handle = open("${gff3FilesToLiftover}","r")
for line in read_handle.readlines():
    if line[0] == "#" or len(line.strip()) == 0:
        continue
    splitrow = re.split("\t",line)
    protein_id = splitrow[0]
    transcript_id = protein_to_transcript[protein_id]
    feature_start = int(splitrow[3].strip())
    feature_stop = int(splitrow[4].strip())
    feature_strand = splitrow[6].strip()
    ##For peptide features, as far as I know they are always specificed either as "+" or "."

    if cds_dict[protein_id]["cds_strand"] == "+":
        new_start = cds_dict[protein_id]["cds_start"] + (feature_start * 3)-3
        new_stop = cds_dict[protein_id]["cds_start"] + (feature_stop * 3)-3
        new_strand = "+"
    elif cds_dict[protein_id]["cds_strand"] == "-":
        new_start = cds_dict[protein_id]["cds_end"] - (feature_stop * 3)+3 
        new_stop =  cds_dict[protein_id]["cds_end"] - (feature_start * 3)+3 
        new_strand = "-"
    else:
        ##Note the extra backslash for nextflow
        sys.stderr.write("Unknown strand\\n")
        exit(1)
        ##Note the extra backslash for nextflow
    outstring = "\\t".join([protein_to_transcript[protein_id],splitrow[1],splitrow[2],str(new_start),str(new_stop),splitrow[5],new_strand,splitrow[7],splitrow[8]])
    write_handle.write(outstring)
write_handle.close()
read_handle.close()
"""
}
liftoverResults.collectFile(name: 'liftover_results.gff3').set{ liftoverResult  }


// Collect parallelized annotations
process annotatedFasta {
  publishDir "transXpress_results", mode: "copy"
  input:
    set val(assemblyAnnotation),file(transcriptomeFile) from transcriptomeToAnnotation
    file proteomeFile from proteomeToAnnotation
    file kallistoFile from transcriptExpression
    file blastxResult
    file blastpResult 
    file pfamResult 
    file pfamDomResult 
    file liftoverResult
    file signalp4Result from signalp4Results.collectFile(name: 'signalp4_annotations.txt')
    file signalp5Result from signalp5Results.collectFile(name: 'signalp5_annotations.txt')
    file signalp5ResultGff3ToAnnotate
    file deeplocResult from deeplocResults.collectFile(name: 'deeploc_annotations.txt')
    file tmhmmResultToAnnotate
  output:
    //TODO: Fix this output names, so they are unique across different assemblers
    file "${dateMetadataPrefix}${assemblyAnnotation}.transcripts_annotated.fasta" into transcriptome_annotated_fasta_ch
    file "${dateMetadataPrefix}${assemblyAnnotation}.proteins_annotated.fasta" into transcriptome_annotated_pep_ch
    file "transcriptome_TPM_blast.csv"
    file "${liftoverResult}"
    file "${blastxResult}"
    file "${blastpResult}"
    file "${pfamResult}"
    file "${pfamDomResult}"
    file "${deeplocResult}"
    file "${tmhmmResultToAnnotate}"
  tag { dateMetadataPrefix+"${assemblyAnnotation}" }
  script:
    """
    #!/usr/bin/env python
    
    import re
    import csv
    import Bio.SeqIO
  
    ## Annotation maps: transcript id -> annotation
    expression_annotations = {}
    blastx_annotations = {}
    blastp_annotations = {}
    pfam_annotations = {}
    tmhmm_annotations = {}
    deeploc_annotations = {}

    ## Load kallisto results
    print ("Loading expression values from ${kallistoFile}")
    with open("${kallistoFile}") as input_handle:
      csv_reader = csv.reader(input_handle, delimiter='\t')
      columns = next(csv_reader)
      for row in csv_reader:
        expression_annotations[row[0]] = columns[1] + "=" + str(row[1])
        for i in range(2, len(columns)):
          expression_annotations[row[0]] += " " + columns[i] + "=" + str(row[i])

    ## Load blastx results
    print ("Loading blastx results from ${blastxResult}")
    with open("${blastxResult}") as input_handle:
      csv_reader = csv.reader(input_handle, delimiter='\t')
      for row in csv_reader:
        if (len(row) < 13): continue
        blastx_annotations[row[0]] = row[12] + ", e=" + str(row[10])

    ## Load blastp results
    print ("Loading blastp results from ${blastpResult}")
    with open("${blastpResult}") as input_handle:
      csv_reader = csv.reader(input_handle, delimiter='\t')
      for row in csv_reader:
        if (len(row) < 13): continue
        blastp_annotations[row[0]] = row[12] + ", e=" + str(row[10])

    ## Load pfam results
    print ("Loading pfam predictions from ${pfamResult}")
    with open("${pfamResult}") as input_handle:
      for line in input_handle:
        if (line.startswith("#")): continue
        row = re.split(" +", line, 22)
        if (len(row) < 23): continue
        pfam_annotations[row[2]] = row[1] + " " + row[18] + ", e=" + str(row[7])

    ## Load tmhmm results
    print ("Loading tmhmm predictions from ${tmhmmResultToAnnotate}")
    with open("${tmhmmResultToAnnotate}") as input_handle:
      csv_reader = csv.reader(input_handle, delimiter='\t')
      for row in csv_reader:
        if (len(row) < 6): continue
        tmhmm_annotations[row[0]] = row[2] + ", " + row[3] + ", " + row[4] + ", " + row[5]
    
    ## Load deeploc results
    print ("Loading deeploc predictions from ${deeplocResult}")
    with open("${deeplocResult}") as input_handle:
      csv_reader = csv.reader(input_handle, delimiter="\t")
      for row in csv_reader:
        if (len(row) < 2): continue
        deeploc_annotations[re.split("[\t ]+",row[0])[0]] = str(row[1]) ##re.split tries to pull out the fasta record ID
    
    ## Do the work
    print ("Annotating FASTA file ${transcriptomeFile}")
    with open("${transcriptomeFile}", 'r') as input_fasta_handle, open("${dateMetadataPrefix}${assemblyAnnotation}.transcripts_annotated.fasta", 'w') as output_fasta_handle:
      for record in Bio.SeqIO.parse(input_fasta_handle, "fasta"):
        transcript_id = record.id
        record.description = "TPM: " + expression_annotations.get(transcript_id)
        if transcript_id in blastx_annotations:
          record.description += "; blastx: " + blastx_annotations.get(transcript_id)
        Bio.SeqIO.write(record, output_fasta_handle, "fasta")
    
    print ("Annotating FASTA file ${proteomeFile}")
    with open("${proteomeFile}", 'r') as input_fasta_handle, open("${dateMetadataPrefix}${assemblyAnnotation}.proteins_annotated.fasta", 'w') as output_fasta_handle:
      for record in Bio.SeqIO.parse(input_fasta_handle, "fasta"):
        transcript_id = re.sub("\\.p[0-9]+\$", "", record.id)
        record.description = "transdecoder " + re.search("ORF type:[^,]+,score=[^,]+", record.description).group(0)
        if transcript_id in expression_annotations:
          record.description += "; TPM: " + expression_annotations.get(transcript_id)
        if record.id in blastp_annotations:
          record.description += "; blastp: " + blastp_annotations.get(record.id)
        if record.id in pfam_annotations:
          record.description += "; pfam: " + pfam_annotations.get(record.id)
        if record.id in tmhmm_annotations:
          record.description += "; tmhmm: " + tmhmm_annotations.get(record.id)
        if record.id in deeploc_annotations:
          record.description += "; deeploc: " + deeploc_annotations.get(record.id)
        Bio.SeqIO.write(record, output_fasta_handle, "fasta")

    print ("Generating transcriptome_TPM_blast.csv table")
    with open("${kallistoFile}", 'r') as input_csv_handle, open("transcriptome_TPM_blast.csv", 'w') as output_csv_handle:
      csv_reader = csv.reader(input_csv_handle, delimiter='\t')
      csv_writer = csv.writer(output_csv_handle, delimiter=',')
      csv_columns = next(csv_reader)
      csv_columns[0] = "transcript"
      for i in range(1, len(csv_columns)):
        csv_columns[i] = "TPM(" + csv_columns[i] + ")"
      csv_columns.append("blastx")
      csv_writer.writerow(csv_columns)
      for row in csv_reader:
        row.append(blastx_annotations.get(row[0], ""))
        csv_writer.writerow(row)
   
    """
}

transcriptome_annotated_fasta_ch.into { transcriptome_annotated_fasta_ch1 ; transcriptome_annotated_fasta_ch2 }
transcriptome_annotated_pep_ch.into { transcriptome_annotated_pep_ch1 ; transcriptome_annotated_pep_ch2 }

process final_checksum {
    executor 'local'
    publishDir "transXpress_results", mode: "copy"
    input:
        file "transcriptome_annotated.fasta" from transcriptome_annotated_fasta_ch1
        file "transcriptome_annotated.pep" from transcriptome_annotated_pep_ch1
    output:
        file "assembly_seq-dependent_checksums.txt"
    tag { dateMetadataPrefix }
    script:
    """
    echo -n "transcriptome_annotated.fasta:" > assembly_seq-dependent_checksums.txt
    seqkit sort -s transcriptome_annotated.fasta | grep -v ">" | md5sum | cut -f 1 -d ' '>> assembly_seq-dependent_checksums.txt
    echo -n "transcriptome_annotated.pep:" >> assembly_seq-dependent_checksums.txt
    seqkit sort -s transcriptome_annotated.pep | grep -v ">" | md5sum | cut -f 1 -d ' '>> assembly_seq-dependent_checksums.txt
    """    
}

