#!/usr/bin/env nextflow

/*
 * Notes: 
 * - must avoid spaces in filenames in samples.txt
 * - samples.txt must contain complete (absolute) paths to read files
 */


params.tempdir = "/lab/weng_scratch/tmp/"

theDate = new java.util.Date().format( 'MMdd' )

theText = file(params.species).text
genus = theText.split(" ")[0]
species = theText.split(" ")[1].trim()
assemblyPrefix = theDate+"_"+genus+"_"+species+"_Trinity"

log.info """
 transXpress
 ===================================
 """+assemblyPrefix


/*
 * Step 1. 
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
##just do nothing
"""
}
toParse.splitCsv(sep:'\t',header:false)
     .map{ row ->
     println row
     return tuple(file(row[2]), file(row[3])) }
     .set{ readPairs_ch }


process convertSamplesToRelative {
executor 'local'
input:
    file toRelative
output:
    file "relative_samples.txt" into relative_samples
script:
"""
while read LINE; do
  START=\$(echo "\$LINE" | cut -f 1,2)
  FORWARD=\$(echo "\$LINE" | cut -f 3 | sed -r 's/[\\/\\.].+\\///g')".R1-P.qtrim.fastq.gz"
  REVERSE=\$(echo "\$LINE" | cut -f 4 | sed -r 's/[\\/\\.].+\\///g')".R2-P.qtrim.fastq.gz"
  echo \$FORWARD \$REVERSE
  echo "\${START}\t\${FORWARD}\t\${REVERSE}" >> relative_samples.txt
done < samples.txt
"""
}
relative_samples.into{ samples_file_toTrinity; relative_samples_toTrinityFinish; relative_samples_toKallisto; samples_file_toYAMLConvert}

process trimmomatic {
cpus 4
input:
 set file(R1_reads),file(R2_reads) from readPairs_ch
tag {"$R1_reads"+" and " +"$R2_reads"}
output:
  set file("${R1_reads}.R1-P.qtrim.fastq.gz"), file("${R2_reads}.R2-P.qtrim.fastq.gz") into filteredPairedReads_toChoice,filteredPairedReads_toKallisto,filteredPairedReads_toYAML
  file "*U.qtrim.fastq.gz" into filteredSingleReads
script:
"""
java -jar /lab/solexa_weng/testtube/trinityrnaseq-Trinity-v2.8.4/trinity-plugins/Trimmomatic/trimmomatic.jar PE -threads ${task.cpus}  ${R1_reads} ${R2_reads} ${R1_reads}.R1-P.qtrim.fastq.gz ${R1_reads}.R1-U.qtrim.fastq.gz ${R2_reads}.R2-P.qtrim.fastq.gz ${R2_reads}.R2-U.qtrim.fastq.gz  ILLUMINACLIP:/lab/solexa_weng/testtube/trinityrnaseq-Trinity-v2.8.4/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25 
"""
}
filteredPairedReads_toTrinity = Channel.create()
filteredPairedReads_toRnaspades = Channel.create()
filteredPairedReads_toChoice.choice(filteredPairedReads_toTrinity,filteredPairedReads_toRnaspades) { "trinity" =~ /rinity/ ? 0 : 1 }

filteredPairedReads_toTrinity.collect().into{ trinityInchwormPairedReads ; trinityFinishPairedReads }

process relativeSamplesToYAML {
executor 'local'
input:
    file samples_file_toYAMLConvert
    file filteredPairedReads from filteredPairedReads_toYAML.collect() //collect flattens the tuple. This input ensures the process waits until trimmomatic is all done, and also allows for assertions as a sanity check
output:
    file "samples_trimmed.yaml" into yaml_samples_rnaspades_ch

script:
    """
    #!/usr/bin/env python
    import os
    import os.path
    import pprint
    print("yaml conversion running")
    read_handle = open("${samples_file_toYAMLConvert}","r")
    lines = read_handle.readlines()
    
    sample_list = []
    for l in lines:
        splitline = l.strip().split("\t")
        if len(splitline) == 1:
            splitline = l.strip().split(" ")
        f = splitline[2]
        r = splitline[3]
        assert os.path.isfile(f)
        assert os.path.isfile(r)
        sample_dict = dict()
        sample_dict['orientation'] = 'fr'
        sample_dict['type'] = 'paired-end'
        sample_dict['right reads'] = [f]
        sample_dict['left reads'] = [r]
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
        sample_list.append(sample_dict)
    write_handle = open("samples_trimmed.yaml","w")
    write_handle.write(pprint.pformat(sample_list))
    write_handle.close()
    read_handle.close()
    """
}


 


process trinityInchwormChrysalis {
  cache 'lenient'
  label = "nf_"+assemblyPrefix+"_trinityInchwormChrysalis"

  cpus 12
  memory "200 GB"

  tag { assemblyPrefix }

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
    #chmod -R a-w ./trinity_out_dir
    """
}

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
  tag { assemblyPrefix+"-"+dir[0]+"/"+dir[1]}
  script:
    """
    ##Have to recreate the directory structure for the read_parition files
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
    file relative_samples_txt from relative_samples_toTrinityFinish
    file finishedCommands from trinityFinishedCmds.collectFile(name: "recursive_trinity.cmds.completed")
    file trinityCmdsCollected from trinityCmds.collectFile(name: "recursive_trinity.cmds")
    file trinityFinishPairedReads
  output:
    file "./trinity_out_dir/Trinity.fasta.gene_trans_map" into originalGeneTransMap
    file "./trinity_out_dir/Trinity.fasta" into Trinity_fasta_ch
    file "./trinity_out_dir/recursive_trinity.cmds.completed"
  memory "1 GB"
  tag { assemblyPrefix }
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
    Trinity --samples_file ${relative_samples_txt} --max_memory ${task.memory.toGiga()}G ${params.TRINITY_PARAMS}
    """ 
}


process runSPAdes {
cpus 16
memory "200 GB"

input:
   file filteredPairedReads from filteredPairedReads_toRnaspades.collect()
   file datasets_YAML from yaml_samples_rnaspades_ch
   //file filteredSingleReads from filteredSingleReads_ch2.collect()
output:
   file assemblyPrefix+"/transcripts.fasta" into spadesAssembly_ch

script:
"""
rnaspades.py --dataset ${datasets_YAML} -t ${task.cpus} -m ${task.memory.toGiga()} --fast -o ${assemblyPrefix} --only-assembler -k 47
"""
}


process renameTrinityAssembly {
   publishDir "transXpress_results", mode: "copy"
   tag { assemblyPrefix }
   input:
    file "Trinity.fasta" from Trinity_fasta_ch
    file "species.txt" from file(params.species) //Just a dummy input
    file "Trinity.fasta.gene_trans_map" from originalGeneTransMap
   output:
    file assemblyPrefix+".fasta" into transcriptomeKallisto, transcriptomeTransdecoder, transcriptomeStats, transcriptomeSplit, transcriptomeAnnotation
    file "Trinity_renamed.fasta.gene_trans_map" into geneTransMap

   script:
   """
   seqkit replace -p '^TRINITY' -r '${assemblyPrefix}' Trinity.fasta > ${assemblyPrefix}.fasta 
   cat Trinity.fasta.gene_trans_map | sed 's/^TRINITY/${assemblyPrefix}/' | sed 's/\tTRINITY/\t${assemblyPrefix}/g' > Trinity_renamed.fasta.gene_trans_map
   """
}

process transdecoderLongOrfs {
  publishDir "transXpress_results", mode: "copy"
  tag { assemblyPrefix }
  input:
    file transcriptomeTransdecoder
  output:
    file "${transcriptomeTransdecoder}.transdecoder_dir/*.pep" into longOrfsProteomeSplit
    file "${transcriptomeTransdecoder}.transdecoder_dir/*" into transdecoderLongOrfsDirFiles
    file "${transcriptomeTransdecoder}.transdecoder_dir.__checkpoints_longorfs/*" into longOrfsCheckpointsFiles
    file "*.cmds" into longOrfsRootCmds
    set val("${transcriptomeTransdecoder}"), file(transcriptomeTransdecoder) into transcriptomeTransdecoderPredict
  script:
    """
    TransDecoder.LongOrfs -t ${transcriptomeTransdecoder}
    chmod -R a-w ${transcriptomeTransdecoder}.transdecoder_dir/ ##write protect the output to troubleshoot downstream accessing.
    """
}

process kallisto {
  publishDir "transXpress_results", mode: "copy"
  cpus 10
  tag { assemblyPrefix }
  input:
    file filteredReadsFromPairs from filteredPairedReads_toKallisto.collect() //This flattens the tuples
    file transcriptomeKallisto
    file geneTransMap
    file relative_samples_txt from relative_samples_toKallisto
  output:
    file "kallisto.isoform.TPM.not_cross_norm" into rawKallistoTable
    file "kallisto.isoform.TMM.EXPR.matrix" optional true into normalizedKallistoTable
    file "kallisto.gene.TPM.not_cross_norm" 
    file "kallisto.gene.TMM.EXPR.matrix" optional true
  script:
    """
    export TRINITY_HOME=\$(dirname `which Trinity`)
    echo TRINITY_HOME set to \${TRINITY_HOME}
    \${TRINITY_HOME}/util/align_and_estimate_abundance.pl --transcripts ${transcriptomeKallisto} ${params.STRAND_SPECIFIC} --seqType fq --samples_file ${relative_samples_txt} --prep_reference --thread_count ${task.cpus} --est_method kallisto --gene_trans_map ${geneTransMap}
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
  tag { assemblyPrefix }
  input:
    file transcriptomeStats
    file expressionStats
  output:
    file "transcriptome_stats.txt"
    file "transcriptome_exN50.plot.pdf"
  script:
    """
    export TRINITY_HOME=\$(dirname `which Trinity`)
    echo TRINITY_HOME set to \${TRINITY_HOME}
    \${TRINITY_HOME}/util/TrinityStats.pl ${transcriptomeStats} > transcriptome_stats.txt
    \${TRINITY_HOME}/util/misc/contig_ExN50_statistic.pl ${expressionStats} ${transcriptomeStats} > transcriptome_exN50
    \${TRINITY_HOME}/util/misc/plot_ExN50_statistic.Rscript transcriptome_exN50
    """
}


transcriptomeSplit
  .splitFasta(by: 100, file: true)
  .into { sprotBlastxChunks; rfamChunks }

longOrfsProteomeSplit
  .splitFasta(by: 100, file: true)
  .into { sprotBlastpChunks; pfamChunks }

process downloadPfam {
  executor 'local'
  storeDir '/lab/solexa_weng/tmp/db'
  output:
    set file("Pfam-A.hmm"), file("Pfam-A.hmm.h??") into pfamDb
  script:
    """
    wget "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
    gunzip Pfam-A.hmm.gz
    hmmpress Pfam-A.hmm
    """
}

process downloadVirusesUniref50 {
  executor 'local'
  storeDir '/lab/solexa_weng/tmp/db'
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
  storeDir '/lab/solexa_weng/tmp/db'
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
  storeDir '/lab/solexa_weng/tmp/db'
  output:
    set file("uniprot_sprot.fasta"), file("uniprot_sprot.fasta.p??") into sprotDb
  script:
    """
    wget "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
    gunzip uniprot_sprot.fasta.gz
    makeblastdb -in uniprot_sprot.fasta -dbtype prot
    """
}


process sprotBlastxParallel {
  cpus 2
  input:
    file chunk from sprotBlastxChunks
    set sprotDb, sprotDbIndex from sprotDb
  output:
    file "blastx_out" into sprotBlastxResults
  tag { assemblyPrefix+"-"+chunk }
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
  tag { assemblyPrefix+"-"+chunk }
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
    file chunk from pfamChunks
    set pfamDb, pfamDbIndex from pfamDb
  output:
    file "pfam_out" into pfamResults
    file "pfam_dom_out" into pfamDomResults
  tag { assemblyPrefix+"-"+chunk }
  script:
    """
    echo pfam ${chunk} using database ${pfamDb}
    hmmscan --cpu ${task.cpus} --domtblout pfam_dom_out --tblout pfam_out ${pfamDb} ${chunk}
    """
}
pfamResults.collectFile(name: 'pfam_annotations.txt').set { pfamResult }
pfamDomResults.collectFile(name: 'pfam_dom_annotations.txt').into { pfamDomResult ; pfamForTransdecoder }

process rfamParallel {
  cpus 2
  input:
    file chunk from rfamChunks
    set rfamDb, rfamDbIndex from rfamDb
  output:
    file "rfam_out" into rfamResults
    //file "rfam_dom_out" into rfamDomResults
  tag { assemblyPrefix+"-"+chunk }
  script:
    """
    echo rfam ${chunk} using database ${rfamDb}
    cmscan -E 0.00001 --incE 0.00001 --rfam --cpu ${task.cpus} --tblout rfam_out ${rfamDb} ${chunk}
    """
}

process publishRfamResults {
  publishDir "transXpress_results", mode: "copy"
  input:
    file rfamResult from rfamResults.collectFile(name: 'rfam_annotations_unsorted.txt')
  output:
    file "rfam_annotations.txt" into rfamResultPub
  tag { assemblyPrefix }
  script:
  """
  cat ${rfamResult} | head -n 2 > header.txt
  cat ${rfamResult} | tail -n 9 > footer.txt
  cat header.txt ${rfamResult} footer.txt | grep -v "#" | sort -k3,3 -k15nr,15 | sort -u -k3,3 --merge | sort -k15nr,15 > rfam_annotations.txt
  """
}

process transdecoderPredict {
  publishDir "transXpress_results", mode: "copy" // , saveAs: { filename -> "transcriptome_after_predict.pep" }
  tag { assemblyPrefix }
  input:
    set val(transcriptomeName),file(transcriptomeFile) from transcriptomeTransdecoderPredict
    file "${transcriptomeName}.transdecoder_dir/*" from transdecoderLongOrfsDirFiles
    file "${transcriptomeName}.transdecoder_dir.__checkpoints_longorfs/*" from longOrfsCheckpointsFiles
    file longOrfsRootCmds 
    file blastpForTransdecoder
    file pfamForTransdecoder
  output:
    file "${transcriptomeName}.transdecoder.pep" into predictProteome, predictProteomeSplitBy100,predictProteomeSplitBy10
  script:
    """
    TransDecoder.Predict -t ${transcriptomeFile} --retain_pfam_hits ${pfamForTransdecoder} --retain_blastp_hits ${blastpForTransdecoder}
    """
}

predictProteomeSplitBy100
  .splitFasta(by: 100, file: true)
  .set { tmhmmChunks }

predictProteomeSplitBy10
  .splitFasta(by: 10, file: true)
  .set{ deeplocChunks }

process deeplocParallel {
  maxForks 11 //DELETE WHEN DONE TESTING
  input:
    file chunk from deeplocChunks
  output:
    file "${chunk}.out.txt" into deeplocResults
  tag { assemblyPrefix+"-"+chunk }
  script:
    """
    export MKL_THREADING_LAYER=GNU
    export PATH="/lab/solexa_weng/testtube/miniconda3/bin:$PATH"
    ##deeploc -f ${chunk} -o ${chunk}.out
    touch ${chunk}.out.txt
    """
}
deeplocResults.collectFile(name: 'deeploc_annotations.txt').set { deeplocResult }

process tmhmmParallel {
  cpus 1
  input:
    file chunk from tmhmmChunks
  output:
    file "tmhmm_out" into tmhmmResults
  tag { assemblyPrefix+"-"+chunk }
  script:
    """
    echo tmhmm ${chunk}
    tmhmm --short < ${chunk} > tmhmm_out
    """
}
tmhmmResults.collectFile(name: 'tmhmm_annotations.tsv').set { tmhmmResult }

// Collect parallelized annotations
process annotatedFasta {
  publishDir "transXpress_results", mode: "copy"
  input:
    file transcriptomeFile from transcriptomeAnnotation
    file proteomeFile from predictProteome
    file kallistoFile from transcriptExpression
    file blastxResult 
    file blastpResult
    file pfamResult 
    file pfamDomResult 
    file deeplocResult
    file tmhmmResult 
  output:
    file assemblyPrefix+"_annotated.fasta" into transcriptome_annotated_fasta_ch
    file assemblyPrefix+"_annotated.pep" into transcriptome_annotated_pep_ch
    file "transcriptome_TPM_blast.csv"
    file "${blastxResult}"
    file "${blastpResult}"
    file "${pfamResult}"
    file "${pfamDomResult}"
    file "${deeplocResult}"
    file "${tmhmmResult}"
  tag { assemblyPrefix }
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
    print ("Loading tmhmm predictions from ${tmhmmResult}")
    with open("${tmhmmResult}") as input_handle:
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
        deeploc_annotations[row[0]] = str(row[1])
    
    ## Do the work
    print ("Annotating FASTA file ${transcriptomeFile}")
    with open("${transcriptomeFile}", 'r') as input_fasta_handle, open("${assemblyPrefix}_annotated.fasta", 'w') as output_fasta_handle:
      for record in Bio.SeqIO.parse(input_fasta_handle, "fasta"):
        transcript_id = record.id
        record.description = "TPM: " + expression_annotations.get(transcript_id)
        if transcript_id in blastx_annotations:
          record.description += "; blastx: " + blastx_annotations.get(transcript_id)
        Bio.SeqIO.write(record, output_fasta_handle, "fasta")
    
    print ("Annotating FASTA file ${proteomeFile}")
    with open("${proteomeFile}", 'r') as input_fasta_handle, open("${assemblyPrefix}_annotated.pep", 'w') as output_fasta_handle:
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
    publishDir "transXpress_results", mode: "copy"
    input:
        file "transcriptome_annotated.fasta" from transcriptome_annotated_fasta_ch1
        file "transcriptome_annotated.pep" from transcriptome_annotated_pep_ch1
    output:
        file "assembly_seq-dependent_checksums.txt"
    tag { assemblyPrefix }
    script:
    """
    echo -n "transcriptome_annotated.fasta:" > assembly_seq-dependent_checksums.txt
    seqkit sort -s transcriptome_annotated.fasta | grep -v ">" | md5sum | cut -f 1 -d ' '>> assembly_seq-dependent_checksums.txt
    echo -n "transcriptome_annotated.pep:" >> assembly_seq-dependent_checksums.txt
    seqkit sort -s transcriptome_annotated.pep | grep -v ">" | md5sum | cut -f 1 -d ' '>> assembly_seq-dependent_checksums.txt
    """    
}

