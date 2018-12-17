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

params.SIGNALP_ORGANISMS = "euk"

/*
 * Step 1. 
 */


///
/// Load read files
///
Channel.fromPath(params.samples)
     .splitCsv(sep:'\t',header:false)
     .map{ row ->
     println row 
     return tuple(file(row[2]), file(row[3])) }
     .set{ readPairs_ch }
 
process trimmomatic {
cpus 4
tag {"$R1_reads"+" and " +"$R2_reads"}
input:
set file(R1_reads),file(R2_reads) from readPairs_ch

output:
  file "${R1_reads}.R1-P.qtrim.fastq.gz" into filteredForwardReads_ch
  file "${R2_reads}.R2-P.qtrim.fastq.gz" into filteredReverseReads_ch
  file "*U.qtrim.fastq.gz" into filteredSingleReads_ch1,filteredSingleReads_ch2
script:
"""
java -jar /lab/solexa_weng/testtube/trinityrnaseq-Trinity-v2.8.4/trinity-plugins/Trimmomatic/trimmomatic.jar PE -threads ${task.cpus}  ${R1_reads} ${R2_reads} ${R1_reads}.R1-P.qtrim.fastq.gz ${R1_reads}.R1-U.qtrim.fastq.gz ${R2_reads}.R2-P.qtrim.fastq.gz ${R2_reads}.R2-U.qtrim.fastq.gz  ILLUMINACLIP:/lab/solexa_weng/testtube/trinityrnaseq-Trinity-v2.8.4/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25 
"""
}
filteredForwardReads_ch.into{ filteredForwardReads_ch1; filteredForwardReads_ch2; filteredForwardReads_ch3; filteredForwardReads_ch4; filteredForwardReads_ch5 }
filteredReverseReads_ch.into{ filteredReverseReads_ch1; filteredReverseReads_ch2; filteredReverseReads_ch3; filteredReverseReads_ch4; filteredReverseReads_ch5 }


process convertSamplesToRelative {
input:
    file "samples.txt" from file(params.samples)
output:
    file "relative_samples.txt" into relative_samples_txt_ch

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
relative_samples_txt_ch.into{ relative_samples_txt_ch1; relative_samples_txt_ch2; relative_samples_txt_ch3; relative_samples_txt_ch4}

process convertReadsToYAML {

input:
  file forwardReads from filteredForwardReads_ch1.collect()
  file reverseReads from filteredReverseReads_ch1.collect()
  file unpairedReads from filteredSingleReads_ch1.collect()
output:
  file "datasets.yaml" into datasets_YAML_ch
script:
"""
   echo "[{
        orientation: "fr",
        type: "paired-end",
        right reads: [" >datasets.yaml
   ls -1 ./*.R1-P.qtrim.fastq.gz | sed 's/^/\\"/g' | sed 's/\$/\\"/g' | sed 's/@\"/\\"/g' | tr "\\n" "," | sed 's/,\$//g' >>datasets.yaml
   echo "]," >> datasets.yaml
   echo "left reads: [" >> datasets.yaml
   ls -1 ./*.R2-P.qtrim.fastq.gz | sed 's/^/\\"/g' | sed 's/\$/\\"/g' | sed 's/@\"/\\"/g' | tr "\\n" "," | sed 's/,\$//g' >>datasets.yaml
   echo "]
      }," >>datasets.yaml
   echo "{
        type: "single",
        single reads: [" >> datasets.yaml
   ls -1 ./*U.qtrim.fastq.gz | sed 's/^/\\"/g' | sed 's/\$/\\"/g' | sed 's/@\"/\\"/g' | tr "\\n" "," | sed 's/,\$//g' >>datasets.yaml
   echo "]}]" >> datasets.yaml
"""  
}

process trinityInchwormChrysalis {
  echo = true
  label = "nf_"+assemblyPrefix+"_trinityInchwormChrysalis"
  stageInMode="copy" 

  cpus 12
  memory "200 GB"

  tag { assemblyPrefix }

  afterScript 'echo \"(Above completion message is from Trinity. transXpress will continue the pipeline execution.)\"'
  //afterScript 'exit(1)'
  input:
    file filteredForwardReads from filteredForwardReads_ch2.collect()
    file filteredReverseReads from filteredReverseReads_ch2.collect()
    file relative_samples_txt from relative_samples_txt_ch1
  output:
    //These commented out lines load all the files into channels.
    //file "./trinity_out_dir/*.fa" into trinityPhase1RootFiles_ch1
    //file "./trinity_out_dir/*.ok" into trinityPhase1RootFiles_ch2
    //file "./trinity_out_dir/*.read_count" into trinityPhase1RootFiles_ch3
    //file "./trinity_out_dir/*.kmer_count" into trinityPhase1RootFiles_ch4
    //file "./trinity_out_dir/*.histo" into trinityPhase1RootFiles_ch5
    //file "./trinity_out_dir/*.cmds" into trinityPhase1RootFiles_ch6
    //file "./trinity_out_dir/*.timing" into trinityPhase1RootFiles_ch7
    //file "./trinity_out_dir/*.sam" into trinityPhase1RootFiles_ch8
    //file "./trinity_out_dir/chrysalis/*.fasta" into trinityPhase1ChrysalisFiles_ch1
    //file "./trinity_out_dir/chrysalis/*.ok" into trinityPhase1ChrysalisFiles_ch2
    //file "./trinity_out_dir/chrysalis/*.out" into trinityPhase1ChrysalisFiles_ch3
    //file "./trinity_out_dir/chrysalis/*.min100" into trinityPhase1ChrysalisFiles_ch4
    //file "./trinity_out_dir/chrysalis/*.bt2" into trinityPhase1ChrysalisFiles_ch5
    //file "./trinity_out_dir/chrysalis/*.bam" into trinityPhase1ChrysalisFiles_ch6
    //file "./trinity_out_dir/chrysalis/*.sorted" into trinityPhase1ChrysalisFiles_ch7
    //file "./trinity_out_dir/chrysalis/*.txt" into trinityPhase1ChrysalisFiles_ch8
    //file "./trinity_out_dir/chrysalis/*.wIwormNames" into trinityPhase1ChrysalisFiles_ch9
    //file "./trinity_out_dir/chrysalis/*.sort" into trinityPhase1ChrysalisFiles_ch10
    //file "./trinity_out_dir/read_partitions/*/*/*.trinity.reads.fa" trinityPhase1ReadPartitionsFiles_ch
    file "./trinity_out_dir" into trinityWorkDir
    file "./trinity_out_dir/recursive_trinity.cmds" into trinityCmds

  script:
    """
    Trinity --no_distributed_trinity_exec --max_memory ${task.memory.toGiga()}G --CPU ${task.cpus} --samples_file ${relative_samples_txt} ${params.TRINITY_PARAMS}
    """
}

process trinityButterflyParallel {
  input:
  //  file trinityWorkDir
    file parallelCommand from trinityCmds.splitText(by: 10, file: "trinityCmd")
  output:
    file "${parallelCommand}.completed" into trinityFinishedCmds
  tag { assemblyPrefix+"-"+parallelCommand }
  script:
    """
    sh ${parallelCommand}
    cp ${parallelCommand} ${parallelCommand}.completed
    """ 
}

process trinityFinish {
   publishDir "transXpress_results", mode: "copy", saveAs: { filename -> filename.replaceAll("trinity_out_dir/Trinity", "transcriptome") }
  input:
    //file "./trinity_out_dir/*" from trinityPhase1RootFiles_ch1.mix(trinityPhase1RootFiles_ch2,trinityPhase1RootFiles_ch3,trinityPhase1RootFiles_ch4,trinityPhase1RootFiles_ch5,trinityPhase1RootFiles_ch6,trinityPhase1RootFiles_ch7,trinityPhase1RootFiles_ch8).collect()
    //file "./trinity_out_dir/chrysalis/*" from trinityPhase1ChrysalisFiles_ch1.mix(trinityPhase1ChrysalisFiles_ch2,trinityPhase1ChrysalisFiles_ch3,trinityPhase1ChrysalisFiles_ch4,trinityPhase1ChrysalisFiles_ch5,trinityPhase1ChrysalisFiles_ch6,trinityPhase1ChrysalisFiles_ch7,trinityPhase1ChrysalisFiles_ch8,trinityPhase1ChrysalisFiles_ch9,trinityPhase1ChrysalisFiles_ch10).collect()
    //file "trinity_out_dir/read_partitions/*/*/*.trinity.reads.fa" from trinityPhase1ChrysalisFiles_ch
    file trinityWorkDir
    file relative_samples_txt from relative_samples_txt_ch2
    file finishedCommands from trinityFinishedCmds.collectFile(name: "recursive_trinity.cmds.completed")
    file filteredForwardReads from filteredForwardReads_ch3.collect()
    file filteredReverseReads from filteredReverseReads_ch3.collect()
  output:
    file "./trinity_out_dir/Trinity.fasta.gene_trans_map" into originalGeneTransMap
    file "./trinity_out_dir/Trinity.fasta" into Trinity_fasta_ch
    file "./trinity_out_dir/recursive_trinity.cmds.completed"
  memory "1 GB"
  script:
    """
    cp ${finishedCommands} ${trinityWorkDir}/recursive_trinity.cmds.completed
    Trinity --samples_file ${relative_samples_txt} --max_memory ${task.memory.toGiga()}G ${params.TRINITY_PARAMS}
    """ 
}

process renameTrinityAssembly {
   publishDir "transXpress_results", mode: "copy"

   input: 
    file "Trinity.fasta" from Trinity_fasta_ch
    file "species.txt" from file(params.species) //Just a dummy input
    file "Trinity.fasta.gene_trans_map" from originalGeneTransMap
   output:
    file assemblyPrefix+".fasta" into transcriptomeKallisto, transcriptomeTransdecoder, transcriptomeTransdecoderPredict, transcriptomeStats, transcriptomeTransrate, transcriptomeSplit, transcriptomeAnnotation
    file "Trinity_renamed.fasta.gene_trans_map" into geneTransMap

   script:
   """
   seqkit replace -p '^TRINITY' -r '${assemblyPrefix}' Trinity.fasta > ${assemblyPrefix}.fasta 
   cat Trinity.fasta.gene_trans_map | sed 's/^TRINITY/${assemblyPrefix}/' | sed 's/\tTRINITY/\t${assemblyPrefix}/g' > Trinity_renamed.fasta.gene_trans_map
   """
}

process transdecoderLongOrfs {
  publishDir "transXpress_results", mode: "copy"
  input:
    file transcriptomeTransdecoder
  output:
    file "${transcriptomeTransdecoder}.transdecoder_dir/longest_orfs.pep" into longOrfsProteomeSplit
    file "${transcriptomeTransdecoder}.transdecoder_dir" into transdecoderWorkDir
  script:
    """
    TransDecoder.LongOrfs -t ${transcriptomeTransdecoder}
    """
}

process kallisto {
  publishDir "transXpress_results", mode: "copy"
  cpus 10
  input:
    file forwardReads from filteredForwardReads_ch5.collect()
    file reverseReads from filteredReverseReads_ch5.collect()
    file transcriptomeKallisto
    file geneTransMap
    file relative_samples_txt from relative_samples_txt_ch4
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

process transrate {
  publishDir "transXpress_results", mode: "copy", saveAs: { filename -> "transrate_results.csv" }
  //Can't be implemented until the absolute vs relative path in samples.txt is figured out
  //Would be nice if there were an upstream module that pre-filtered all the read naming things.
  //scratch params.scratch_dir
  cpus 8
  input:
    file filteredForwardReads from filteredForwardReads_ch4.collect()
    file filteredReverseReads from filteredReverseReads_ch4.collect()
    file relative_samples_txt from relative_samples_txt_ch3
    file transcriptomeTransrate 
  output:
    file "transrate_results/"+assemblyPrefix+"/contigs.csv" into transrateResults

/// This strips out the odd read naming from SRA, that screws up transrate
/// seqkit replace -p '_forward/1' -r '' \$LEFT | gzip > F_reads.fq.gz
/// seqkit replace -p '_reverse/2' -r '' \$RIGHT | gzip > R_reads.fq.gz

  script:
    """
    LEFT=`cut -f 3 < ${relative_samples_txt} | tr '\n' ' ' | sed 's/ *\$//g'`
    RIGHT=`cut -f 4 < ${relative_samples_txt} | tr '\n' ' ' | sed 's/ *\$//g'`
    seqkit replace -p '_forward/1' -r '' \$LEFT | gzip > F_reads.fq.gz
    seqkit replace -p '_reverse/2' -r '' \$RIGHT | gzip > R_reads.fq.gz
    transrate --threads ${task.cpus} --assembly=${transcriptomeTransrate} --left=F_reads.fq.gz --right=R_reads.fq.gz
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
  storeDir 'db'
  output:
    set file("Pfam-A.hmm"), file("Pfam-A.hmm.h??") into pfamDb
  script:
    """
    wget "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
    gunzip Pfam-A.hmm.gz
    hmmpress Pfam-A.hmm
    """
}

process downloadRfam {
  executor 'local'
  storeDir 'db'
  output:
    set file("Rfam_with_desc.cm"), file("Rfam_with_desc.cm.???") into rfamDb
  script:
    """
    wget "ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam_with_desc.cm.gz"
    gunzip Rfam_with_desc.cm.gz
    cmpress Rfam_with_desc.cm
    """
}

process downloadSprot {
  executor 'local'
  storeDir 'db'
  output:
    set file("uniprot_sprot.fasta"), file("uniprot_sprot.fasta.p??") into sprotDb
  script:
    """
    wget "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
    gunzip uniprot_sprot.fasta.gz
    makeblastdb -in uniprot_sprot.fasta -dbtype prot
    """
}

/*
This takes a looong time, even when it is working
process downloadVirusesUniref50 {
  executor 'local'
  storeDir 'db'
  errorStrategy 'ignore'
  output:
    set file("virusesUniref50.fasta"), file("virusesUniref50.pep.fasta.p??") into virusDb
  script:
    """
    wget -t 3 -O virusesUniref50.pep.fasta.gz "https://www.uniprot.org/uniref/?query=uniprot%3A%28taxonomy%3A%22Viruses+%5B10239%5D%22%29+AND+identity%3A0.5&format=fasta&compress=yes"
    makeblastdb -in virusesUniref50.pep.fasta -dbtype prot
    """
}
*/

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
    file transdecoderWorkDir
    file transcriptomeTransdecoderPredict
    file blastpForTransdecoder
    file pfamForTransdecoder
  output:
    file "${transcriptomeTransdecoderPredict}.transdecoder.pep" into predictProteome, predictProteomeSplit //This seems a bit weird. Referring to it indirectly, rather than directly
    file "${transcriptomeTransdecoderPredict}.transdecoder.*"
  script:
    """
    TransDecoder.Predict -t ${transcriptomeTransdecoderPredict} --retain_pfam_hits ${pfamForTransdecoder} --retain_blastp_hits ${blastpForTransdecoder}
    """
}

predictProteomeSplit
  .splitFasta(by: 100, file: true)
  .into { signalpChunks; tmhmmChunks }

process signalpParallel {
  cpus 1
  input:
    file chunk from signalpChunks
  output:
    file "signalp_out" into signalpResults
  tag { assemblyPrefix+"-"+chunk }
  script:
    """
    echo signalp ${chunk}
    signalp -t ${params.SIGNALP_ORGANISMS} -f short ${chunk} > signalp_out
    """
}
signalpResults.collectFile(name: 'signalp_annotations.txt').set { signalpResult }

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
    file transrateFile from transrateResults
    file kallistoFile from transcriptExpression
    file blastxResult 
    file blastpResult
    file pfamResult 
    file pfamDomResult 
    file signalpResult
    file tmhmmResult 
  output:
    file assemblyPrefix+"_annotated.fasta" into transcriptome_annotated_fasta_ch
    file assemblyPrefix+"_annotated.pep" into transcriptome_annotated_pep_ch
    file "transcriptome_TPM_blast.csv"
    file "${blastxResult}"
    file "${blastpResult}"
    file "${pfamResult}"
    file "${pfamDomResult}"
    file "${signalpResult}"
    file "${tmhmmResult}"
  script:
    """
    #!/usr/bin/env python
    
    import re
    import csv
    import Bio.SeqIO
  
    ## Annotation maps: transcript id -> annotation
    expression_annotations = {}
    transrate_annotations = {}
    blastx_annotations = {}
    blastp_annotations = {}
    pfam_annotations = {}
    tmhmm_annotations = {}
    signalp_annotations = {}

    ## Load kallisto results
    print ("Loading expression values from ${kallistoFile}")
    with open("${kallistoFile}") as input_handle:
      csv_reader = csv.reader(input_handle, delimiter='\t')
      columns = next(csv_reader)
      for row in csv_reader:
        expression_annotations[row[0]] = columns[1] + "=" + str(row[1])
        for i in range(2, len(columns)):
          expression_annotations[row[0]] += " " + columns[i] + "=" + str(row[i])

    ## Load transrate results
    print ("Loading transrate results from ${transrateFile}")
    with open("${transrateFile}") as input_handle:
      csv_reader = csv.reader(input_handle, delimiter=',')
      columns = next(csv_reader)
      for row in csv_reader:
        if (len(row) < 18): continue
        transrate_annotations[row[0]] = columns[5] + "=" + str(row[5]) + " " + columns[7] + "=" + str(row[7])

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
    
    ## Load signalp results
    print ("Loading signalp predictions from ${signalpResult}")
    with open("${signalpResult}") as input_handle:
      for line in input_handle:
        if (line.startswith("#")): continue
        row = re.split(" +", line)
        if (len(row) < 9): continue
        signalp_annotations[row[0]] = str(row[9]) + ", D=" + str(row[8]) + ", pos=" + str(row[2])
    
    ## Do the work
    print ("Annotating FASTA file ${transcriptomeFile}")
    with open("${transcriptomeFile}", 'r') as input_fasta_handle, open("${assemblyPrefix}_annotated.fasta", 'w') as output_fasta_handle:
      for record in Bio.SeqIO.parse(input_fasta_handle, "fasta"):
        transcript_id = record.id
        record.description = "TPM: " + expression_annotations.get(transcript_id)
        if transcript_id in transrate_annotations:
          record.description += "; transrate: " + transrate_annotations.get(transcript_id)
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
        if transcript_id in transrate_annotations:
          record.description += "; transrate: " + transrate_annotations.get(transcript_id)
        if record.id in blastp_annotations:
          record.description += "; blastp: " + blastp_annotations.get(record.id)
        if record.id in pfam_annotations:
          record.description += "; pfam: " + pfam_annotations.get(record.id)
        if record.id in tmhmm_annotations:
          record.description += "; tmhmm: " + tmhmm_annotations.get(record.id)
        if record.id in signalp_annotations:
          record.description += "; signalp: " + signalp_annotations.get(record.id)
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
      csv_columns.append("transrate")
      csv_writer.writerow(csv_columns)
      for row in csv_reader:
        row.append(blastx_annotations.get(row[0], ""))
        row.append(transrate_annotations.get(row[0], ""))
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
    script:
    """
    echo -n "transcriptome_annotated.fasta:" > assembly_seq-dependent_checksums.txt
    seqkit sort -s transcriptome_annotated.fasta | grep -v ">" | md5sum | cut -f 1 -d ' '>> assembly_seq-dependent_checksums.txt
    echo -n "transcriptome_annotated.pep:" >> assembly_seq-dependent_checksums.txt
    seqkit sort -s transcriptome_annotated.pep | grep -v ">" | md5sum | cut -f 1 -d ' '>> assembly_seq-dependent_checksums.txt
    """    
}


Channel.fromPath(["/lab/solexa_weng/Seq_data/Projects/Tim_Fallon/BUSCO_profiles/BUSCO_v2_profiles/bacteria_odb9/",\
"/lab/solexa_weng/Seq_data/Projects/Tim_Fallon/BUSCO_profiles/BUSCO_v2_profiles/fungi_odb9/",\
"/lab/solexa_weng/Seq_data/Projects/Tim_Fallon/BUSCO_profiles/BUSCO_v2_profiles/metazoa_odb9/",\
"/lab/solexa_weng/Seq_data/Projects/Tim_Fallon/BUSCO_profiles/BUSCO_v2_profiles/eukaryota_odb9/"]).set{ BUSCO_lineages }

BUSCO_lineages.into { BUSCO_lineages_ch1 ; BUSCO_lineages_ch2 }

BUSCO_cmds_pep = BUSCO_lineages_ch1.combine(transcriptome_annotated_pep_ch2)
BUSCO_cmds_trans = BUSCO_lineages_ch2.combine(transcriptome_annotated_fasta_ch2)

BUSCO_cmds_mixed = BUSCO_cmds_pep.mix(BUSCO_cmds_trans)

process do_BUSCO {
 publishDir "transXpress_results", mode: "copy"
 cpus 21
 //scratch params.scratch_dir
 input:
     set file(BUSCO_lineage), file(inputFasta) from BUSCO_cmds_mixed
 output:
     file 'run_'+assemblyPrefix+'*/*'
 tag { inputFasta+"_"+BUSCO_lineage }
 """
 #! /bin/bash
 ##A new script that uses BUSCO_v3

 ##Note that BUSCO_v3 
 #In short, I’ve cloned from the BUSCO_v3 github (https://gitlab.com/ezlab/busco) to -> /lab/solexa_weng/testtube/busco
 #I’ve installed it to my user directory (python setup.py install --user)
 #I’ve updated the paths & other config things in (/lab/solexa_weng/testtube/busco/config/config.ini), to the Tak paths for the appropriate software, & other settings

 export AUGUSTUS_CONFIG_PATH=/lab/solexa_weng/testtube/busco/augustus_config/ ##I copied the augustus config directory myself.

 BUSCO_SCRIPTS=/lab/solexa_weng/testtube/busco/scripts/

 NAME=\$(basename ${inputFasta})
 LINEAGE_NAME=\$(basename ${BUSCO_lineage})
 
 if [[ \$NAME == *'.pep'* ]]; then
  TYPE=prot
 else
  TYPE=transcriptome
 fi

 \${BUSCO_SCRIPTS}/run_BUSCO.py -z -t /mnt/ramdisk -f -i ${inputFasta} -l ${BUSCO_lineage} -o \${NAME}_\${LINEAGE_NAME} -m \$TYPE --cpu ${task.cpus} >./BUSCO.\${NAME}_\${LINEAGE_NAME}.stdout.log 2>./BUSCO.\${NAME}_\${LINEAGE_NAME}.stderr.log
 find /mnt/ramdisk -name \"*\${NAME}*\${LINEAGE_NAME}*\" | xargs rm -f
 """
 
}

