#!/usr/bin/env nextflow

/*
 * Notes: 
 * - must avoid spaces in filenames in samples.txt
 * - samples.txt must contain complete (absolute) paths to read files
 */


params.samples = "samples.txt"
samplesFile = file(params.samples)


log.info """
 transXpress
 ===================================
 Input samples file: ${params.samples}
 """


TRINITY_PARAMS = " --seqType fq"
TRINITY_PARAMS += " --trimmomatic --quality_trimming_params \"ILLUMINACLIP:${workflow.projectDir}/adapters.fasta:3:30:10 SLIDINGWINDOW:4:20 LEADING:20 TRAILING:20 MINLEN:25\""

// Used for both trinity and kallisto jobs!
STRAND_SPECIFIC = "" // --SS_lib_type=RF
TRINITY_PARAMS += " --min_glue 2"
TRINITY_PARAMS += " --min_kmer_cov 10"
TRINITY_PARAMS += " --no_normalize_reads"

SIGNALP_ORGANISM = "euk"


/*
 * Step 1. 
 */
process trinityInchwormChrysalis {
  input:
    file samplesFile
  output:
    file "trinity_out_dir" into trinityWorkDir
    file "trinity_out_dir/recursive_trinity.cmds" into trinityCmds
  cpus 12
  memory "200 GB"
  script:
    """
    Trinity --no_distributed_trinity_exec --max_memory ${task.memory.toGiga()}G --CPU ${task.cpus} --samples_file ${samplesFile} ${TRINITY_PARAMS}
    """
}

trinityCmds.splitText(by: 10, file: "trinityCmd").set { trinityParallelCmds }

process trinityButterflyParallel {
  input:
  //  file trinityWorkDir
    file parallelCommand from trinityParallelCmds
  output:
    file "${parallelCommand}.completed" into trinityFinishedCmds
  tag { parallelCommand }
  script:
    """
    sh ${parallelCommand}
    cp ${parallelCommand} ${parallelCommand}.completed
    """ 
}

trinityFinishedCmds.collectFile(name: "recursive_trinity.cmds.completed").set { trinityAllFinishedCmds }

process trinityFinish {
  publishDir ".", mode: "copy", saveAs: { filename -> filename.replaceAll("trinity_out_dir/Trinity", "transcriptome") }
  input:
    file samplesFile
    file trinityWorkDir
    file finishedCommands from trinityAllFinishedCmds
  output:
    file "trinity_out_dir/Trinity.fasta" into transcriptomeKallisto, transcriptomeTransdecoder, transcriptomeTransdecoderPredict, transcriptomeStats, transcriptomeTransrate, transcriptomeSplit, transcriptomeAnnotation
    file "trinity_out_dir/Trinity.fasta.gene_trans_map" into geneTransMap
  memory "1 GB"
  script:
    """
    cp ${finishedCommands} ${trinityWorkDir}/recursive_trinity.cmds.completed
    Trinity --samples_file ${samplesFile} --max_memory ${task.memory.toGiga()}G ${TRINITY_PARAMS}
    """ 
}

process transdecoderLongOrfs {
  publishDir ".", mode: "copy"
  input:
    file transcriptomeTransdecoder
  output:
    file "${transcriptomeTransdecoder}.transdecoder_dir/longest_orfs.pep" into proteomeSplit
    file "${transcriptomeTransdecoder}.transdecoder_dir" into transdecoderWorkDir
  script:
    """
    TransDecoder.LongOrfs -t ${transcriptomeTransdecoder}
    """
}

process transdecoderPredict {
  publishDir ".", mode: "copy" // , saveAs: { filename -> "transcriptome_after_predict.pep" }
  input:
    file transdecoderWorkDir
    file transcriptomeTransdecoderPredict
    file blastpForTransdecoder
    file pfamForTransdecoder
  output:
    file "${transcriptomeTransdecoderPredict}.transdecoder.pep" into proteomeAnnotation
    file "${transcriptomeTransdecoderPredict}.transdecoder.*"
  script:
    """
    TransDecoder.Predict -t ${transcriptomeTransdecoderPredict} --retain_pfam_hits ${pfamForTransdecoder} --retain_blastp_hits ${blastpForTransdecoder}
    """
}


process kallisto {
  publishDir ".", mode: "copy"
  cpus 10
  input:
    file transcriptomeKallisto
    file geneTransMap
    file samplesFile
  output:
    file "kallisto.isoform.TPM.not_cross_norm" into rawKallistoTable
    file "kallisto.isoform.TMM.EXPR.matrix" optional true into normalizedKallistoTable
    file "kallisto.gene.TPM.not_cross_norm" 
    file "kallisto.gene.TMM.EXPR.matrix" optional true
  script:
    """
    export TRINITY_HOME=\$(dirname `which Trinity`)
    echo TRINITY_HOME set to \${TRINITY_HOME}
    \${TRINITY_HOME}/util/align_and_estimate_abundance.pl --transcripts ${transcriptomeKallisto} ${STRAND_SPECIFIC} --seqType fq --samples_file ${samplesFile} --prep_reference --thread_count ${task.cpus} --est_method kallisto --gene_trans_map ${geneTransMap}
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
  publishDir ".", mode: "copy"
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
  publishDir ".", mode: "copy", saveAs: { filename -> "transrate_results.csv" }
  cpus 8
  input:
    file samplesFile 
    file transcriptomeTransrate 
  output:
    file "transrate_results/Trinity/contigs.csv" into transrateResults
  script:
    """
    LEFT=`cut -f 3 < ${samplesFile} | tr '\n' ',' | sed 's/,*\$//g'`
    RIGHT=`cut -f 4 < ${samplesFile} | tr '\n' ',' | sed 's/,*\$//g'`
    transrate --threads ${task.cpus} --assembly=${transcriptomeTransrate} --left=\$LEFT --right=\$RIGHT
    """
}

transcriptomeSplit
  .splitFasta(by: 100, file: true)
  .set { sprotBlastxChunks }

proteomeSplit
  .splitFasta(by: 100, file: true)
  .into { sprotBlastpChunks; pfamChunks; signalpChunks; tmhmmChunks }


process downloadPfam {
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

process downloadSprot {
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
  tag { chunk }
  script:
    """
    echo blastp ${chunk} using database ${sprotDb}
    blastp -query ${chunk} -db ${sprotDb} -num_threads ${task.cpus} -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt "6 std stitle" -out blastp_out
    """
}

process pfamParallel {
  cpus 2
  input:
    file chunk from pfamChunks
    set pfamDb, pfamDbIndex from pfamDb
  output:
    file "pfam_out" into pfamResults
    file "pfam_dom_out" into pfamDomResults
  tag { chunk }
  script:
    """
    echo pfam ${chunk} using database ${pfamDb}
    hmmscan --cpu ${task.cpus} --domtblout pfam_dom_out --tblout pfam_out ${pfamDb} ${chunk}
    """
}

process signalpParallel {
  cpus 1
  input:
    file chunk from signalpChunks
  output:
    file "signalp_out" into signalpResults
  tag { chunk }
  script:
    """
    echo signalp ${chunk}
    signalp -t ${SIGNALP_ORGANISM} -f short ${chunk} > signalp_out
    """
}

process tmhmmParallel {
  cpus 1
  input:
    file chunk from tmhmmChunks
  output:
    file "tmhmm_out" into tmhmmResults
  tag { chunk }
  script:
    """
    echo tmhmm ${chunk}
    tmhmm --short < ${chunk} > tmhmm_out
    """
}

// Collect parallelized annotations
sprotBlastxResults.collectFile(name: 'blastx_annotations.tsv').set { blastxResult }
sprotBlastpResults.collectFile(name: 'blastp_annotations.tsv').into { blastpForTransdecoder; blastpResult }
pfamResults.collectFile(name: 'pfam_annotations.txt').set { pfamResult }
pfamDomResults.collectFile(name: 'pfam_dom_annotations.txt').into { pfamDomResult ; pfamForTransdecoder }
signalpResults.collectFile(name: 'signalp_annotations.txt').set { signalpResult }
tmhmmResults.collectFile(name: 'tmhmm_annotations.tsv').set { tmhmmResult }


process annotatedFasta {
  publishDir ".", mode: "copy"
  input:
    file transcriptomeFile from transcriptomeAnnotation
    file proteomeFile from proteomeAnnotation
    file transrateFile from transrateResults
    file kallistoFile from transcriptExpression
    file blastxResult 
    file blastpResult
    file pfamResult 
    file pfamDomResult 
    file signalpResult
    file tmhmmResult 
  output:
    file "transcriptome_annotated.fasta"
    file "transcriptome_annotated.pep"
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
    with open("${transcriptomeFile}", 'r') as input_fasta_handle, open("transcriptome_annotated.fasta", 'w') as output_fasta_handle:
      for record in Bio.SeqIO.parse(input_fasta_handle, "fasta"):
        transcript_id = record.id
        record.description = "TPM: " + expression_annotations.get(transcript_id)
        if transcript_id in transrate_annotations:
          record.description += "; transrate: " + transrate_annotations.get(transcript_id)
        if transcript_id in blastx_annotations:
          record.description += "; blastx: " + blastx_annotations.get(transcript_id)
        Bio.SeqIO.write(record, output_fasta_handle, "fasta")
    
    print ("Annotating FASTA file ${proteomeFile}")
    with open("${proteomeFile}", 'r') as input_fasta_handle, open("transcriptome_annotated.pep", 'w') as output_fasta_handle:
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



