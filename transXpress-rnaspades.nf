#!/usr/bin/env nextflow

/*
 * Notes: 
 * - must avoid spaces in filenames in samples.txt
 * - samples.txt must contain complete (absolute) paths to read files
 */


params.TRINITY_PARAMS += " --trimmomatic --quality_trimming_params \"ILLUMINACLIP:${workflow.projectDir}/adapters.fasta:3:30:10 SLIDINGWINDOW:4:20 LEADING:20 TRAILING:20 MINLEN:25\""
params.tempdir = "/lab/weng_scratch/tmp/"

theDate = new java.util.Date().format( 'MMdd' )

theText = file(params.species).text
genus = theText.split(" ")[0]
species = theText.split(" ")[1].trim()
assemblyPrefix = theDate+"_"+genus+"_"+species+"_rnaSPAdes"

log.info """
 transXpress
 ===================================
 """+assemblyPrefix

params.SIGNALP_ORGANISMS = "euk"

/*
 * Step 1. 
 */

Channel.fromPath(params.samples)
     .splitCsv(sep:'\t',header:false)
     .map{ row ->
     println row 
     return tuple(file(row[2]), file(row[3])) }
     .into{ readPairs_ch ; readPairs_ch2 }

readPairs_ch2
       .map{ item ->
       item0 = '\"'+item[0]+'\"' //Double quoting the string
       item1 = '\"'+item[1]+'\"' //Double quoting the string
       return tuple(item0,item1) }
       .set{ readPairsQuoted_ch1 }

process trimmomatic {

cpus 4
input:
 set file(R1_reads),file(R2_reads) from readPairs_ch
tag {"$R1_reads"+" and " +"$R2_reads"}
output:
  set file("${R1_reads}.R1-P.qtrim.fastq.gz"), file("${R2_reads}.R2-P.qtrim.fastq.gz") into filteredPairedReads_ch1,filteredPairedReads_ch2,filteredPairedReads_ch5
  file "*U.qtrim.fastq.gz" into filteredSingleReads_ch1,filteredSingleReads_ch2
script:
"""
java -jar /lab/solexa_weng/testtube/trinityrnaseq-Trinity-v2.8.4/trinity-plugins/Trimmomatic/trimmomatic.jar PE -threads ${task.cpus}  ${R1_reads} ${R2_reads} ${R1_reads}.R1-P.qtrim.fastq.gz ${R1_reads}.R1-U.qtrim.fastq.gz ${R2_reads}.R2-P.qtrim.fastq.gz ${R2_reads}.R2-U.qtrim.fastq.gz  ILLUMINACLIP:/lab/solexa_weng/testtube/trinityrnaseq-Trinity-v2.8.4/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25 
"""
}


process convertReadsToYAML {
input:
  val readPairTuples from readPairsQuoted_ch1.toList() //toList doesn't flatten the tuple structure
  file pairedReads from filteredPairedReads_ch1.collect() //collect does flatten the tuple structure
  file unpairedReads from filteredSingleReads_ch1.collect()
output:
  file "datasets.yaml" into datasets_YAML_ch
script:
    """
    #!/usr/bin/env python
    import glob
    import pprint
    
    readPairs = ${readPairTuples} ##Little bit of hackery. Since groovy and python are similar, quoting the string values in the groovy datastructure makes a python datastructure
    
    sample_list = []
    for p in readPairs:
        f = p[0].split("/")[-1:][0] ##Just want the last part of the path.
        r = p[1].split("/")[-1:][0] ##Just want the last part of the path.
        assert len(glob.glob(f+"*.R1-P.qtrim.fastq.gz")) == 1
        assert len(glob.glob(r+"*.R2-P.qtrim.fastq.gz")) == 1
        f_trimmed = glob.glob(f+"*.R1-P.qtrim.fastq.gz") ##Use the last part of the path to find the trimmed version
        r_trimmed = glob.glob(r+"*.R2-P.qtrim.fastq.gz") ##Use the last part of the path to find the trimmed version
        sample_dict = dict()
        sample_dict['orientation'] = 'fr'
        sample_dict['type'] = 'paired-end'
        sample_dict['right reads'] = f_trimmed
        sample_dict['left reads'] = r_trimmed
        sample_list.append(sample_dict)

    unpaired_reads = glob.glob("*U.qtrim.fastq.gz")

    for u in unpaired_reads:
         sample_dict = dict()
         sample_dict['type'] = 'single'
         sample_dict['single reads'] = [u]
         sample_list.append(sample_dict)
    handle = open("datasets.yaml","w")
    handle.write(pprint.pformat(sample_list))
    handle.close()
    """  
}

process runSPAdes {
cpus 16
memory "200 GB"

input:
   file datasets_YAML from datasets_YAML_ch
   file pairedReads from filteredPairedReads_ch2.collect() //collect flattens the tuple structure
   //file filteredForwardReads from filteredForwardReads_ch2.collect()
   //file filteredReverseReads from filteredReverseReads_ch2.collect()
   file filteredSingleReads from filteredSingleReads_ch2.collect()
output:
   file assemblyPrefix+"/transcripts.fasta" into spadesAssembly_ch

script:
"""
rnaspades.py --dataset ${datasets_YAML} -t ${task.cpus} -m ${task.memory.toGiga()} --fast -o ${assemblyPrefix} --only-assembler -k 47
"""

}

process renameAssembly {
   publishDir "transXpress_results", mode: "copy"

   input: 
    file spadesAssembly_ch
   output:
    file assemblyPrefix+".fasta" into transcriptomeKallisto, transcriptomeTransdecoder, transcriptomeTransdecoderPredict, transcriptomeSplit, transcriptomeAnnotation

   script:
   """
   seqkit replace -p '^' -r '${assemblyPrefix}_' ${spadesAssembly_ch} > ${assemblyPrefix}.fasta 
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
  tag { assemblyPrefix+"-"+chunk }
  script:
    """
    echo rfam ${chunk} using database ${rfamDb}
    cmscan --incE 0.00001 --rfam --cpu ${task.cpus} --tblout rfam_out ${rfamDb} ${chunk}
    """
}

process publishRfamResults {
  publishDir "transXpress_results", mode: "copy"
  input:
    file rfamResult from rfamResults.collectFile(name: 'rfam_annotations_unsorted.txt')
  output:
    file "rfam_annotations.txt" into rfamResultPub1
  script:
  """
  cat ${rfamResult} |  grep -v "#" | sort -k3,3 -k15nr,15 | sort -u -k3,3 --merge | sort -k15nr,15 > rfam_annotations_sorted.txt
  cat ${rfamResult} | head -n 2 > header.txt
  cat ${rfamResult} | tail -n 9 > footer.txt
  cat header.txt rfam_annotations_sorted.txt footer.txt > rfam_annotations.txt
  """
}

process transdecoderPredict {
  publishDir "transXpress_results", mode: "copy" // , saveAs: { filename -> "transcriptome_after_predict.pep" }
  stageInMode="copy"
  input:
    file transdecoderWorkDir
    file transcriptomeTransdecoderPredict
    file blastpForTransdecoder
    file pfamForTransdecoder
  output:
    file "${transcriptomeTransdecoderPredict}.transdecoder.pep" into predictProteome, predictProteomeSplitBy100,predictProteomeSplitBy10 //This seems a bit weird. Referring to it indirectly, rather than directly
    file "${transcriptomeTransdecoderPredict}.transdecoder.*"
  script:
    """
    TransDecoder.Predict -t ${transcriptomeTransdecoderPredict} --retain_pfam_hits ${pfamForTransdecoder} --retain_blastp_hits ${blastpForTransdecoder}
    """
}

predictProteomeSplitBy100
  .splitFasta(by: 100, file: true)
  .into { signalpChunks; tmhmmChunks }

predictProteomeSplitBy10
  .splitFasta(by: 10, file: true)
  .set{deeplocChunks}

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


process deeplocParallel {
  maxForks 11 //DELETE_THIS_AFTER_TESTING
  input:
    file chunk from deeplocChunks
  output:
    file "${chunk}.out.txt" into deeplocResults
  tag { assemblyPrefix+"-"+chunk }
  script:
    """
    export MKL_THREADING_LAYER=GNU
    export PATH="/lab/solexa_weng/testtube/miniconda3/bin:$PATH"
    
    deeploc -f ${chunk} -o ${chunk}.out
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

process kallistoIndex {
 input:
  file transcriptFasta from transcriptomeKallisto
 output:
  file "*.kalidx" into kallistoIndex
 tag {assemblyPrefix}
 script:
 """
 kallisto index -i ${transcriptFasta}.kalidx ${transcriptFasta}
 """
}

process kallistoDirect {
 cpus 16
 input:
  file kallistoIndex
  set file(forwardReads), file(reverseReads) from filteredPairedReads_ch5
 output:
  file "output*" into kallistoResults
 tag {assemblyPrefix + forwardReads}
 script:
 """
 kallisto quant -t ${task.cpus} -i ${kallistoIndex} -o output -b 100 ${forwardReads} ${reverseReads}
 """
}


// Collect parallelized annotations
process annotatedFasta {
  publishDir "transXpress_results", mode: "copy"
  input:
    file transcriptomeFile from transcriptomeAnnotation
    file proteomeFile from predictProteome
    //file kallistoFile from transcriptExpression
    file blastxResult 
    file blastpResult
    file pfamResult 
    file pfamDomResult 
    file signalpResult
    file tmhmmResult
    file deeplocResult 
  output:
    file assemblyPrefix+"_annotated.fasta" into transcriptome_annotated_fasta_ch
    file assemblyPrefix+"_annotated.pep" into transcriptome_annotated_pep_ch
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
    deeploc_annotations = {}

    ## Load kallisto results
    ##Unplugged transtrate, add a dollar sign in front of kallistoFile to plug back in
    ##print ("Loading expression values from !{kallistoFile}")
    ##with open("!{kallistoFile}") as input_handle:
    ##  csv_reader = csv.reader(input_handle, delimiter='\t')
    ##  columns = next(csv_reader)
    ##  for row in csv_reader:
    ##    expression_annotations[row[0]] = columns[1] + "=" + str(row[1])
    ##    for i in range(2, len(columns)):
    ##      expression_annotations[row[0]] += " " + columns[i] + "=" + str(row[i])

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
        if transcript_id in blastx_annotations:
          record.description += "; blastx: " + blastx_annotations.get(transcript_id)
        Bio.SeqIO.write(record, output_fasta_handle, "fasta")
    
    print ("Annotating FASTA file ${proteomeFile}")
    with open("${proteomeFile}", 'r') as input_fasta_handle, open("${assemblyPrefix}_annotated.pep", 'w') as output_fasta_handle:
      for record in Bio.SeqIO.parse(input_fasta_handle, "fasta"):
        transcript_id = re.sub("\\.p[0-9]+\$", "", record.id)
        record.description = "transdecoder " + re.search("ORF type:[^,]+,score=[^,]+", record.description).group(0)
        if record.id in blastp_annotations:
          record.description += "; blastp: " + blastp_annotations.get(record.id)
        if record.id in pfam_annotations:
          record.description += "; pfam: " + pfam_annotations.get(record.id)
        if record.id in tmhmm_annotations:
          record.description += "; tmhmm: " + tmhmm_annotations.get(record.id)
        if record.id in deeploc_annotations:
          record.description += "; deeploc: " + deeploc_annotations.get(record.id)
        if record.id in signalp_annotations:
          record.description += "; signalp: " + signalp_annotations.get(record.id)
        Bio.SeqIO.write(record, output_fasta_handle, "fasta")

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

