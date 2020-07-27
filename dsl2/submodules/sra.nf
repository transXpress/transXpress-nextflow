nextflow.enable.dsl=2

process fasterq_dump {
  //conda "sratools"
  scratch 'ram-disk'
  cpus 2
  input: 
   val SRA
  output:
    tuple path("./*_1.fastq.gz"),path("./*_2.fastq.gz")
  script:
    """
    fasterq-dump --print-read-nr --split-files --temp ./ ${SRA}
    ls -1 | grep ".fastq" | xargs -P 2 -n 1 gzip
    """
}

workflow {
 fasterq_dump("SRR931174")
}
