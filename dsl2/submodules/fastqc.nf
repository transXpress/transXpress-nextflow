nextflow.enable.dsl=2

process fastqc {
    tag "FASTQC on $sample_id"
    input:
     tuple sampleid,reads

    output:
     file("fastqc_${sample_id}_logs")


    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}
