nextflow.preview.dsl=2

process prepareTrimmomaticAdapters {
output:
 path "adapters.fa"
shell:
'''
seqkit seq -w 0 ${CONDA_PREFIX}/share/trimmomatic-*/adapters/*.fa > adapters.fa
'''
}

process fixReadNaming {
input:
 tuple path(R1_reads),path(R2_reads)
output:
 tuple path("fixed/${R1_reads}"),path("fixed/${R2_reads}")
script:
"""
mkdir fixed
seqkit replace -p "_forward" -r "" ${R1_reads} -o fixed/${R1_reads}
seqkit replace -p "_reverse" -r "" ${R2_reads} -o fixed/${R2_reads}
"""
}

process trimmomatic {
input:
 tuple file(R1_reads),file(R2_reads)
 file "adapters.fasta"
tag {"$R1_reads"+" and " +"$R2_reads"}
output:
  set file("${R1_reads}.R1-P.qtrim.fastq.gz"), file("${R2_reads}.R2-P.qtrim.fastq.gz") //filteredPairedReads
  //file "*U.qtrim.fastq.gz" // filteredSingleReads
script:
"""
trimmomatic PE -threads ${task.cpus} ${R1_reads} ${R2_reads} ${R1_reads}.R1-P.qtrim.fastq.gz ${R1_reads}.R1-U.qtrim.fastq.gz ${R2_reads}.R2-P.qtrim.fastq.gz ${R2_reads}.R2-U.qtrim.fastq.gz ${params.TRIMMOMATIC_PARAMS} 
"""
}

workflow trim {
take: reads
main:
prepareTrimmomaticAdapters()
fixReadNaming(reads)
trimmomatic(fixReadNaming.out,prepareTrimmomaticAdapters.out)
emit:
 trimmomatic.out
}

