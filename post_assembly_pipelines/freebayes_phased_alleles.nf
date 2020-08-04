nextflow.enable.dsl=2

//See here for prototype implementation:
///lab/weng_scratch/Tim/PPYR_OGS/Supporting_non-OGS_data/Alleles

process bowtie2_index {
publishDir "results", mode:'link',overwrite:'true'
module 'bowtie2'
 input:
  path fasta
 output:
  tuple val("$fasta"), path("${fasta}.*.bt2")

script:
"""
bowtie2-build ${fasta} ${fasta}
"""
}

process align_PE {
publishDir "results", mode:'link',overwrite:'true'
module 'bowtie2:samtools'
//conda "bowtie2 samtools"
cpus 4
 input:
  tuple val(read_id), path(reads)
  tuple val(index_name), path(index_files)
 output:
  path "${read_id}-${index_name}.aligned_reads.wrg.bam"
tag "${read_id}-${index_name}"
script:
"""
##bowtie2 --threads ${task.cpus} -X 2000 --fr -x ${index_name} -1 ${reads[0]} -2 ${reads[1]} | samtools view -b -h -F 4 | samtools sort -@ 4 -T ${read_id}-${index_name} - | bamaddrg -b /dev/stdin -s ${read_id}-${index_name}.aligned_reads.wrg.bam -r ${read_id}-${index_name}.aligned_reads.wrg.bam > ${read_id}-${index_name}.aligned_reads.wrg.bam

bowtie2 --threads ${task.cpus} -X 2000 --fr -x ${index_name} -1 ${reads[0]} -2 ${reads[1]} --rg-id ${read_id} --rg "SM:${read_id}" | samtools view -b -h -F 4 | samtools sort -@ 4 -T ${read_id}-${index_name} -o ${read_id}-${index_name}.aligned_reads.wrg.bam
"""
}

process freebayes_call {
publishDir "results", mode:'link',overwrite:'true'
 conda 'freebayes'
 input:
  path reference_fasta
  path bam_files
 output:
  tuple path("all.vcf.gz"),path("all2.vcf.gz.tbi")
script:
"""
freebayes --max-complex-gap 75 --min-coverage 3 --fasta-reference ${reference_fasta} -p 6 ${bam_files} | bgzip > all.vcf.gz
tabix all.vcf.gz
"""
}

process bcf_conv {
publishDir "results", mode:'link',overwrite:'true'
 input:
  tuple path(vcf),path(vcf_index)
 output:
  path "all.bcf"
script:
"""
bcftools view ${vcf} -Ob -o all.bcf
"""
}

workflow {
 reads = Channel.fromFilePairs('../raw_data/SRR*_{1,2}.fastq.gz')
 transcripts = Channel.fromPath('../transXpress_results/2*.transcripts.fasta')
 
 bowtie2_index(transcripts)
 align_PE(reads,bowtie2_index.out.collect())
 freebayes_call(transcripts,align_PE.out.collect())
 bcf_conv(freebayes_call.out)
}
