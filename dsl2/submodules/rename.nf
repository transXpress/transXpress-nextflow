process renameAssembly {
   executor 'local'
   storeDir "transXpress_results"
   input:
    tuple val(assembler), path(geneTransMap), path(transcriptome_fasta)
    val dateMetadataPrefix
   output:
    tuple assembler, path("${assembler}_renamed.fasta.gene_trans_map"), path("${dateMetadataPrefix}${assembler}.transcripts.fasta")
   tag { dateMetadataPrefix+"${assembler}" }
   script:
   """
   seqkit replace -p 'TRINITY_' -r '' ${transcriptome_fasta} | seqkit replace -p '^' -r '${dateMetadataPrefix}${assembler}_' > ${dateMetadataPrefix}${assembler}.transcripts.fasta
   cat ${geneTransMap} | sed 's/TRINITY_//g' | sed 's/^/${dateMetadataPrefix}${assembler}_/g' | sed 's/\t/\t${dateMetadataPrefix}${assembler}_/g' > ${assembler}_renamed.fasta.gene_trans_map
   """
}
