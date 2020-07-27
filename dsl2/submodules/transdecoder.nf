nextflow.enable.dsl=2

process transdecoderLongOrfs {
  publishDir "transXpress_results", mode: "copy"
  queue params.queue_standard_nodes
  time params.queue_stdtime
  cpus 1
  memory "5 GB"
  input:
    val(assemblerName)
    path(transcriptomeTransdecoder) // transcriptomeToTransdecoder
  output:
    path "${transcriptomeTransdecoder}.transdecoder_dir/*.pep", emit:longOrfsProteomeSplit
    path "${transcriptomeTransdecoder}.transdecoder_dir/*", emit:transdecoderLongOrfsDirFiles
    path "${transcriptomeTransdecoder}.transdecoder_dir.__checkpoints_longorfs/*", emit:longOrfsCheckpointsFiles
    path "*.cmds", emit:longOrfsRootCmds
  tag { dateMetadataPrefix+"${assemblerName}" }
  script:
    """
    TransDecoder.LongOrfs -t ${transcriptomeTransdecoder} -m 30 ##Minimum protein length = -m amino acids
    #chmod -R a-w ${transcriptomeTransdecoder}.transdecoder_dir/ ##write protect the output to troubleshoot downstream accessing.
    """
}

process transdecoderPredict {
  publishDir "transXpress_results", mode: "copy" // , saveAs: { filename -> "transcriptome_after_predict.pep" }
  cpus 1
  queue params.queue_standard_nodes
  time { task.attempt > 1 ? params.queue_stdtime : params.queue_shorttime }
  input:
    val(assembler)
    path(transcriptomeFile) // transcriptomeTransdecoderPredict
    path "${transcriptomeFile}.transdecoder_dir/*" // transdecoderLongOrfsDirFiles
    path "${transcriptomeFile}.transdecoder_dir.__checkpoints_longorfs/*" // longOrfsCheckpointsFiles
    path longOrfsRootCmds
    //path blastpForTransdecoder
    //path pfamForTransdecoder
  output:
    tuple path("${transcriptomeFile}.transdecoder.pep"),path("${transcriptomeFile}.transdecoder.bed"),path("${transcriptomeFile}.transdecoder.cds"),path("${transcriptomeFile}.transdecoder.gff3")
  tag { dateMetadataPrefix+"${assembler}" }
  script:
    """
    ##TransDecoder.Predict -t ${transcriptomeFile} --retain_pfam_hits !{pfamForTransdecoder} --retain_blastp_hits !{blastpForTransdecoder}
    TransDecoder.Predict -t ${transcriptomeFile}
    """
}

workflow transdecoderFull {
take: transcriptomeTuple
main:
 transdecoderLongOrfs(transcriptomeTuple.map{it[0]},transcriptomeTuple.map{it[2]})

 transdecoderPredict(transcriptomeTuple.map{it[0]},transcriptomeTuple.map{it[2]},\
   transdecoderLongOrfs.out.transdecoderLongOrfsDirFiles,\
   transdecoderLongOrfs.out.longOrfsCheckpointsFiles,\
   transdecoderLongOrfs.out.longOrfsRootCmds)   
emit:
 transdecoderPredict.out 
}

