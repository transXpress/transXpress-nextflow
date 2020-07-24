nextflow.preview.dsl=2

process rnaSPAdes {

    //label = "conda"

    memory = 20.GB
    cpus = 10

    input:
      file reads
      file datasets

    output:
      tuple val("rnaSPAdes"), \
          file("rnaSPAdes.gene_trans_map"), \
          file("rnaSPAdes/transcripts.fasta") \

    script:
      """
      rnaspades.py --dataset ${datasets} \
                   -t ${task.cpus} \
                   -m ${task.memory.toGiga()} \
                   -o rnaSPAdes \
                   --only-assembler -k 47 \
                   ${params.rnaSPAdes_params}
      # Make a fake gene to transcript file:
      cat "rnaSPAdes/transcripts.fasta" | grep ">" | tr -d ">" | cut -f 1 -d " " > tmp.txt
      paste tmp.txt tmp.txt > rnaSPAdes.gene_trans_map
      """
}
