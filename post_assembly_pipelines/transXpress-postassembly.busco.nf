nextflow.enable.dsl=2
params.peptides = ""
params.transcripts = ""

process downloadBUSCOLineage {
input:
 val lineage
output:
 file "*"

script:
"""
wget https://busco-data.ezlab.org/v4/data/lineages/${lineage}.tar.gz
tar -xvf ${lineage}.tar.gz
rm -f ${lineage}.tar.gz
"""
}

process do_BUSCO {
 conda 'busco>=4.1'
 publishDir "results/busco", mode: "link"
 cpus 6
 input:
     tuple path(BUSCO_lineage), path(inputFasta)
 output:
     path 'run_*/*'
 tag { inputFasta+"_"+BUSCO_lineage }
 """
 NAME=\$(basename ${inputFasta})
 LINEAGE_NAME=\$(basename ${BUSCO_lineage})
 
 if [[ \$NAME == *'.pep'* ]]; then
  TYPE=prot
 else
  TYPE=transcriptome
 fi

 run_busco -z -f -i ${inputFasta} -l ${BUSCO_lineage} -o \${NAME}_\${LINEAGE_NAME} -m \$TYPE --cpu ${task.cpus}
 """ 
}

workflow {

Channel.from("eukaryota_odb10.2019-11-20").set{ BUSCO_lineages }

//peptides_ch = Channel.fromPath(params.peptides)
transcripts_ch = Channel.fromPath(params.transcripts)

downloadBUSCOLineage(BUSCO_lineages)

mixed = downloadBUSCOLineage.out.combine(transcripts_ch)

do_BUSCO(mixed)

}
