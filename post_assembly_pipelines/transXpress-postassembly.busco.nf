

Channel.from("bacteria_odb9",\
"fungi_odb9",\
"metazoa_odb9",\
"eukaryota_odb9").set{ BUSCO_lineages }

Channel.fromPath(params.peptides).set{peptides_ch}
Channel.fromPath(params.transcripts).set{transcripts_ch)

process downloadBUSCOProfile {
storeDir params.storeDB
input:
 val lineage from BUSCO_lineages
output:
 file "${lineage}" into BUSCO_profiles

script:
"""
wget https://busco.ezlab.org/datasets/${lineage}.tar.gz
tar -xvf ${lineage}.tar.gz
"""
}

BUSCO_profiles.into { BUSCO_lineages_ch1 ; BUSCO_lineages_ch2 }

BUSCO_cmds_pep = BUSCO_lineages_ch1.combine(peptides_ch)
BUSCO_cmds_trans = BUSCO_lineages_ch2.combine(transcripts_ch)

BUSCO_cmds_mixed = BUSCO_cmds_pep.mix(BUSCO_cmds_trans)

process do_BUSCO {
 conda "busco"
 publishDir "transXpress_results", mode: "copy"
 cpus 6
 input:
     set file(BUSCO_lineage), file(inputFasta) from BUSCO_cmds_mixed
 output:
     file 'run_'+assemblyPrefix+'*/*'
 tag { inputFasta+"_"+BUSCO_lineage }
 """
 NAME=\$(basename ${inputFasta})
 LINEAGE_NAME=\$(basename ${BUSCO_lineage})
 
 if [[ \$NAME == *'.pep'* ]]; then
  TYPE=prot
 else
  TYPE=transcriptome
 fi

 run_BUSCO.py -z -t /mnt/ramdisk -f -i ${inputFasta} -l ${BUSCO_lineage} -o \${NAME}_\${LINEAGE_NAME} -m \$TYPE --cpu ${task.cpus}
 """
 
}
