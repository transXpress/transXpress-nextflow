Channel.fromPath(["/lab/solexa_weng/Seq_data/Projects/Tim_Fallon/BUSCO_profiles/BUSCO_v2_profiles/bacteria_odb9/",\
"/lab/solexa_weng/Seq_data/Projects/Tim_Fallon/BUSCO_profiles/BUSCO_v2_profiles/fungi_odb9/",\
"/lab/solexa_weng/Seq_data/Projects/Tim_Fallon/BUSCO_profiles/BUSCO_v2_profiles/metazoa_odb9/",\
"/lab/solexa_weng/Seq_data/Projects/Tim_Fallon/BUSCO_profiles/BUSCO_v2_profiles/eukaryota_odb9/"]).set{ BUSCO_lineages }

BUSCO_lineages.into { BUSCO_lineages_ch1 ; BUSCO_lineages_ch2 }

BUSCO_cmds_pep = BUSCO_lineages_ch1.combine(transcriptome_annotated_pep_ch2)
BUSCO_cmds_trans = BUSCO_lineages_ch2.combine(transcriptome_annotated_fasta_ch2)

BUSCO_cmds_mixed = BUSCO_cmds_pep.mix(BUSCO_cmds_trans)

process do_BUSCO {
 publishDir "transXpress_results", mode: "copy"
 cpus 21

 input:
     set file(BUSCO_lineage), file(inputFasta) from BUSCO_cmds_mixed
 output:
     file 'run_'+assemblyPrefix+'*/*'
 tag { inputFasta+"_"+BUSCO_lineage }
 """
 #! /bin/bash
 ##A new script that uses BUSCO_v3

 ##Note that BUSCO_v3 
 #In short, I’ve cloned from the BUSCO_v3 github (https://gitlab.com/ezlab/busco) to -> /lab/solexa_weng/testtube/busco
 #I’ve installed it to my user directory (python setup.py install --user)
 #I’ve updated the paths & other config things in (/lab/solexa_weng/testtube/busco/config/config.ini), to the Tak paths for the appropriate software, & other settings

 export AUGUSTUS_CONFIG_PATH=/lab/solexa_weng/testtube/busco/augustus_config/ ##I copied the augustus config directory myself.

 BUSCO_SCRIPTS=/lab/solexa_weng/testtube/busco/scripts/

 NAME=\$(basename ${inputFasta})
 LINEAGE_NAME=\$(basename ${BUSCO_lineage})
 
 if [[ \$NAME == *'.pep'* ]]; then
  TYPE=prot
 else
  TYPE=transcriptome
 fi

 \${BUSCO_SCRIPTS}/run_BUSCO.py -z -t /mnt/ramdisk -f -i ${inputFasta} -l ${BUSCO_lineage} -o \${NAME}_\${LINEAGE_NAME} -m \$TYPE --cpu ${task.cpus} >./BUSCO.\${NAME}_\${LINEAGE_NAME}.stdout.log 2>./BUSCO.\${NAME}_\${LINEAGE_NAME}.stderr.log
 find /mnt/ramdisk -name \"*\${NAME}*\${LINEAGE_NAME}*\" | xargs rm -f
 """
 
}
