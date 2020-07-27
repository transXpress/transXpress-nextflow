nextflow.enable.dsl=2

process trinityInchwormChrysalis {
  //See here for an alternative approach:
  //https://github.com/biocorecrg/transcriptome_assembly/blob/564f6af2e4db9625ae9de6884a6524b4ec57cece/denovo_assembly/denovo_assembly.nf#L157
  //label "conda"
  cache 'lenient'
  cpus params.assembly_CPUs
  memory params.assembly_MEM+" GB"
  queue params.queue_highmemory_nodes
  time params.queue_longtime
  stageInMode "copy"
  scratch "ram-disk"
  tag { dateMetadataPrefix+"Trinity" }

  afterScript 'echo \"(Above completion message is from Trinity. transXpress will continue the pipeline execution.)\"'
  input:
     path trinityInchwormPairedReads
     path samples_file
  output:
    path "trinity_out_dir/[!Tcr]*", emit: trinityWorkDirRootFiles //Not files starting with c or r, so not chrysalis, read_partitions, recursive trinity cmds, Trinity.timing
    path "trinity_out_dir/chrysalis/*", emit: trinityWorkDirChrysalisFiles
    path "trinity_out_dir/read_partitions/**.trinity.reads.fa", emit: trinityPhase1ReadPartitionsFiles
    //

  script:
    """
    Trinity --no_distributed_trinity_exec --max_memory ${task.memory.toGiga()}G --CPU ${task.cpus}\
         --samples_file ${samples_file}\
         --seqType ${params.TRINITY_SEQTYPE}\
         --min_glue ${params.TRINITY_MIN_GLUE}\
         --min_kmer_cov ${params.TRINITY_MIN_KMER_COV}\
         ${params.TRINITY_OTHER}
    #chmod -R a-w ./trinity_out_dir ##<- make the results read only, to troubleshoot Trinity processes writing into directories they shouldn't.
    """
}

process trinityButterflyParallel {
  //TODO This process has hardcoded parameters.  It should really be getting them from the TRINITY params...
  //See here for an alternative approach to this node:
  //https://github.com/biocorecrg/transcriptome_assembly/blob/564f6af2e4db9625ae9de6884a6524b4ec57cece/denovo_assembly/denovo_assembly.nf#L183
  //See also here, for more details on the parallelization of Trinity "Stage 2" https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity#optional-adapting-trinity-to-a-computing-grid-for-massively-parallel-processing-of-embarrassingly-parallel-steps
  cache 'lenient'
  cpus params.general_CPUs
  queue params.queue_standard_nodes
  time params.queue_stdtime
  stageInMode "copy"
  scratch "ram-disk"
  input:
    path("trinity_out_dir/*") // trinityWorkDirRootFiles_ch1
    path("trinity_out_dir/chrysalis/*") // trinityWorkDirChrysalisFiles_ch1
    tuple val(dir),path(fastaFiles)   // partitionedReadGroups_ch
  output:
   path "commands.completed", emit:trinityFinishedCmds
   path "commands.txt", emit:trinityCmds
   path "trinity_out_dir/read_partitions/${dir[0]}/${dir[1]}/*out.Trinity.fasta", emit:butterflyTrinityFiles
  tag { dateMetadataPrefix+"Trinity-"+dir[0]+"/"+dir[1]}
  script:
    """
    ##Have to recreate the directory structure for the read_partition files
    mkdir "trinity_out_dir/read_partitions"
    mkdir "trinity_out_dir/read_partitions/${dir[0]}"
    mkdir "trinity_out_dir/read_partitions/${dir[0]}/${dir[1]}"
    for f in ./*.fa
    do
     ##Once the directory structure is made, link the FASTA file back into it.
     ##Note the dollar sign escaping for nextflow
     ln -s "../../../../\$f" "trinity_out_dir/read_partitions/${dir[0]}/${dir[1]}/\$f"
    done
    for f in ./trinity_out_dir/read_partitions/${dir[0]}/${dir[1]}/*.fa
    do
     ##Make our own Trinity commands, using relative paths
     echo "Trinity --single '\$f' --output '\$f.out' --CPU 1 --max_memory 1G --run_as_paired --seqType fa --trinity_complete --full_cleanup --no_distributed_trinity_exec --min_glue 2 --min_kmer_cov 2" >> commands.txt
    done
    ##Execute in parallel
    parallel --jobs ${task.cpus} < commands.txt
    cp commands.txt commands.completed
    #chmod -R a-w ./trinity_out_dir
    sleep 15 ##Try to prevent filesystem latency errors
    """ 
}

process trinityFinish {
   cache "lenient"
   executor "local"
   stageInMode "copy"
   scratch "ram-disk"
   publishDir "transXpress_results", mode: "link", saveAs: { filename -> filename.replaceAll("trinity_out_dir/Trinity", "transcriptome") }
  input:
    path "trinity_out_dir/*" // trinityWorkDirRootFiles_ch2 //An attempt to relativize the butterfly processes
    path "trinity_out_dir/chrysalis/*" // trinityWorkDirChrysalisFiles_ch2 //An attempt to relativize the butterfly processes
    path butterflyTrinityFilesCollected // butterflyTrinityFiles.collect() 
    path relativeSamples // relativeSamples_toTrinityFinish
    path finishedCommands // trinityFinishedCmds.collectFile(name: "recursive_trinity.cmds.completed")
    path trinityCmdsCollected // trinityCmds.collectFile(name: "recursive_trinity.cmds")
    //path trinityFinishPairedReads
  output:
    tuple val("Trinity"), path("./trinity_out_dir/Trinity.fasta.gene_trans_map"),path("./trinity_out_dir/Trinity.fasta") // trinityFinalOutput
    //path "./trinity_out_dir/recursive_trinity.cmds.completed"
  memory "1 GB"
  cpus 1
  //tag { dateMetadataPrefix+"Trinity" }
  script:
    """
    ##Make some dummy read files, so we don't have to stageIn the actual (large) files
    cat ${relativeSamples} | cut -f 3,4 | tr "\t" "\n" | while read line; do echo ">dummy\nATCG" | gzip > \$line; done

    mkdir trinity_out_dir/read_partitions/
    mkdir trinity_out_dir/read_partitions/Fb_0/
    mkdir trinity_out_dir/read_partitions/Fb_0/CBin_0 ##This is just a dummy directory to fool Trinity
    sleep 15 ##Filesystem latency errors?
    ls -1L | grep ".fasta" > fasta_files.txt
    while read f; do
     ##Link the files into a directory structure that Trinity can deal with correctly, even if it isn't 100% right
     ln -s "../../../../\$f" "trinity_out_dir/read_partitions/Fb_0/CBin_0/."
    done < fasta_files.txt

    ##Have to produce these files to trick Trinity into thinking things are done
    cp ${trinityCmdsCollected} trinity_out_dir/recursive_trinity.cmds
    cp ${finishedCommands} trinity_out_dir/recursive_trinity.cmds.completed
    touch trinity_out_dir/recursive_trinity.cmds.ok
    Trinity --samples_file ${relativeSamples} --max_memory ${task.memory.toGiga()}G --seqType ${params.TRINITY_SEQTYPE}\
         --min_glue ${params.TRINITY_MIN_GLUE}\
         --min_kmer_cov ${params.TRINITY_MIN_KMER_COV}\
         ${params.TRINITY_OTHER}
    """ 
}

workflow fullTrinity {
take: reads; samples_file
main:
 trinityInchwormChrysalis(reads,samples_file)

//This is a workaround, of sorts, to parallelize the Trinity butterfly steps
//As Trinity usually expects to be writing into a single directory over time
//which doesn't play nice with Nextflow, have to do some tricks
//to recreate the structure of the Trinity work directory
//so each process ends up with a symbolic linked copy of the Trinity work directory
//and can write real files in the proper directory structure.
//But that directory structure can't be remade by Nextflow using Queue Channels, that I am aware of
//So the below mapping trick makes the directory structure into variables
//that the trinityButterflyParallel & trinityFinish processes can use to recreate the structure
trinityInchwormChrysalis.out.trinityPhase1ReadPartitionsFiles.flatten().map{ file ->
                                        def filePath = file.toAbsolutePath().toString()
                                        dir1Match = filePath =~ /Fb_[0-9]+/ ///Regex matching
                                        dir2Match = filePath =~ /CBin_[0-9]+/ ///Regex matching
					dir1String = dir1Match[0]
					dir2String = dir2Match[0]
 
                                        return tuple( [dir1String,dir2String] ,file)
                                        }
                                        .groupTuple()
                                        .set{ partitionedReadGroups_ch }
	//Add .groupTuple() to execute in groups by the directories

trinityButterflyParallel(trinityInchwormChrysalis.out.trinityWorkDirRootFiles.collect(),\
                         trinityInchwormChrysalis.out.trinityWorkDirChrysalisFiles.collect(),\
                         partitionedReadGroups_ch)

trinityFinish(trinityInchwormChrysalis.out.trinityWorkDirRootFiles.collect(),\
   trinityInchwormChrysalis.out.trinityWorkDirChrysalisFiles.collect(),\
   trinityButterflyParallel.out.butterflyTrinityFiles.collect(),\
   samples_file,\
   trinityButterflyParallel.out.trinityFinishedCmds.collectFile(name:"recursive_trinity.cmds.completed"),\
   trinityButterflyParallel.out.trinityCmds.collectFile(name: "recursive_trinity.cmds"))

emit:
 trinityFinish.out
}
