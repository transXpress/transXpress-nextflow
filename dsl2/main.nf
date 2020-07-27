nextflow.enable.dsl=2

//params.reads = "$baseDir/raw_data/reads/*_{1,2}.fq*"

include downloadTransxpressDbs as downloadDbs from './submodules/databases.nf'
include parseSamplesTSV from './submodules/sample_parsing.nf'
include trim from './submodules/trimmomatic.nf'

include rnaSPAdes from './submodules/rnaspades.nf' 
include fullTrinity from './submodules/trinity.nf'
include renameAssembly from './submodules/rename.nf'

include transdecoderFull from './submodules/transdecoder.nf'
include transdecoderLongOrfs from './submodules/transdecoder.nf'


workflow {

//The below code ensures that if the workflow is rerun on a different day from the same input files, then 
//it doesn't recalculate the assembly_prefix
def outFile = new File(params.output_dir+'/assembly_prefix.txt')
if (!outFile.exists()) {

    def resultsDir = new File(params.output_dir)
    if (!resultsDir.exists()) {
        resultsDir.mkdirs() //Make the transXpress_results directory
    }


    theDate = ""
    if (params.prefix_add_date == true) {
        theDate = new java.util.Date().format(params.prefix_add_date_formatting) //yyMMdd by default
    }

    metadata = ""
    if (params.prefix_add_metadata_file != "") {
        theText = file(params.prefix_add_metadata_file).text
        metadata = theText.replace(" ", "_").trim()
    }

    dateMetadataToJoin = []
    if (theDate != "") {
        dateMetadataToJoin.add(theDate)
    }
    if (metadata != "") {
        dateMetadataToJoin.add(metadata)
    }

    dateMetadataPrefix = dateMetadataToJoin.join("_")
    if (dateMetadataPrefix != "") {
        dateMetadataPrefix += "_"
    }

    outFile.withWriter('UTF-8') {
        writer ->
            writer.write(dateMetadataPrefix)
    }
}

//Convert into a global string variable
dateMetadataPrefix = outFile.text

dateMetadataPrefix_ch = Channel.from(dateMetadataPrefix)

println(dateMetadataPrefix)

 downloadDbs()
 rfam_ch = downloadDbs.out.rfam
 pfam_ch = downloadDbs.out.pfam
 sprot_ch = downloadDbs.out.sprot

 parseSamplesTSV(params.samples)
 parseSamplesTSV.out.reads.view()
 parseSamplesTSV.out.yaml.view()

 trim(parseSamplesTSV.out.reads) 
 
 //Run one assembler.
 assemblerString = params.assembler.toLowerCase().trim()
 if (assemblerString == "trinity") {
   fullTrinity(trim.out.collect(),parseSamplesTSV.out.tsv)
   fullTrinity.out.view()
   finishedAssemblies = fullTrinity.out
 }
 else if (assemblerString == "rnaspades") {
   rnaSPAdes(trim.out.collect(),parseSamplesTSV.out.yaml)
   finishedAssemblies = rnaSPAdes.out
 }
 else if (assemblerString == "both") {
  fullTrinity(trim.out.collect(),parseSamplesTSV.out.tsv)
  rnaSPAdes(trim.out.collect(),parseSamplesTSV.out.yaml)
  finishedAssemblies = fullTrinity.out.mix(rnaSPAdes.out)
 }
 else {
   println "Error:"+params.assembler+" doesn't match an implemented assembler"
   return -1 //Error status
 }
 
 renameAssembly(finishedAssemblies,dateMetadataPrefix_ch)
 //renameAssembly.out.view()
 //transdecoderLongOrfs(renameAssembly.out[0],renameAssembly.out[2])
 transdecoderFull(renameAssembly.out) 
}

