Channel.fromPath(params.query).set{queryFasta}

params.store_dir = "./db"

process downloadUniprotViruses {
  executor 'local'
  storeDir params.store_dir
  errorStrategy 'finish'
  output:
    file "uniprot_total_viruses.dat.gz" into uniprot_total_viruses_swissFmt_ch1, uniprot_total_viruses_swissFmt_ch2
  script:
    """
    ##See here for parsing awk oneliner: https://www.biostars.org/p/153531/#154847
    ##wget -t 3 -O virusesUniprot.pep.fasta.gz "https://www.uniprot.org/uniprot/?query=taxonomy%3A%22Viruses+%5B10239%5D%22&format=fasta&compress=yes"
    wget -t 3 -O uniprot_sprot_viruses.dat.gz "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_viruses.dat.gz"
    wget -t 3 -O uniprot_trembl_viruses.dat.gz "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_viruses.dat.gz"
    cat uniprot_sprot_viruses.dat.gz >> uniprot_trembl_viruses.dat.gz
    rm -f uniprot_sprot_viruses.dat.gz
    mv uniprot_trembl_viruses.dat.gz uniprot_total_viruses.dat.gz
    """
}

process dumpUniprotFasta {
  executor 'local'
  storeDir params.store_dir
input: 
 file swissFmtDb from uniprot_total_viruses_swissFmt_ch1
output:
    file "virusesUniprot.pep.fasta" into uniref_virusFasta
script:
"""
#!/usr/bin/env python
##Decided to do it with Python, as the Awk script
##doesn't deal with multi-accessions well.
import gzip
import Bio
import Bio.SeqIO

write_handle = open("virusesUniprot.pep.fasta","w")
read_handle = gzip.open("${swissFmtDb}")
for r in Bio.SeqIO.parse(read_handle,"swiss"):
    Bio.SeqIO.write(r,write_handle,"fasta")
write_handle.close()
read_handle.close()
##zcat ${swissFmtDb} | awk '{if (/^ /) {gsub(/ /, ""); print} else if (/^AC/) print ">" \$2}' | sed 's/;\$//g' > virusesUniprot.pep.fasta
"""
}
uniref_virusFasta.into{virusFastas_ch1;virusFastas_ch2}


process downloadNCBItaxdump {
executor 'local'
storeDir params.store_dir
output:
 file "data/taxdump.tar.gz" into ncbi_taxonomy_dump_ch1,ncbi_taxonomy_dump_ch2
 file "data/nodes.dmp" into ncbi_nodes_ch1, ncbi_nodes_ch2
 file "data/names.dmp" into ncbi_names_ch1, ncbi_names_ch2
 file "data/merged.dmp" into ncbi_tax_merged_ch1
 file "data/delnodes.dmp" into ncbi_delnodes_ch1
script:
"""
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P data/
tar zxf data/taxdump.tar.gz -C data/ nodes.dmp names.dmp merged.dmp delnodes.dmp
"""
}

process setupETE3TaxonomyDatabase {
conda "ete3 biopython"
storeDir params.store_dir
input:
 file ncbi_taxonomy from ncbi_taxonomy_dump_ch1
output:
 file "ete3_taxa.sqlite" into ete3_taxonomy_database
script:
"""
#!/usr/bin/env python
from ete3 import Tree
from ete3 import NCBITaxa
from pathlib import Path
Path('ete3_taxa.sqlite').touch()
ncbi = NCBITaxa("ete3_taxa.sqlite")
ncbi.update_taxonomy_database("${ncbi_taxonomy}")
"""
}

//This takes about 35 GB of ram, for the full Uniprot viral swiss '.dat' database
process splitSwissFmtByTen {
scratch 'ram-disk'
stageOutMode 'move'
cpus 10
memory "2.5 MB" //I think LSF might be interpreting MB as GB currently?
input:
 file swissDb from uniprot_total_viruses_swissFmt_ch2
output:
 file "*.swiss.gz" into splitUniprotVirusSwissFmtDb
script:
"""
csplit --digits=16  --quiet --prefix=outfile --elide-empty-files <(zcat ${swissDb}) "/^ID/" "{*}"
ls -1 | grep "outfile" > files.txt 
grep -P "outfile...............0" files.txt | xargs cat | gzip > 0.${swissDb}.swiss.gz &
grep -P "outfile...............1" files.txt | xargs cat | gzip > 1.${swissDb}.swiss.gz &
grep -P "outfile...............2" files.txt | xargs cat | gzip > 2.${swissDb}.swiss.gz &
grep -P "outfile...............3" files.txt | xargs cat | gzip > 3.${swissDb}.swiss.gz &
grep -P "outfile...............4" files.txt | xargs cat | gzip > 4.${swissDb}.swiss.gz &
grep -P "outfile...............5" files.txt | xargs cat | gzip > 5.${swissDb}.swiss.gz &
grep -P "outfile...............6" files.txt | xargs cat | gzip > 6.${swissDb}.swiss.gz &
grep -P "outfile...............7" files.txt | xargs cat | gzip > 7.${swissDb}.swiss.gz &
grep -P "outfile...............8" files.txt | xargs cat | gzip > 8.${swissDb}.swiss.gz &
grep -P "outfile...............9" files.txt | xargs cat | gzip > 9.${swissDb}.swiss.gz &
wait
cat files.txt | xargs rm -f
"""
}

//This was too slow in my experience
//The new parallelized parseUniprotTaxonMapping is respectably fast
//But does have some high memory requirements
//process downloadUniprotTaxonMapping {
//storeDir params.store_dir
//output:
// file "uniprot2taxonid.tsv" into taxon_lookup_file
//script:
//"""
//wget -O - "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz" | zless | grep -P "\tNCBI_TaxID\t" > uniprot2taxonid.tsv
//"""
//}

splitUniprotVirusSwissFmtDb.flatten().combine(ete3_taxonomy_database).set{taxonParseTuples}

process parseUniprotTaxonMapping {
conda "ete3 biopython"
cache "lenient"
scratch 'ram-disk'
stageInMode 'copy'
stageOutMode 'move'
memory "2.5 MB" //I think LSF might be interpreting MB as GB  currently?
input:
 set file(virusSwissFile),file(ete3_taxonomy) from taxonParseTuples
output:
 file "*taxon_lookup.tsv" into taxon_lookup_files
script:
"""
#!/usr/bin/env python
import gzip
import Bio
from Bio import SeqIO
from ete3 import Tree
from ete3 import NCBITaxa
import os

ncbi = NCBITaxa("${ete3_taxonomy}")
write_handle = open("${virusSwissFile}_taxon_lookup.tsv","w",1000000) ##1 MB write buffer
read_handle = gzip.open("${virusSwissFile}")
for record in SeqIO.parse(read_handle, "swiss"):
    assert len(record.annotations['ncbi_taxid']) == 1
    
    taxid = record.annotations['ncbi_taxid']
    virus_obj = ncbi.get_taxid_translator(taxid)
    keys = list(virus_obj)
    if len(keys) > 0:
        theKey = keys[0]
        virus_name = virus_obj[theKey]
    else:
        virus_name = "virus name parse error"
        
    if 'host_ncbi_taxid' in record.annotations.keys():
        host_taxid = record.annotations['host_ncbi_taxid']
        ##print(record.id,taxid,host_taxid,record.name)
    else:
        ##Derive our own host species identifier, from the name of the virus
        #putative_genus = virus_name.split(" ")[0]
        #putative_species = virus_name.split(" ")[1]
        #genus_species =[putative_genus+" "+putative_species]
        #looked_up_taxon = ncbi.get_name_translator(genus_species)
        looked_up_taxon = False ######<<<< This disables the lookup, basically
        ##It is disabled, as it works for simple cases "Drosophila melanogaster", but not others
        if bool(looked_up_taxon) == False:
            ##No hit
            host_taxid = ['N/A']
        else:
            ##Found a hit
            theName = list(looked_up_taxon)[0]
            host_taxid = looked_up_taxon[theName]
    theString = "\t".join([record.id,taxid[0],str(host_taxid),record.name,virus_name])
    write_handle.write(theString+os.linesep)
write_handle.close()
read_handle.close()
"""
}

process mergeTaxonFiles {
storeDir params.store_dir
input:
 file taxonFiles from taxon_lookup_files.collect()
output:
 file "uniprot-virus2taxid.tsv" into merged_taxon_lookup_file_ch1,merged_taxon_lookup_file_ch2
script:
"""
ls -1L | grep "taxon_lookup.tsv" | xargs cat | sort > uniprot-virus2taxid.tsv
"""
}


process clone_BlobTools {
conda "anaconda matplotlib docopt tqdm wget pyyaml git pysam"
storeDir params.store_dir
input:
 file nodes from ncbi_nodes_ch1
 file names from ncbi_names_ch1
output:
 file "blobtools" into blobtools_dir
script:
"""
git clone https://github.com/DRL/blobtools.git
./blobtools/blobtools nodesdb --nodes ${nodes} --names ${names}
"""
}

process downloadNCBIRefSeqViruses {
  executor 'local'
  storeDir params.store_dir
  errorStrategy 'finish'
  output:
    file "ncbi_refseq_virus.all.protein.faa" into ncbi_virusFasta
  script:
    """
    wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/*protein.faa.gz
    zcat *.protein.faa.gz | seqkit sort -n > ncbi_refseq_virus.all.protein.faa
    rm -f *.protein.faa.gz
    """
}

//ncbi_virusFasta.mix(uniref_virusFasta).into{virusFastas_ch1;virusFastas_ch2}

process filterOutPhages {
  conda "seqkit"
  executor 'local'
  storeDir params.store_dir
  errorStrategy 'finish'
  input:
    file inputFasta from virusFastas_ch1
  output:
    file("no-phages_${inputFasta}") into notPhagesFastas
  script:
    """
    seqkit grep -vnrp "[Pp]hage" ${inputFasta} > no-phages_${inputFasta}
    """
}

process filterToPhages {
  conda "seqkit"
  executor 'local'
  storeDir params.store_dir
  errorStrategy 'finish'
  input:
   file inputFasta from virusFastas_ch2
  output:
    file "only-phages_${inputFasta}"
  script:
    """
    seqkit grep -nrp "[Pp]hage" ${inputFasta} > only-phages_${inputFasta}
    """
}

notPhagesFastas.into{ notPhagesFastas_ch1 ; notPhagesFastas_ch2 ; notPhagesFastas_ch3 ;notPhagesFastas_ch4 }

process makeDiamondDatabases {
  conda "diamond"
  storeDir params.store_dir
input:
 file(inputFasta) from notPhagesFastas_ch1
output:
 file "${inputFasta}.dmnd" into dmnd_Dbs
script:
"""
diamond makedb --in ${inputFasta} -d ${inputFasta}
"""
}

process makeBlastpDatabases {
  conda "blast"
  storeDir params.store_dir
input:
 file(inputFasta) from notPhagesFastas_ch2
output:
  file "${inputFasta}.???" into blast_Dbs
script:
"""
makeblastdb -in ${inputFasta} -dbtype prot
"""
}

process makeMMSeqsDatabases {
 conda "mmseqs2"
 storeDir params.store_dir
input:
 file(inputFasta) from notPhagesFastas_ch3
 file ncbiNodes from ncbi_nodes_ch2
 file ncbiNodes from ncbi_names_ch2
 file ncbiMerged from ncbi_tax_merged_ch1
 file ncbiDelnodes from ncbi_delnodes_ch1
 file customTaxIDFile from merged_taxon_lookup_file_ch2
output:
 file("mmseqs_db") into mmseqs_Dbs
script:
"""
mkdir mmseqs_db tmp ncbi_taxonomy_files
ln -s ../nodes.dmp ncbi_taxonomy_files/nodes.dmp
ln -s ../names.dmp ncbi_taxonomy_files/names.dmp
ln -s ../merged.dmp ncbi_taxonomy_files/merged.dmp
ln -s ../delnodes.dmp ncbi_taxonomy_files/delnodes.dmp

mmseqs createdb ${inputFasta} mmseqs_db/${inputFasta}.mmseqsdb

mmseqs createtaxdb mmseqs_db/${inputFasta}.mmseqsdb tmp \
   --ncbi-tax-dump ncbi_taxonomy_files \
   --tax-mapping-file <(cat ${customTaxIDFile} | cut -f1,2 | sed 's/\t/\tNCBI_TaxID\t/g')

mmseqs prefilter mmseqs_db/${inputFasta}.mmseqsdb  mmseqs_db/${inputFasta}.mmseqsdb mmseqs_db/${inputFasta}.mmseqsdb_pref -s 1.0

"""
}

process clusterMMSeqsDatabase {
 conda "mmseqs2"
 storeDir params.store_dir
input:
 file mmseqs_Db from mmseqs_Dbs
output:
 file "clustered-75id_${mmseqs_Db}" into clustered_mmseqs_Dbs
script:
"""
mmseqs align ${mmseqs_Db}/*.mmseqsdb ${mmseqs_Db}/*.mmseqsdb ${mmseqs_Db}/*.mmseqsdb_pref resultDB_aln
mmseqs clust ${mmseqs_Db}/*.mmseqsdb clustered-75id_${mmseqs_Db} resultDB_aln tmp --min-seq-id 0.75 -c 0.8 --cov-mode 1
"""
}

//Should do a map/join to linkup the Dbs and the fasta files.

process splitQueryFasta {
conda "ucsc-fasplit"
input:
 file(inputFasta) from queryFasta
output:
 file("split/*.fa") into splitQueryFastas
script:
"""
mkdir split
faSplit about ${inputFasta} 10000000 split/
"""
}

splitQueryFastas.flatten().combine(dmnd_Dbs).set{dmnd_files}

process doDiamondSearch {
conda "diamond"
cpus 4
input:
 set file(query_file),file(dmnd_db_file) from dmnd_files
output:
 set val("${dmnd_db_file}"), file("*.tsv") into dmnd_results
script:
"""
DMND_DB_NAME=${dmnd_db_file}
DMND_DB_NAME=\${DMND_DB_NAME%.dmnd}
diamond blastx --more-sensitive -d \${DMND_DB_NAME} -q ${query_file} --outfmt 6 --threads ${task.cpus} -o ${query_file}-\${DMND_DB_NAME}.dmnd-matches.tsv
sleep 15
"""
}

dmnd_results.groupTuple().into{mergeableResults_ch;printme_ch}
printme_ch.subscribe{ println it }

process mergeSimilaritySearchResults {
publishDir "results"
input:
 set val(key), file(theFiles) from mergeableResults_ch
output:
 file "${key}_merged.tsv" into mergedResults_ch
script:
"""
ls -1L | grep ".tsv" > result_files.txt 
##cat `head -n 1 result_files.txt` | grep -P "^@" > ${key}_merged.header
cat result_files.txt | xargs cat | grep -vP "^@" | sort -g -k 11,11 >> ${key}_merged.tsv
"""
}

mergedResults_ch.combine(blobtools_dir).combine(merged_taxon_lookup_file_ch1).set{taxifyTuples}

process taxifyResults {
conda "anaconda matplotlib docopt tqdm wget pyyaml git pysam"
publishDir "results"
input:
 set file(results),file(blobtools),file(taxid_mapping) from taxifyTuples
script:
"""
${blobtools}/blobtools taxify --hit_file ${results} \
  --hit_column_qseqid 0 \
  --hit_column_sseqid 1 \
  --hit_column_score 11 \
  --taxid_mapping_file ${taxid_mapping} \
  --map_col_sseqid 0 \
  --map_col_taxid 1
"""
}

process bootstrapPavianContainer {
 storeDir params.store_dir
output:
 file "pavian-singularity.img" into pavian_singularity_img
script:
"""
export SINGULARITY_CACHEDIR=${params.store_dir}
singularity build --force --sandbox pavian-singularity.img docker://rocker/shiny-verse
singularity exec --cleanenv --writable pavian-singularity.img Rscript <(echo 'if (!require(remotes)) { install.packages("remotes") }')
singularity exec --cleanenv --writable pavian-singularity.img Rscript <(echo 'remotes::install_github("fbreitwieser/pavian")')
"""
}

process runPavian {
 input:
  file psi from pavian_singularity_img
script:
"""
export SINGULARITY_CACHEDIR=${params.store_dir}
##Note the extra backslashes for nextflow escaping
HOSTIP=`ip -4 addr show eno1 | grep -oP '(?<=inet\\s)\\d+(\\.\\d+){3}'`
#singularity exec --cleanenv --writable ${psi} Rscript <(echo 'pavian::runApp(port=13371,host="\${HOSTIP}")')
"""
}
