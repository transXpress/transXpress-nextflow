Channel.fromPath(params.query).set{queryFasta}

params.store_dir = "./db"

process downloadUniprotViruses {
  executor 'local'
  storeDir params.store_dir
  errorStrategy 'finish'
  output:
    file "virusesUniprot.pep.fasta" into uniref_virusFasta
    file "uniprot_total_viruses.dat.gz" into uniprot_total_viruses
  script:
    """
    ##See here for parsing awk oneliner: https://www.biostars.org/p/153531/#154847
    ##wget -t 3 -O virusesUniprot.pep.fasta.gz "https://www.uniprot.org/uniprot/?query=taxonomy%3A%22Viruses+%5B10239%5D%22&format=fasta&compress=yes"
    wget -t 3 -O uniprot_sprot_viruses.dat.gz "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_viruses.dat.gz"
    wget -t 3 -O uniprot_trembl_viruses.dat.gz "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_viruses.dat.gz"
    cat uniprot_sprot_viruses.dat.gz >> uniprot_trembl_viruses.dat.gz
    rm -f uniprot_sprot_viruses.dat.gz
    mv uniprot_trembl_viruses.dat.gz uniprot_total_viruses.dat.gz
    zcat uniprot_total_viruses.dat.gz | awk '{if (/^ /) {gsub(/ /, ""); print} else if (/^ID/) print ">" \$2}' > virusesUniprot.pep.fasta
    """
}

process downloadNCBItaxdump {
storeDir params.store_dir
output:
 file "data/taxdump.tar.gz" into ncbi_taxonomy_dump
 file "data/nodes.dmp" into ncbi_nodes
 file "data/names.dmp" into ncbi_names
script:
"""
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P data/
tar zxf data/taxdump.tar.gz -C data/ nodes.dmp names.dmp
"""
}

process setupETE3TaxonomyDatabase {
conda "ete3 biopython"
storeDir params.store_dir
input:
 file ncbi_taxonomy from ncbi_taxonomy_dump
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

process downloadUniprotTaxonMapping {
storeDir params.store_dir
output:
 file "uniprot2taxonid.tsv" into uniprot_taxon_mapping
script:
"""
wget -O - "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz" | zless | grep -P "\tNCBI_TaxID\t" > uniprot2taxonid.tsv
"""
}

process clone_BlobTools {
conda "anaconda matplotlib docopt tqdm wget pyyaml git pysam"
storeDir params.store_dir
input:
 file nodes from ncbi_nodes
 file names from ncbi_names
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

uniref_virusFasta.into{virusFastas_ch1;virusFastas_ch2}
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
output:
 file("${inputFasta}.mmseqsdb*") into mmseqs_Dbs
script:
"""
mmseqs createdb ${inputFasta} ${inputFasta}.mmseqsdb
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

process mergeResults {
publishDir "results"
input:
 set val(key), file(theFiles) from mergeableResults_ch
output:
 file "${key}_merged.tsv" into mergedResults_ch
script:
"""
ls -1 | grep ".tsv" > result_files.txt 
##cat `head -n 1 result_files.txt` | grep -P "^@" > ${key}_merged.header
cat result_files.txt | xargs cat | grep -vP "^@" | sort -g -k 11,11 >> ${key}_merged.tsv
"""
}

mergedResults_ch.combine(blobtools_dir).combine(taxon_lookup_file).set{taxifyTuples}

process taxifyResults {
conda "anaconda matplotlib docopt tqdm wget pyyaml git pysam"
publishDir "results"
input:
 set file(results),file(blobtools),file(taxid_mapping) from taxifyTuples
script:
"""
TAXID_NAME="${taxid_mapping}"
TAXID_NAME_NOGZ="\${TAXID_NAME%.gz}"
zcat ${taxid_mapping} | cut -f 2,13 > \${TAXID_NAME_NOGZ}
${blobtools}/blobtools taxify --hit_file ${results} \
  --hit_column_qseqid 0 \
  --hit_column_sseqid 1 \
  --hit_column_score 11 \
  --taxid_mapping_file \${TAXID_NAME_NOGZ} \
  --map_col_sseqid 0 \
  --map_col_taxid 1

rm -f \${TAXID_NAME_NOGZ}
"""
}
