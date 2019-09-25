process downloadUniprotViruses {
  executor 'local'
  storeDir "./db"
  errorStrategy 'finish'
  output:
    file "virusesUniprot.pep.fasta" into uniref_virusFasta
  script:
    """
    wget -t 3 -O virusesUniprot.pep.fasta.gz "https://www.uniprot.org/uniprot/?query=taxonomy%3A%22Viruses+%5B10239%5D%22&format=fasta&compress=yes"
    zcat virusesUniprot.pep.fasta.gz | seqkit sort -n > virusesUniprot.pep.fasta
    """
}

process downloadNCBIRefSeqViruses {
  executor 'local'
  storeDir "./db"
  errorStrategy 'finish'
  output:
    file "ncbi_refseq_virus.all.protein.faa" into ncbi_virusFasta
  script:
    """
    wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/*protein.faa.gz
    zcat *.protein.faa.gz | seqkit sort -n > ncbi_refseq_virus.all.protein.faa
    """
}

ncbi_virusFasta.mix(uniref_virusFasta).into{virusFastas_ch1;virusFastas_ch2}

process filterOutPhages {
  executor 'local'
  storeDir "./db"
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
  executor 'local'
  storeDir "./db"
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
  storeDir "./db"
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
  storeDir "./db"
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
 storeDir "./db"
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

