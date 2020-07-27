nextflow.enable.dsl=2
/*
 * Step 0 Download reference databases
 */

process downloadEggNOG {
  label 'stored'
  output:
    path "NOG.annotations.tsv" into eggNOGDb
  script:
    """
    wget "http://eggnogdb.embl.de/download/latest/data/NOG/NOG.annotations.tsv.gz"
    gunzip NOG.annotations.tsv.gz
    """
}

process downloadRfam {
  label 'stored'
  output:
    set path("Rfam.cm"), path("Rfam.cm.???")
  script:
    """
    wget "ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz"
    gunzip Rfam.cm.gz
    cmpress Rfam.cm
    """
}

process downloadSprot {
  label 'stored'
  output:
    set path("uniprot_sprot.fasta"), path("uniprot_sprot.fasta.p??")
  script:
    """
    wget "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
    gunzip uniprot_sprot.fasta.gz
    makeblastdb -in uniprot_sprot.fasta -dbtype prot
    """
}

process downloadPfam {
  label 'stored'
  output:
    set path("Pfam-A.hmm"), path("Pfam-A.hmm.h??")
  script:
    """
    wget "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
    gunzip Pfam-A.hmm.gz
    hmmpress Pfam-A.hmm
    """
}

workflow downloadTransxpressDbs {
main:
 downloadRfam()
 downloadSprot()
 downloadPfam()
 //downloadEggNOG()
emit:
 rfam=downloadRfam.out;sprot=downloadSprot.out;pfam=downloadPfam.out
}

workflow {
 downloadTransxpressDbs() 
}
