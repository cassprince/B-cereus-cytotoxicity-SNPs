library(Biostrings)
library(msa)
library(bios2mds)
library(tidyverse)
library(seqinr)

setwd("C:/Users/cassp/OneDrive/Documents/Kovac Lab/Biomarkers paper/genes_030326/fastas")

files = list.files(path = "C:/Users/cassp/OneDrive/Documents/Kovac Lab/Biomarkers paper/genes_030326/fastas", pattern = "*.fasta", full.names = FALSE)

alignment2Fasta <- function(alignment, filename) {
  sink(filename)
  
  n <- length(rownames(alignment))
  for(i in seq(1, n)) {
    cat(paste0('>', rownames(alignment)[i]))
    cat('\n')
    the.sequence <- toString(unmasked(alignment)[[i]])
    cat(the.sequence)
    cat('\n')  
  }
  
  sink(NULL)
}

lapply(files, function(x) {
  name = strsplit(x, "_all")[[1]]
  print(name[1])
  
  gene_seqs = readDNAStringSet(x)
  gene_msa_o = msa(gene_seqs, method = c("ClustalOmega"), type = "dna")
  alignment2Fasta(gene_msa_o, paste0("alignments/", name[1], "_alignment.fasta"))
  
  #gene_fasta = msaConvert(gene_msa_o, "bios2mds::align")
  #export.fasta(gene_fasta, outfile = paste0("alignments/", name[1], "_alignment.fasta"), ncol = 80, open = "w")
})

