library(treeio)
library(ggtree)
library(readxl)
library(ggplot2)
library(ggnewscale)
library(castor)
library(phangorn)
library(ggprism)
library(tidytree)
library(tidyverse)
library(viridis)
library(RColorBrewer)

setwd("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Biomarkers paper")

tree = read.newick("core_SNPs_matrix.biomarker.contree")
tree$tip.label = gsub("_contigs.fasta", "", tree$tip.label)
tree$tip.label = gsub(".fasta", "", tree$tip.label)
tree_mid = midpoint(tree)

df = read_csv("genes_030326\\df_all_genes_030526.csv", col_types = c("ccnncfffffffff")) %>%
  select(-...1)
  
genes = df %>%
  select(acc, nheA:cytK2) %>%
  column_to_rownames(var = "acc") 

tox = df %>%
  select(acc, cytotoxicity) %>%
  column_to_rownames(var = "acc") 

p = ggtree(tree_mid) + 
  geom_tiplab(size=1) + 
  geom_treescale(linesize=0.75, offset = 1.2, x=0, y=200, width=0.05, fontsize=3.75) +
  geom_nodepoint(aes(fill = as.numeric(label)), size = 1, shape = 21) + 
  scale_fill_gradient(low = "white", high = "black", name = "Bootstrap\npercentage") +
  new_scale_fill() 

gheatmap(p, genes, width=0.04, offset=0.07, font.size=1.5, color=NA, colnames_angle = 45, colnames_offset_y = -0.75)

+
  scale_fill_manual(values = c("0" = "gray50", "1" = "#b1c6f0" , "2" = "#5b7bba", "3" = "#1e3c75"), na.value = "white", name = "Gene copy number") +
  new_scale_fill()

gheatmap(p, tox, width=0.03, offset=0.05, font.size=1.5, color=NA, colnames = FALSE) +
  scale_fill_gradient(low = "white", high = "red", name = "Cytotoxicity", na.value = "white")

