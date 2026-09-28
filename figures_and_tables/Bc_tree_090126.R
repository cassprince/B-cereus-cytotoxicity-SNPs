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


setwd("C:/Users/cassp/OneDrive/Documents/GitHub/B-cereus-cytotoxicity-SNPs")

tree = read.newick("C:/Users/cassp/OneDrive - Cornell University/Biomarkers paper/Biomarkers_RaxML_Tree.txt")
tree_mid = midpoint(tree)

df = read_csv("df_all_genes_090126.csv", col_types = c("cncfffffffffnnn"))
  
genes = df %>%
  select(acc, nheA:cytK2) %>%
  column_to_rownames(var = "acc") 

quast = df %>%
  select(acc, num_contigs) %>%
  column_to_rownames(var = "acc") 

tox = df %>%
  select(acc, cytotoxicity) %>%
  column_to_rownames(var = "acc") 

phylo = df %>%
  select(acc, panC_group) %>%
  column_to_rownames(var = "acc") 

p = ggtree(tree_mid) + 
  geom_tiplab(size=1) + 
  geom_treescale(linesize=0.75, offset = 1.2, x=0, y=200, width=0.05, fontsize=3.75) +
  geom_nodepoint(aes(fill = as.numeric(label)), size = 1, shape = 21) + 
  scale_fill_gradient(low = "white", high = "black", name = "Bootstrap\npercentage") +
  new_scale_fill() 

p1 = gheatmap(p, genes, width=0.04, offset=0.05, font.size=1.5, color=NA, colnames_angle = 45, colnames_offset_y = -0.75)+
  scale_fill_manual(values = c("0" = "gray50", "1" = "#b1c6f0" , "2" = "#5b7bba", "3" = "#1e3c75"), na.value = "white", name = "Gene copy number") +
  new_scale_fill()

p2 = gheatmap(p1, quast, width=0.03, offset=0.03, font.size=1.5, color=NA, colnames = FALSE) +
  scale_fill_gradient(low = "green", high = "red", name = "number of contigs", na.value = "white") + 
  new_scale_fill()


p2 = gheatmap(p1, tox, width=0.03, offset=0.03, font.size=1.5, color=NA, colnames = FALSE) +
  scale_fill_gradient(low = "white", high = "red", name = "Cytotoxicity", na.value = "white") + 
  new_scale_fill()

gheatmap(p2, phylo, width=0.009, offset=0.08, font.size=1.5, color=NA, colnames_angle = 45, colnames_offset_y = -0.75)

ggsave("quast_tree_090126.png", units = "in", width = 7, height = 8, dpi = 600)
