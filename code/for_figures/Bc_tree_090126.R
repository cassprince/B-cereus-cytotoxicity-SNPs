#library(treeio)
library(ggtree)
library(ggnewscale)
library(castor)
library(phangorn)
library(tidytree)
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(phytools)


setwd("C:/Users/princecr/Documents/GitHub/B-cereus-cytotoxicity-SNPs/data")

tree = read.newick("Biomarkers_RaxML_tree_081526.txt")
tree_mid = midpoint(tree)

df = read_csv("df_all_genes_090126.csv", col_types = c("cncfffffffffnnn"))
df_snps = read_csv("snps_all_signif_pres_abs_102126.csv") %>%
  select(Isolate, nheA_up_SNP_t_312, nheA_up_SNP_c_418, hblC_up_SNP_g_561,  hblC_up_SNP_t_355, hblA_SNP_E_210, hblB_SNP_T_270) %>%
  mutate(across(where(is.numeric), as.character))

  
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

snps = df_snps %>%
  column_to_rownames(var = "Isolate") 

p = ggtree(tree_mid) + 
  geom_tiplab(size=1) + 
  geom_treescale(linesize=0.75, offset = 1.2, x=0, y=200, width=0.05, fontsize=3.75) +
  geom_nodepoint(aes(fill = as.numeric(label)), size = 1, shape = 21) + 
  scale_fill_gradient(low = "white", high = "black", name = "Bootstrap\npercentage") +
  new_scale_fill() 

p1 = gheatmap(p, genes, width=0.09, offset=0.05, font.size=1.5, color=NA, colnames = FALSE)+
  scale_fill_manual(values = c("0" = "gray50", "1" = "#b1c6f0" , "2" = "#5b7bba", "3" = "#1e3c75"), na.value = "white", name = "Gene copy number") +
  new_scale_fill()


p2 = gheatmap(p1, tox, width=0.03, offset=0.03, font.size=1.5, color=NA, colnames = FALSE) +
  scale_fill_gradient(low = "white", high = "red4", name = "Cytotoxicity", na.value = "white") + 
  new_scale_fill()

p3 = gheatmap(p2, snps, width=0.06, offset=0.09, font.size=1.5, color=NA, colnames = FALSE) +
  scale_fill_manual(values = c("0" = "#d5efd5", "1" = "#5a935a"), na.value = "gray50", name = "SNP presence") +
  new_scale_fill()

p3

ggsave("quast_tree_110426.png", units = "in", width = 7, height = 8, dpi = 600)
