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

snps = read_csv("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Biomarkers paper\\SNP hits\\4_21_25\\SNP_hits_4_22_25.csv")
snps_nonsyn = snps %>% 
  select(Isolate, HblA_D194E:nhe_promoter_544T) 
snps_nonsyn = snps_nonsyn %>%
  mutate_at(colnames(snps_nonsyn), factor) %>%
  column_to_rownames(var = "Isolate") %>% 
  select(starts_with(c("H", "N"), ignore.case = FALSE)) 

snps_syn = snps %>% 
  select(Isolate, HblA_D194E:nhe_promoter_544T) 
snps_syn = snps_syn %>%
  mutate_at(colnames(snps_syn), factor) %>% 
  column_to_rownames(var = "Isolate") %>% 
  select(-starts_with(c("H", "N"), ignore.case = FALSE)) %>%
  as.data.frame()

genes = snps %>% 
  select(Isolate, nhe:cytk2) 
genes = genes %>%
  mutate_at(colnames(genes), factor) %>% 
  column_to_rownames(var = "Isolate") 

cytotoxicity = snps %>%
  select(Isolate, Average.Cell.Viability) %>% 
  column_to_rownames(var = "Isolate")

tree = read.newick("core_SNPs_matrix.biomarker.contree")
tree$tip.label = gsub("_contigs.fasta", "", tree$tip.label)
tree$tip.label = gsub(".fasta", "", tree$tip.label)
tree_mid = midpoint(tree)



###---phylogram rectangluar---###
p = ggtree(tree_mid, size = 0.5) + geom_tiplab(size=1) + geom_treescale(linesize=0.75, offset = 1.2, x=0, y=200, width=0.05, fontsize=3.75) +
  geom_nodepoint(aes(fill = as.numeric(label)), size = 1, shape = 21) + 
  scale_fill_gradient(low = "white", high = "black", name = "Bootstrap\npercentage") +
  new_scale_fill() 

p_genes = gheatmap(p, genes, width=0.04, offset=0.07, font.size=1.5, color=NA, colnames_angle = 45, colnames_offset_y = -0.75) +
  scale_fill_manual(values = c("1" = "#125c8a", "0" = "#c5d9ed"), na.value = "white", labels=c('Absent','Present'), name = "Gene Presence") +
  new_scale_fill() 

p_snps_syn = gheatmap(p_genes, snps_syn, width=0.16, offset=0.09, font.size=1.5, color=NA, colnames_angle = 45, colnames_offset_y = -0.75) +
  scale_fill_manual(values = c("1" = "#125c8a", "0" = "#c5d9ed", "NA" = "white"), na.value = "white", labels=c('Absent','Present','Gene absent'), name = "SNP Presence") +
  new_scale_fill()

p_snps_nonsyn = gheatmap(p_snps_syn, snps_nonsyn, width=0.08, offset=0.157, font.size=1.5, color=NA, colnames_angle = 45, colnames_offset_y = -0.75) +
  scale_fill_manual(values = c("1" = "#125c8a", "0" = "#c5d9ed", "NA" = "white"), na.value = "white", labels=c('Absent','Present','Gene absent'), name = "SNP Presence") +
  new_scale_fill()

p_tox = gheatmap(p_snps_nonsyn, cytotoxicity, width=0.03, offset=0.05, font.size=1.5, color=NA, colnames = FALSE) + 
  scale_fill_distiller(palette = "Reds", direction = 1, na.value = "white", "Cytotoxicity")  

ggsave("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Figures\\Biomarkers\\gene_snp_tree_4_24_25.png", p_tox, units = "in", width = 14, height = 11, dpi = 600)



########stacked bar chart of groups per snp

snps_df = full_join(nheA_snps, nheB_snps, by = c("Isolate", "Average.Cell.Viability.Reversed", "Adjusted.panC.Group..predicted.species.")) %>%
  full_join(nheC_snps, by = c("Isolate", "Average.Cell.Viability.Reversed", "Adjusted.panC.Group..predicted.species.")) %>%
  full_join(hblC_snps, by = c("Isolate", "Average.Cell.Viability.Reversed", "Adjusted.panC.Group..predicted.species.")) %>%
  full_join(hblD_snps, by = c("Isolate", "Average.Cell.Viability.Reversed", "Adjusted.panC.Group..predicted.species."))

snps_df <- snps_df[ -c(1, 11, 14, 21, 28) ]
test = data.frame(prop.table(table(snps_df$Adjusted.panC.Group..predicted.species., snps_df$nheC_SNP_34), margin = 2))
phylo = data.frame(prop.table(table(snps_df$Adjusted.panC.Group..predicted.species.)))
phylo$Var2 = 0

ggplot(phylo, aes(x = Var2, y = Freq, fill = Var1)) +
  geom_bar(position = "fill", stat = "identity", color = "black") +
  theme_prism() +
  theme_prism(axis_text_angle = 45) +
  scale_y_continuous(expand= c(0,0)) +
  xlab("Presence/absence") +
  ylab("Proportion of groups") 


#table of stats per group
group_by(df, Adjusted.panC.Group..predicted.species.) %>%
  summarise(number = n(), avg_cyt = mean(Average.Cell.Viability....0.7.is.cytotoxic.), sd = sd(Average.Cell.Viability....0.7.is.cytotoxic.))


#####How many isolates have two or more gene operons/genes?

genes2 = genes %>%
  mutate(nhe = recode(nhe, "TRUE" = 1, "FALSE" = 0)) %>%
  mutate(hbl = recode(hbl, "TRUE" = 1, "FALSE" = 0)) %>%
  mutate(cytk1 = recode(cytk1, "TRUE" = 1, "FALSE" = 0)) %>%
  mutate(cytk2 = recode(cytk2, "TRUE" = 1, "FALSE" = 0))

test = data.frame(gene_count = rowSums(genes2))

test %>%
  filter(gene_count >= 2) %>%
  summarize(n = n())
