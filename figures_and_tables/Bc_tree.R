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
df = data.frame(read_excel("Mastersheet_no_clones.xlsx"))
df[,"Adjusted.panC.Group..predicted.species."] = gsub("[*].*$", "", df[,"Adjusted.panC.Group..predicted.species."])

snps = data.frame(read_excel("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Biomarkers paper\\SNP hits\\Covars 8_14_23\\SNP_hits_sheet.xlsx"))
snps_for_tree = snps[4:24]
rownames(snps_for_tree) = snps[,"Isolate"]

cytotoxicity = data.frame(df[,"Average.Cell.Viability....0.7.is.cytotoxic."])
clade = data.frame(clade = df[,"Adjusted.panC.Group..predicted.species."])

nhe = data.frame(nhe = as.character(grepl("3/3", df$vir.nhe)))
hbl = data.frame(hbl = as.character(grepl("hblA.*hblC.*hblD", df$vir.hbl)))
cytk1 = data.frame(cytk1 = as.character(grepl("cytK-1", df$vir.cytK)))
cytk2 = data.frame(cytk2 = as.character(grepl("cytK-2", df$vir.cytK)))
genes = data.frame(cbind(nhe = nhe, hbl = hbl, cytk1 = cytk1, cytk2 = cytk2))

rownames(genes) = df[,"Isolate"]
rownames(cytotoxicity) = df[,"Isolate"]
rownames(clade) = df[,"Isolate"]
clade2 = data.frame(cbind(label = rownames(clade), clade = clade))

tree = read.newick("core_SNPs_matrix.biomarker.contree")
tree$tip.label = gsub("_contigs.fasta", "", tree$tip.label)
tree$tip.label = gsub(".fasta", "", tree$tip.label)
tree_mid <- midpoint(tree)



###---phylogram rectangluar---###
p = ggtree(tree_mid, size = 0.5) + geom_tiplab(size=1) + geom_treescale(linesize=0.75, offset = 1.2, x=0, y=200, width=0.05, fontsize=3.75) +
  geom_nodepoint(aes(fill = as.numeric(label)), size = 1, shape = 21) + 
  scale_fill_gradient(low = "white", high = "black", name = "Bootstrap\npercentage") +
  new_scale_fill() 


p_phylo = gheatmap(p, clade, offset=0.1, width=0.015, font.size=1, colnames = FALSE, color=NA) +
  scale_fill_discrete(name = "group", na.value = "white") + 
  ylim(-0.75, NA) + 
  coord_cartesian(clip = "off")

p_genes = gheatmap(p, genes, width=0.05, offset=0.107, font.size=1.5, color=NA, colnames_angle = 45, colnames_offset_y = -0.75) +
  scale_fill_manual(values = c("TRUE" = "#125c8a", "FALSE" = "#c5d9ed"), na.value = "white", labels=c('Absent','Present'), name = "Gene Presence") +
  new_scale_fill() 

p_snps = gheatmap(p_genes, snps_for_tree, width=0.3, offset=0.14, font.size=1.5, color=NA, colnames_angle = 45, colnames_offset_y = -0.75) +
  scale_fill_manual(values = c("1" = "#125c8a", "0" = "#c5d9ed", "NA" = "white"), na.value = "white", labels=c('Absent','Present','Gene absent'), name = "SNP Presence") +
  new_scale_fill()

p_tox = gheatmap(p_snps, cytotoxicity, width=0.03, offset=0.08, font.size=1.5, color=NA, colnames = FALSE) + 
  scale_fill_distiller(palette = "Reds", direction = 1, na.value = "white", "Cytotoxicity")  
#scale_fill_viridis(option = "B", direction = -1)
#scale_fill_gradient(low = "white", high = "black", name = "Cytotoxicity")

p_tox

#"1" = "#532196", "0" = "#d2c3e6"

png(file = "C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Figures\\Biomarkers\\gene_snp_tree_9_13_24.png", units = "in", width = 14, height = 11, res = 600)
#svg(file = "C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Figures\\Biomarkers\\gene_tree.svg")

p_tox

dev.off()




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
