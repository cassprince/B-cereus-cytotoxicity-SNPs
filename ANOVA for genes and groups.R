library(emmeans)
library(multcomp)
library(readxl)
library(tidyverse)

setwd("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Biomarkers paper")
df<-data.frame(read_excel("Mastersheet_no_clones.xlsx"))

nhe = data.frame(nhe = as.character(grepl("3/3", df$vir.nhe)))
hbl = data.frame(hbl = as.character(grepl("hblA.*hblC.*hblD", df$vir.hbl)))
cytk1 = data.frame(cytk1 = as.character(grepl("cytK-1", df$vir.cytK)))
cytk2 = data.frame(cytk2 = as.character(grepl("cytK-2", df$vir.cytK)))

##Clean species names.
species<-df[,"final.taxon.names"]
species <- gsub(";.*", "", species)
species <- recode(species, "B. mosaicus subsp. cereus biovar Emeticus*" = "B. mosaicus subsp. cereus biovar Emeticus")

Group = df[,"Adjusted.panC.Group..predicted.species."]
Group = gsub("\\*", "", Group)
cytotoxicity = df[,"Average.Cell.Viability....0.7.is.cytotoxic."]

data = data.frame(cbind(Group, hbl, nhe, cytk1, cytk2, cytotoxicity))


### Hbl

model = lm(cytotoxicity ~ Group*hbl, data=data) 
summary(model)

hbl_group = emmeans(model, pairwise ~ hbl|Group) #compare gene presence within a group. 
df_hbl_group = summary(hbl_group)$contrasts

group_hbl = emmeans(model, pairwise ~ Group|hbl) #between group comparisons within gene presence 
df_group_hbl = summary(group_hbl)$contrasts %>%
  mutate(hbl, gene = "hbl") %>%
  rename(gene_presence = hbl)


### Nhe

model = lm(cytotoxicity ~ Group*nhe, data=data) 
summary(model)

nhe_group = emmeans(model, pairwise ~ nhe|Group) #compare gene presence within a group. 
df_nhe_group = summary(nhe_group)$contrasts

group_nhe = emmeans(model, pairwise ~ Group|nhe) #between group comparisons within gene presence 
df_group_nhe = summary(group_nhe)$contrasts %>%
  mutate(nhe, gene = "nhe") %>%
  rename(gene_presence = nhe)


### CytK-2

model = lm(cytotoxicity ~ Group*cytk2, data=data) 
summary(model)

cytk2_group = emmeans(model, pairwise ~ cytk2|Group) #compare gene presence within a group. 
df_cytk2_group = summary(cytk2_group)$contrasts

group_cytk2 = emmeans(model, pairwise ~ Group|cytk2) #between group comparisons within gene presence 
df_group_cytk2 = summary(group_cytk2)$contrasts %>%
  mutate(cytk2, gene = "cytK2") %>%
  rename(gene_presence = cytk2)
  
## Final dfs

group_gene = bind_rows(df_group_hbl, df_group_nhe, df_group_cytk2) #between group comparisons within gene presence 

write.csv(group_gene, "Table_S3_group_gene.csv", row.names = FALSE, col.names = TRUE, na = "NA")

gene_group = bind_rows(list(hbl = df_hbl_group, nhe = df_nhe_group, cytK2 = df_cytk2_group), .id = "gene") #compare gene presence within a group. 
  
write.csv(gene_group, "Table_S3_gene_group.csv", row.names = FALSE, col.names = TRUE, na = "NA")


##########


anova_group = aov(cytotoxicity ~ Group, df) #ANOVA of groups and cytotoxicity.
summary(anova_group)
df_anova_group = data.frame(TukeyHSD(anova_group)$Group) %>%
  rownames_to_column(var = "groups_compared")

write.csv(df_anova_group, "Table_S2.csv", row.names = FALSE)




##### ANOVA for all interaction effects

test = aov(lm(cytotoxicity ~ hbl + nhe + cytk1 + cytk2 + Group + hbl:nhe:cytk1:cytk2:Group, data=data))
test
df_anova_group = data.frame(TukeyHSD(test)$Group) %>%
  rownames_to_column(var = "groups_compared")

write.csv(anova(test), "Table_S2_multiway.csv", row.names = FALSE)
write.csv(df_anova_group, "Table_S3_groups.csv", row.names = FALSE)
