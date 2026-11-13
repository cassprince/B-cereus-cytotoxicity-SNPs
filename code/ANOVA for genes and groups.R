library(emmeans)
library(multcomp)
library(tidyverse)
library(ggdist)
library(ggsignif)

setwd("C:/Users/cassp/Documents/GitHub/B-cereus-cytotoxicity-SNPs/data")
df_genes = read_csv("df_all_genes_090126.csv") %>%
  mutate(nheABC = ifelse(nheA>0 & nheB>0 & nheC >0, "1", "0")) %>%
  mutate(hblCDAB = ifelse(hblC>0 & hblD>0 & hblA>0 & hblB >0, "1", "0")) %>%
  mutate(across(c(nheABC, hblCDAB, cytK1, cytK2), ~ as.factor(.)))


### Hbl

model = lm(cytotoxicity ~ panC_group*hblCDAB, data=df_genes) 
summary(model)

hbl_group = emmeans(model, pairwise ~ hblCDAB|panC_group) #compare gene presence within a group. 
df_hbl_group = summary(hbl_group)$contrasts

group_hbl = emmeans(model, pairwise ~ panC_group|hblCDAB) #between group comparisons within gene presence 
df_group_hbl = summary(group_hbl)$contrasts


### Nhe

model = lm(cytotoxicity ~ panC_group*nheABC, data=df_genes) 
summary(model)

nhe_group = emmeans(model, pairwise ~ nheABC|panC_group) #compare gene presence within a group. 
df_nhe_group = summary(nhe_group)$contrasts

group_nhe = emmeans(model, pairwise ~ panC_group|nheABC) #between group comparisons within gene presence 
df_group_nhe = summary(group_nhe)$contrasts 

### CytK-2

model = lm(cytotoxicity ~ panC_group*cytK2, data=df_genes) 
summary(model)

cytk2_group = emmeans(model, pairwise ~ cytK2|panC_group) #compare gene presence within a group. 
df_cytk2_group = summary(cytk2_group)$contrasts

group_cytk2 = emmeans(model, pairwise ~ panC_group|cytK2) #between group comparisons within gene presence 
df_group_cytk2 = summary(group_cytk2)$contrasts
  
## Final dfs

group_gene = bind_rows(df_group_hbl, df_group_nhe, df_group_cytk2) #between group comparisons within gene presence 

write.csv(group_gene, "Table_S3_group_gene_111126.csv", row.names = FALSE, na = "NA")

gene_group = bind_rows(list(hblCDAB = df_hbl_group, nheABC = df_nhe_group, cytK2 = df_cytk2_group), .id = "gene") #compare gene presence within a group. 
  
write.csv(gene_group, "Table_S3_gene_group_111126.csv", row.names = FALSE,  na = "NA")


##########


anova_group = aov(cytotoxicity ~ panC_group, df_genes) #ANOVA of groups and cytotoxicity.
summary(anova_group)
df_anova_group = data.frame(TukeyHSD(anova_group)$panC_group) %>%
  rownames_to_column(var = "groups_compared")

write.csv(df_anova_group, "Table_S2_111126.csv", row.names = FALSE)




##### ANOVA for all interaction effects

test = aov(lm(cytotoxicity ~ hblCDAB + nheABC + cytK1 + cytK2 + panC_group + hblCDAB:nheABC:cytK1:cytK2:panC_group, data= df_genes))
test
summary(test)

df_anova_group = data.frame(TukeyHSD(test)$panC_group) %>%
  rownames_to_column(var = "groups_compared")

write.csv(anova(test), "Table_S2_multiway_111126.csv", row.names = FALSE)
write.csv(df_anova_group, "Table_S3_groups_111126.csv", row.names = FALSE)





