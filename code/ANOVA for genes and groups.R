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


df_plot = df_genes %>%
  pivot_longer(cols = c(nheABC, hblCDAB, cytK1, cytK2), names_to = "genes", values_to = "presence") %>%
  mutate(presence = gsub("1", "present", presence)) %>%
  mutate(presence = gsub("0", "absent", presence))%>%
  mutate(genes = factor(genes, levels = c("nheABC", "hblCDAB", "cytK1", "cytK2"))) 
df_summ = df_plot %>%
  group_by(genes, presence) %>%
  summarize(n = n())


ggplot(df_plot, aes(x = presence, y = cytotoxicity)) +
  geom_hline(yintercept = 0.51, color = "gray", linetype = 1) +
  stat_halfeye(justification = -0.9, width = 0.3, .width = 0) +
  geom_boxplot(width = 0.4, aes(fill = presence)) +
  geom_signif(comparisons = list(c("present", "absent")), y_position = 1.35) +
  geom_text(data = df_summ, aes(x = presence, y = 1.65, label = paste0("n=", n)), size = 4) +
  theme_classic() +
  theme(text = element_text(color = "black", size = 14), 
        axis.text = element_text(color = "black", size = 14), 
        axis.title.y = element_text(color = "black", size = 15), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        strip.text = element_text(face = "italic")) +
  labs(x = "", y = "Normalized Cytotoxicity") +
  scale_fill_manual(values = c("#FCC983", "#B5DCFF")) +
  facet_wrap(~genes, nrow = 1, )
  
ggsave("../fig3_raincloud_genes_111126.png", width = 7.4, height = 4.3, units = "in")


#scale_color_manual(values = c("gray", "#FC8383", "#FCC983", "#FCF083", "#B3FC83", "#B5DCFF", "#A0A5F7", "#C1A0F7", "#F7A0E7")) +


