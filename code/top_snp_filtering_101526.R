library(tidyverse)
library(readxl)

setwd("C:/Users/princecr/Documents/GitHub/B-cereus-cytotoxicity-SNPs/data/mafft")

files = list.files(path="C:/Users/princecr/Documents/GitHub/B-cereus-cytotoxicity-SNPs/data/mafft", pattern="*.csv", full.names=FALSE, recursive=FALSE)

files_lr_unfilt = files[grepl("unfilt_logreg",files)]

metadata = data.frame(read_excel("../Mastersheet_082026.xlsx")) %>%
  select(Isolate, cytotoxicity = Cytotoxicity....0.7.is.cytotoxic., panC_group = Adjusted_panC_Group.predicted_species.) %>%
  mutate(panC_group = gsub("\\s*\\([^\\)]+\\)","",panC_group)) %>%
  mutate(panC_group = gsub("\\*", "", panC_group))

df_lr_unfilt = data.frame()  
for (x in files_lr_unfilt) {
  df_lr_unfilt = df_lr_unfilt %>%
    rbind(read_csv(x))
}

df_lr_unfilt = df_lr_unfilt %>%
  mutate(lr_p_val_tox_bonf = p.adjust(`lr_p-val_tox`, method = "bonferroni")) %>%
  mutate(lr_p_val_phylo_bonf = p.adjust(`lr_p-val_phylo`, method = "bonferroni"))

df_lr_one_unfilt = df_lr_unfilt %>%
  arrange(lr_p_val_tox_bonf) %>%
  distinct(gene, position, .keep_all = TRUE) %>%
  mutate(phylo_signif = lr_p_val_phylo_bonf < 0.05) %>%
  mutate(tox_signif = lr_p_val_tox_bonf < 0.05) %>%
  mutate(avg_tox_diff = abs(avg_tox_w - avg_tox_wo)) %>%
  mutate(med_tox_diff = abs(median_tox_w - median_tox_wo)) %>%
  mutate(total_seqs = total_w+total_wo) %>%
  mutate(gene = gsub("_prot", "", gene)) %>%
  mutate(gene = gsub("600bp", "", gene)) %>%
  mutate(name = gsub("_prot", "", name)) %>%
  mutate(name = gsub("600bp", "", name)) 

df_prefinal_filt = df_lr_one_unfilt %>%
  filter(tox_signif == TRUE) 


df_low_p = df_prefinal_filt %>%
  slice_min(lr_p_val_tox_bonf, n = 80)

df_high_tox = df_prefinal_filt %>%
  slice_max(med_tox_diff, n = 80)

df_final = df_low_p %>%
  inner_join(df_high_tox) %>%
  rbind(df_prefinal_filt %>%
          slice_max(med_tox_diff, n = 5)) %>%
  rbind(df_prefinal_filt %>%
          slice_min(lr_p_val_tox_bonf, n = 5)) %>%
  distinct(.keep_all = TRUE)

level_order = c("nheA up", "NheA", "NheB", "NheC", "hblC up", "HblC", "HblD", "HblA", "HblB", "CytK1", "cytK2 up", "CytK2")

df_lr_one_unfilt %>%
  mutate(gene = gsub("nheA", "NheA", gene)) %>%
  mutate(gene = gsub("nheB", "NheB", gene)) %>%
  mutate(gene = gsub("nheC", "NheC", gene)) %>%
  mutate(gene = gsub("NheA_up", "nheA up", gene)) %>%
  mutate(gene = gsub("hblA", "HblA", gene)) %>%
  mutate(gene = gsub("hblB", "HblB", gene)) %>%
  mutate(gene = gsub("hblC", "HblC", gene)) %>%
  mutate(gene = gsub("hblD", "HblD", gene)) %>%
  mutate(gene = gsub("HblC_up", "hblC up", gene)) %>%
  mutate(gene = gsub("cytK1", "CytK1", gene)) %>%
  mutate(gene = gsub("cytK2", "CytK2", gene)) %>%
  mutate(gene = gsub("CytK2_up", "cytK2 up", gene)) %>%
  ggplot(aes(x = position, y = -log10(`lr_p-val_tox`))) +
  geom_hline(yintercept = 1.3, color = "gray70", linetype = 2) +
  geom_hline(yintercept = -log10(0.05/dim(df_lr_unfilt)[1]), color = "gray30", linetype = 1) +
  geom_point(aes(color = avg_tox_diff, shape = phylo_signif)) +
  theme_classic() +
  labs(x = "Position", y = expression(-log["10"]("p-value"))) +
  theme(text = element_text(color = "black", size = 13),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14)) +
  scale_y_continuous(expand = c(0,0), n.breaks = 7, limits = c(0,15)) +
  scale_color_gradient(name = "Difference in \nmean cytotoxicity", low = "gray90", high = "navy", limits = c(0,0.7)) +
  scale_shape_manual(name = "Significantly correlated\nwith panC group?", values = c(16, 17)) +
  facet_wrap(~factor(gene, levels = level_order), nrow = 1, scales = "free_x") 


ggsave("../../figures/fig4_manhattan_110426.png", width = 12, height = 4, units = "in", dpi = 350)








df_plot = df_lr_one_unfilt %>%
  mutate(top_SNP = (name=="hblB_SNP_T_270" | name=="hblA_SNP_E_210"| name=="hblC_up_SNP_g_561"| name=="nheA_up_SNP_t_312"| name=="hblC_up_SNP_t_355"| name=="nheA_up_SNP_c_418"))

df_plot %>%
  ggplot(aes(x = position, y = -log10(`lr_p-val_tox`))) +
  geom_hline(yintercept = -log10(0.05), color = "gray70", linetype = 2) +
  geom_hline(yintercept = -log10(0.05/dim(df_lr_unfilt)[1]), color = "gray30", linetype = 1) +
  geom_point(data = df_plot %>% filter(top_SNP == "FALSE"),
             aes(color = top_SNP, shape = phylo_signif), alpha = 0.5) +
  geom_point(data = df_plot %>% filter(top_SNP == "TRUE"),
             aes(color = top_SNP, shape = phylo_signif)) +
  theme_classic() +
  labs(x = "Position", y = "-log10(p-value)") +
  theme(
    axis.text.x = element_blank()
  ) +
  scale_y_continuous(expand = c(0,0), n.breaks = 7, limits = c(0,15)) +
  scale_color_manual(values = c("TRUE" = "cornflowerblue", "FALSE" = "grey70")) +
  scale_shape_manual(name = "significantly correlated\nwith panC group?", values = c(16, 17)) +
  facet_wrap(~gene, nrow = 1, scales = "free_x") 


df_final %>%
  select(-`lr_p-val_tox_bonf`, -`lr_p-val_phylo_bonf`) %>%
  rename(identity = nuc) %>%
  write_csv("C:/Users/cassp/OneDrive/Documents/GitHub/B-cereus-cytotoxicity-SNPs/data/top_SNP_set_101526.csv")

df_plot %>%
  select(-`lr_p-val_tox_bonf`, -`lr_p-val_phylo_bonf`) %>%
  rename(identity = nuc) %>%
  write_csv("C:/Users/cassp/OneDrive/Documents/GitHub/B-cereus-cytotoxicity-SNPs/data/full_SNP_set_101526.csv")
  

