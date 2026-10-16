library(tidyverse)
library(readxl)

setwd("C:/Users/cassp/OneDrive/Documents/GitHub/B-cereus-cytotoxicity-SNPs/data/mafft")

files = list.files(path="C:/Users/cassp/OneDrive/Documents/GitHub/B-cereus-cytotoxicity-SNPs/data/mafft", pattern="*.csv", full.names=FALSE, recursive=FALSE)

files_lr_unfilt = files[grepl("unfilt_logreg",files)]

metadata = data.frame(read_excel("C:/Users/cassp/OneDrive - Cornell University/Biomarkers paper/Mastersheet_082026.xlsx")) %>%
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

level_order = c("nheA", "nheB", "nheC", "hblC", "hblD", "hblA", "hblB", "cytK1", "cytK2", "nheA_up", "hblC_up", "cytK2_up")

df_lr_one_unfilt %>%
  ggplot(aes(x = position, y = -log10(`lr_p-val_tox`))) +
  geom_hline(yintercept = 1.3, color = "gray70", linetype = 2) +
  geom_hline(yintercept = -log10(0.05/dim(df_lr_unfilt)[1]), color = "gray30", linetype = 1) +
  geom_point(aes(color = med_tox_diff, shape = phylo_signif)) +
  theme_classic() +
  labs(x = "Position", y = "-log10(p-value)") +
  theme(
    axis.text.x = element_blank()
  ) +
  scale_y_continuous(expand = c(0,0), n.breaks = 7, limits = c(0,15)) +
  scale_color_gradient(name = "difference in \nmedian cytotoxicity", low = "gray90", high = "darkblue", limits = c(0,0.7)) +
  scale_shape_manual(name = "significantly correlated\nwith panC group?", values = c(16, 17)) +
  facet_wrap(~factor(gene, levels = level_order), nrow = 1, scales = "free_x") 



df_plot = df_lr_one_unfilt %>%
  mutate(top_SNP = df_lr_one_unfilt$name %in% df_final$name) 

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
  facet_wrap(~factor(gene, levels = level_order), nrow = 1, scales = "free_x") 


df_final %>%
  select(-`lr_p-val_tox_bonf`, -`lr_p-val_phylo_bonf`) %>%
  rename(identity = nuc) %>%
  write_csv("C:/Users/cassp/OneDrive/Documents/GitHub/B-cereus-cytotoxicity-SNPs/data/top_SNP_set_101526.csv")

df_plot %>%
  select(-`lr_p-val_tox_bonf`, -`lr_p-val_phylo_bonf`) %>%
  rename(identity = nuc) %>%
  write_csv("C:/Users/cassp/OneDrive/Documents/GitHub/B-cereus-cytotoxicity-SNPs/data/full_SNP_set_101526.csv")
  

