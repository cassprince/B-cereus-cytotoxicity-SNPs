library(tidyverse)

setwd("C:/Users/princecr/Documents/GitHub/B-cereus-cytotoxicity-SNPs/data")

df = read_csv("full_SNP_set_101526.csv") %>%
  mutate(snp = paste0(gene, "_", position, identity)) 

final_df = df %>%
  select(snp, total_seqs, total_with = total_w, total_without = total_wo, mean_tox_with = avg_tox_w, sd_tox_with = sd_tox_w, mean_tox_without = avg_tox_wo, sd_tox_without = sd_tox_wo, lr_pval_tox_bonf = lr_p_val_tox_bonf, lr_pval_phylo_bonf = lr_p_val_phylo_bonf)

write_csv(final_df, "Table_S4_111326.csv")
