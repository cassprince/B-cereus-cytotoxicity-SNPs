library(readxl)
library(tidyverse)

setwd("C:/Users/princecr/Documents/GitHub/B-cereus-cytotoxicity-SNPs/data")
df = data.frame(read_excel("Mastersheet_082026.xlsx"))
snps = data.frame(read_csv("top_snps_111326.csv")) 

df_genes = read_csv("df_all_genes_090126.csv") %>%
  select(-cytotoxicity) %>%
  relocate(acc, panC_group, num_contigs, N50, N90)


final_df = df %>% 
  select(isolate = Isolate, ncbi_accession = SRR.Number.or.Assembly, genomospecies = final_taxon_names, isolation_source = Isolation.Source, isolation_year = Isolation.Year, normalized_cytotoxicity = Cytotoxicity....0.7.is.cytotoxic.) %>%
  left_join(df_genes, by = join_by(isolate == acc)) %>%
mutate(ncbi_accession = as.character(gsub("\\*", "", df$SRR.Number.or.Assembly))) %>% 
  left_join(snps, by = join_by(isolate == Isolate))


write.csv(final_df, "Table_S1_111326.csv", row.names = FALSE)
