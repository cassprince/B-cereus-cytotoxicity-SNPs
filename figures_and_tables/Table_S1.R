library(readxl)
library(tidyverse)

setwd("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Biomarkers paper")
df = data.frame(read_excel("Mastersheet_no_clones.xlsx"))
snps = data.frame(read_excel("SNP hits\\Covars 8_14_23\\SNP_hits_sheet.xlsx"))
snp_names = data.frame(read.csv("SNP hits\\Covars 8_14_23\\SNP_names.csv")) %>%
  select(old_names, new_names)


final_df = df %>% 
  select(Isolate, Alt.Name, NCBI.Accenssion.Number...Assembly.only., Adjusted.panC.Group..predicted.species., final.taxon.names, Isolation.Source, Isolation.Year, Average.Cell.Viability....0.7.is.cytotoxic.) %>%
  mutate(nheA = as.character(grepl("nheA", df$vir.nhe))) %>% 
  mutate(nheB = as.character(grepl("nheB", df$vir.nhe))) %>% 
  mutate(nheC = as.character(grepl("nheC", df$vir.nhe))) %>% 
  mutate(nheABC = as.character(grepl("3/3", df$vir.nhe))) %>% 
  mutate(hblA = as.character(grepl("hblA", df$vir.hbl))) %>% 
  mutate(hblB = as.character(grepl("hblB", df$vir.hbl))) %>% 
  mutate(hblC = as.character(grepl("hblC", df$vir.hbl))) %>% 
  mutate(hblD = as.character(grepl("hblD", df$vir.hbl))) %>% 
  mutate(hblACD = as.character(grepl("hblA.*hblC.*hblD", df$vir.hbl))) %>%
  mutate(cytK1 = as.character(grepl("cytK-1", df$vir.cytK))) %>% 
  mutate(cytK2 = as.character(grepl("cytK-2", df$vir.cytK))) %>%
  mutate(NCBI.Accenssion.Number...Assembly.only. = as.character(gsub("\\*", "", df$NCBI.Accenssion.Number...Assembly.only.))) %>%
  rename(NCBI_accession = NCBI.Accenssion.Number...Assembly.only.) %>%
  rename(alternate_name = Alt.Name) %>%
  rename(isolate = Isolate) %>%
  rename(normalized_cytotoxicity = Average.Cell.Viability....0.7.is.cytotoxic.) %>%
  rename(isolation_source = Isolation.Source) %>%
  rename(isolation_year = Isolation.Year) %>%
  rename(genomospecies = final.taxon.names) %>%
  rename(panC_group = Adjusted.panC.Group..predicted.species.) %>%
  left_join(snps[,c(1,4:24)], join_by(isolate == Isolate))

final_df = data.frame(lapply(final_df, function(x) {gsub("TRUE", "1", x)}))
final_df = data.frame(lapply(final_df, function(x) {gsub("FALSE", "0", x)}))
colnames(final_df)[20:40] = snp_names$new_names

write.csv(final_df, "Table_S1.csv", row.names = FALSE)
