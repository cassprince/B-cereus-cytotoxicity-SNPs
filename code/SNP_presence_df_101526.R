library(tidyverse)

setwd("C:/Users/princecr/Documents/GitHub/B-cereus-cytotoxicity-SNPs/data/mafft")

files = list.files(path="C:/Users/princecr/Documents/GitHub/B-cereus-cytotoxicity-SNPs/data/mafft", pattern="*.csv", full.names=FALSE, recursive=FALSE)

files_snps = files[grepl("unfilt_snps",files)]

df_SNPs = read_csv("C:/Users/princecr/Documents/GitHub/B-cereus-cytotoxicity-SNPs/data/full_SNP_set_101526.csv") %>%
  filter(tox_signif == TRUE)

metadata = read_csv("C:/Users/princecr/Documents/GitHub/B-cereus-cytotoxicity-SNPs/Mastersheet_082026.csv") %>%
  select(Isolate, cytotoxicity = "Cytotoxicity ( >0.7 is cytotoxic)", panC_group = "Adjusted_panC_Group(predicted_species)") %>%
  mutate(panC_group = gsub("\\s*\\([^\\)]+\\)","",panC_group)) %>%
  mutate(panC_group = gsub("\\*", "", panC_group))

list = append(df_SNPs$name, c("Isolate"))
df = metadata
for (i in files_snps){
  print(i)
  
  df_raw_long = read_csv(i) %>%
    rename_with(~str_replace(., 'prot_', '')) %>%
    rename_with(~str_replace(., '600bp', '')) 
  
  
  
  df_raw = df_raw_long[, which((names(df_raw_long) %in% list)==TRUE)]
  
  if(length(df_raw)>= 2){
    df_sum = df_raw %>%
      group_by(Isolate) %>%
      summarise(across(colnames(df_raw)[2]:colnames(df_raw)[length(df_raw)], ~ sum(.x, na.rm = TRUE)))
    
    df = df_sum %>% 
      right_join(df, by = join_by(Isolate)) %>%
      drop_na(Isolate) 
  }
   
}

# In file below, NA = the genome doesn't have the gene with that SNP, 0 = the gene is present but none of the genes have the SNP, 1 = one gene in the genome has the SNP, 2 = two genes have the SNP, etc. 
df %>%
  relocate(Isolate, cytotoxicity, panC_group) %>%
  write_csv("C:/Users/princecr/Documents/GitHub/B-cereus-cytotoxicity-SNPs/data/snps_all_signif_summed_102126.csv")

# In file below, NA = the genome doesn't have the gene with that SNP, 0 = the gene is present but none of the genes have the SNP, 1 = AT LEAST ONE gene in the genome has the SNP.
df %>%
  relocate(Isolate, cytotoxicity, panC_group) %>%
  mutate(across(colnames(.)[4]:colnames(.)[length(.)], ~ gsub(2, 1, .))) %>%
  mutate(across(colnames(.)[4]:colnames(.)[length(.)], ~ gsub(3, 1, .))) %>%
  write_csv("C:/Users/princecr/Documents/GitHub/B-cereus-cytotoxicity-SNPs/data/snps_all_signif_pres_abs_102126.csv")



