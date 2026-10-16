library(tidyverse)
library(readxl)

setwd("C:/Users/cassp/OneDrive/Documents/GitHub/B-cereus-cytotoxicity-SNPs/data/mafft")

files = list.files(path="C:/Users/cassp/OneDrive/Documents/GitHub/B-cereus-cytotoxicity-SNPs/data/mafft", pattern="*.csv", full.names=FALSE, recursive=FALSE)

files_snps = files[grepl("snps",files)]

df_SNPs = read_csv("C:/Users/cassp/OneDrive/Documents/GitHub/B-cereus-cytotoxicity-SNPs/data/full_SNP_set_101526.csv") %>%
  filter(top_SNP == TRUE)

metadata = data.frame(read_excel("C:/Users/cassp/OneDrive - Cornell University/Biomarkers paper/Mastersheet_082026.xlsx")) %>%
  select(Isolate, cytotoxicity = Cytotoxicity....0.7.is.cytotoxic., panC_group = Adjusted_panC_Group.predicted_species.) %>%
  mutate(panC_group = gsub("\\s*\\([^\\)]+\\)","",panC_group)) %>%
  mutate(panC_group = gsub("\\*", "", panC_group))

for (i in files_snps){  
  name = paste0(gsub("_snps_101426.csv","",i), "_snps")  
  assign(name, read.csv(i, header=TRUE))
  
}




df = metadata
for (i in files_snps){
  test = read_csv(i)
  
  test %>% filter(colnames)
  
}


list = append(df_SNPs$name, c("Isolate", "cytotoxicity"))
df = metadata
for (i in files_snps){
  if (grepl("nhe", i)){
    df_raw = read_csv(i) %>%
    rename_with(~str_replace(., 'prot_', '')) %>%
    rename_with(~str_replace(., '600bp', '')) 
    
    newDF = df_raw[, which((names(df_raw) %in% list)==TRUE)]
    df = newDF %>% right_join(df, by = join_by(Isolate, cytotoxicity)) }}




