library(tidyverse)
library(Biostrings)
library(seqinr)

setwd("C:/Users/cassp/OneDrive/Documents/GitHub/B-cereus-cytotoxicity-SNPs/blast_results/filtered_megablast_qc50_090126/fasta_files")

# Geneious counts "-" as a mismatch which lowers the % id strongly. Instead, I calculate the % id using the # matches/sequence length.

hblC_seqs = readAAStringSet("translated_seqs/hblC_all_prot.fasta") %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename_with(~"seq", .cols = last_col()) %>%
  mutate(rowname = gsub(" <unknown description>", "", rowname)) %>%
  mutate(full_name = str_trim(rowname)) %>%
  separate(rowname, c("contig", "range"), ":") %>%
  separate(range, c("start", "stop"), "-") %>%
  separate(contig, c("acc", "contig"), "_") %>%
  mutate(start = as.numeric(start), stop = as.numeric(stop)) %>%
  mutate(seq_lengths = str_length(seq))


file_list = c("../perc_id_matrices/hblC_prot_PS03263_hblCDAB.csv", "../perc_id_matrices/hblC_prot_PS03263_hblCDA_ch.csv", "../perc_id_matrices/hblC_prot_PS03263_hblCDA_pl.csv")

df_prot = data.frame()

for (file in file_list){
  print(file)
  
  id_matrix_prot = read_csv(file) %>%
    rename(`...1` = "seq1") %>%
    column_to_rownames(var = "seq1")
  
  upper_idx_prot = upper.tri(id_matrix_prot)
  
  df_temp = data.frame(
    seq1 = rownames(id_matrix_prot)[row(id_matrix_prot)[upper_idx_prot]],
    seq2 = colnames(id_matrix_prot)[col(id_matrix_prot)[upper_idx_prot]],
    identity = id_matrix_prot[upper_idx_prot])
  
  df_prot = rbind(df_prot, df_temp)
}



# Compare all of the protein seqs to a spp that has 3 hblC variants. Variants 1 and 3 are similar to eachother. 

df_ch_prot = df_prot %>%
  filter(seq1 == "PS03263_assembly_1:1288801-1290118" | seq2 == "PS03263_assembly_1:1288801-1290118") %>%
  mutate(comp_seq = ifelse(seq1 == "PS03263_assembly_1:1288801-1290118", seq2, seq1)) %>%
  select(full_name = comp_seq, id_ch =identity) %>%
  left_join(hblC_seqs %>% select(full_name, seq_lengths), by = join_by(full_name)) %>%
  mutate(perc_id_ch = id_ch/seq_lengths*100)

df_pl_prot = df_prot %>%
  filter(seq1 == "PS03263_assembly_4:62756-64075" | seq2 == "PS03263_assembly_4:62756-64075") %>%
  mutate(comp_seq = ifelse(seq1 == "PS03263_assembly_4:62756-64075", seq2, seq1)) %>%
  select(full_name = comp_seq, id_pl = identity) %>%
  left_join(hblC_seqs %>% select(full_name, seq_lengths), by = join_by(full_name)) %>%
  mutate(perc_id_pl = id_pl/seq_lengths*100)

df_CDAB_prot = df_prot %>%
  filter(seq1 == "PS03263_assembly_1:1408214-1409531" | seq2 == "PS03263_assembly_1:1408214-1409531") %>%
  mutate(comp_seq = ifelse(seq1 == "PS03263_assembly_1:1408214-1409531", seq2, seq1)) %>%
  select(full_name = comp_seq, id_cdab = identity) %>%
  left_join(hblC_seqs %>% select(full_name, seq_lengths), by = join_by(full_name)) %>%
  mutate(perc_id_cdab = id_cdab/seq_lengths*100)
  
# Join the three dataframes.

df_join = left_join(df_ch_prot, df_pl_prot, by = join_by(full_name, seq_lengths)) %>%
  left_join(df_CDAB_prot, by = join_by(full_name, seq_lengths)) %>%
  drop_na() %>%
  select(-id_ch, -id_pl, -id_cdab, -seq_lengths) %>%
  pivot_longer(cols = c(perc_id_ch, perc_id_pl, perc_id_cdab), names_to = "variant", values_to = "perc_id") %>%
  group_by(full_name) %>%
  mutate(top_hit = (perc_id == max(perc_id))) %>%
  mutate(variant = gsub("perc_id_", "", variant))
  

df_join_filt = df_join %>%
  filter(top_hit == TRUE) %>%
  left_join(hblC_seqs %>% select(full_name, seq), by = join_by("full_name"))

# Write separate fasta files for the two main variants. 

write.fasta(df_join_filt %>% 
              filter(variant == "ch") %>% 
              pull(seq) %>% 
              as.list(), 
            df_join_filt %>% filter(variant == "ch") %>% pull(full_name) %>% as.list(), 
            "hblC_prot_ch.fasta",
            nbchar = 60)

write.fasta(df_join_filt %>% 
              filter(variant == "pl") %>% 
              pull(seq) %>% 
              as.list(), 
            df_join_filt %>% filter(variant == "pl") %>% pull(full_name) %>% as.list(), 
            "hblC_prot_pl.fasta",
            nbchar = 60)

write.fasta(df_join_filt %>% 
              filter(variant == "cdab") %>% 
              pull(seq) %>% 
              as.list(), 
            df_join_filt %>% filter(variant == "cdab") %>% pull(full_name) %>% as.list(), 
            "hblC_prot_CDAB.fasta",
            nbchar = 60)

##############

# Upload hblC upstream sequences and filter based on hblC variants from above.
# Start of hblC up is the stop of hblC

hblC_up_seqs = readDNAStringSet("upstream_seqs/hblC_up_600bp.fasta") %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename_with(~"seq", .cols = last_col()) %>%
  mutate(rowname = gsub(" <unknown description>", "", rowname)) %>%
  mutate(full_name = str_trim(rowname)) %>%
  separate(rowname, c("contig", "range"), ":") %>%
  separate(range, c("start", "stop"), "-") %>%
  separate(contig, c("acc", "contig"), "_") %>%
  mutate(start = as.numeric(start), stop = as.numeric(stop)) %>%
  mutate(start = ifelse(start < stop, start, stop)) %>%
  mutate(stop = ifelse(start > stop, start, stop)) %>%
  mutate(seq_lengths = str_length(seq))


df_up_filt = df_join_filt %>% 
  ungroup() %>%
  mutate(to_split = full_name) %>%
  separate(to_split, c("contig", "range"), ":") %>%
  separate(range, c("start", "stop"), "-") %>%
  separate(contig, c("acc", "contig"), "_") %>%
  mutate(start = as.numeric(start), stop = as.numeric(stop)) %>%
  mutate(start = ifelse(start < stop, start, stop)) %>%
  mutate(stop = ifelse(start > stop, start, stop)) %>%
  select(-seq, -full_name) %>%
  left_join(hblC_up_seqs, by = join_by(acc, contig, start == stop)) %>%
  drop_na() %>%
  select(-acc, -contig, -start.y, -stop, -start, -seq_lengths)


write.fasta(df_up_filt %>%
              filter(variant == "ch") %>%
              pull(seq) %>% 
              as.list(), 
            df_up_filt %>%
              filter(variant == "ch") %>% 
              pull(full_name) %>% 
              as.list(), 
            "hblC_up_ch.fasta",
            nbchar = 60)

write.fasta(df_up_filt %>%
              filter(variant == "pl") %>%
              pull(seq) %>% 
              as.list(), 
            df_up_filt %>%
              filter(variant == "pl") %>% 
              pull(full_name) %>% 
              as.list(), 
            "hblC_up_pl.fasta",
            nbchar = 60)

write.fasta(df_up_filt %>%
              filter(variant == "cdab") %>%
              pull(seq) %>% 
              as.list(), 
            df_up_filt %>%
              filter(variant == "cdab") %>% 
              pull(full_name) %>% 
              as.list(), 
            "hblC_up_CDAB.fasta",
            nbchar = 60)
