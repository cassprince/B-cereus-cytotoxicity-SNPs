library(tidyverse)
library(readxl)
library(ivs)
library(ggbreak)


setwd("C:/Users/cassp/OneDrive/Documents/Kovac Lab/Biomarkers paper/genes_nano_031826")

columns = c("query", "acc", "perc_id", "ali_len", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "eval", "score")

lengths = read_tsv("contig_lengths_nano.txt", col_names = c("acc", "contig_length")) %>%
  separate(acc, c("acc", "contig"), "_")

plasflow = read_excel("Plasflow_Cassidy_031926.xlsx") %>%
  select(acc = Isolate, contig_length, label)

df_nheA = read_tsv("nheA_blast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_") %>%
  left_join(lengths, by = join_by(acc, contig))

df_nheB = read_tsv("nheB_blast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_") %>%
  left_join(lengths, by = join_by(acc, contig))

df_nheC = read_tsv("nheC_blast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_") %>%
  left_join(lengths, by = join_by(acc, contig))

df_cytK1 = read_tsv("cytK1_blast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_") %>%
  left_join(lengths, by = join_by(acc, contig))

df_cytK2 = read_tsv("cytK2_blast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_") %>%
  left_join(lengths, by = join_by(acc, contig))

df_hblA = read_tsv("hblA_blast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_")  %>%
  filter(query == "hblA_NZ_JBJJSR010000001.1:c3064164-3063037")

df_hblB = read_tsv("hblB_blast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_")  %>%
  filter(query == "hblB_NZ_JBJJSR010000001.1:c3062661-3061261")

df_hblC = read_tsv("hblC_blast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_")  %>%
  left_join(lengths, by = join_by(acc, contig)) %>%
  filter(query == "hblC_NZ_JBJJSR010000001.1:3065483-3066802")

df_hblD = read_tsv("hblD_blast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_")  %>%
  left_join(lengths, by = join_by(acc, contig)) %>%
  filter(query == "hblD_NZ_JBJJSR010000001.1:c3065421-3064201")


df = data.frame(read_excel("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Biomarkers paper\\Mastersheet_no_clones.xlsx")) %>%
  select(acc = Isolate, N50, cytotoxicity = Average.Cell.Viability....0.7.is.cytotoxic., panC_group = Adjusted.panC.Group..predicted.species., iso_source = Isolation.Source, iso_source_spec = Isolation.Source.Specific) 


####################


df_hblA_mult_mega = df_hblA %>%
  drop_na() %>%
  mutate(full_name = paste0(acc, "_", contig))%>%
  mutate(start = ifelse(s_start < s_end, s_start,s_end)) %>%
  mutate(stop = ifelse(s_start > s_end, s_start,s_end))

df_hblA_mega = df_hblA_mult_mega %>% 
  filter(perc_id > 85) %>% 
  distinct(acc, contig, s_start, .keep_all = TRUE)


df_hblB_mult_mega = df_hblB %>%
  drop_na() %>%
  mutate(full_name = paste0(acc, "_", contig))%>%
  mutate(start = ifelse(s_start < s_end, s_start,s_end)) %>%
  mutate(stop = ifelse(s_start > s_end, s_start,s_end))

df_hblB_mega = df_hblB_mult_mega %>% 
  filter(perc_id > 85) %>% 
  distinct(acc, contig, s_start, .keep_all = TRUE)


# Distinguish whether a shared gene hit is more like hblA or hblB and label it as such.
joined = full_join(df_hblA_mega, df_hblB_mega, by = join_by("acc", "contig"), suffix = c("_hblA", "_hblB"))

check = joined %>%
  mutate(overlap = map2_lgl(iv(start_hblA, stop_hblA), iv(start_hblB, stop_hblB), type = "any", iv_overlaps)) %>%
  mutate(best_score = ifelse(overlap == "TRUE", pmax(perc_id_hblA, perc_id_hblB), NA)) %>%
  mutate(best_score = ifelse(overlap == "TRUE", pmax(perc_id_hblA, perc_id_hblB), NA))

check$best_score_gene = ifelse(check$overlap == TRUE, names(check[c("perc_id_hblA", "perc_id_hblB")])[max.col(check[c("perc_id_hblA", "perc_id_hblB")], "first")], NA)
check$best_score_gene = sub("perc_id_", "", check$best_score_gene)


df_hblA_final = check %>%
  filter(best_score_gene == "hblA" | is.na(perc_id_hblB) | (overlap == FALSE & !is.na(perc_id_hblA))) %>% 
  arrange(desc(perc_id_hblA)) %>%
  distinct(acc, contig, start_hblA, .keep_all = TRUE)%>% ### Filtered out any genes with the same genomic start location as another hit. Helped limit the number same hits from different 
  distinct(acc, contig, stop_hblA, .keep_all = TRUE) %>% ### Filtered out any genes with the same genomic stop location
  select(acc, contig, ends_with("hblA")) %>%
  rename_with(~str_remove(., '_hblA')) %>%
  left_join(lengths, by = join_by(acc, contig))


df_hblB_final = check %>%
  filter(best_score_gene == "hblB" | is.na(perc_id_hblA) | (overlap == FALSE & !is.na(perc_id_hblB))) %>% 
  arrange(desc(perc_id_hblB)) %>%
  distinct(acc, contig, start_hblB, .keep_all = TRUE)%>% ### Filtered out any genes with the same genomic start location as another hit. Helped limit the number same hits from different 
  distinct(acc, contig, stop_hblB, .keep_all = TRUE) %>% ### Filtered out any genes with the same genomic stop location
  select(acc, contig, ends_with("hblB")) %>%
  rename_with(~str_remove(., '_hblB')) %>%
  left_join(lengths, by = join_by(acc, contig))

##########

df_all_genes = lengths %>%
  distinct(acc) %>%
  left_join(df_nheA %>% group_by(acc) %>% summarize(nheA = n()), 
            by = join_by(acc)) %>%
  left_join(df_nheB %>% group_by(acc) %>% summarize(nheB = n()), 
            by = join_by(acc)) %>%
  left_join(df_nheC %>% group_by(acc) %>% summarize(nheC = n()), 
            by = join_by(acc)) %>%
  left_join(df_hblC %>% group_by(acc) %>% summarize(hblC = n()), 
            by = join_by(acc)) %>%
  left_join(df_hblD %>% group_by(acc) %>% summarize(hblD = n()), 
            by = join_by(acc)) %>%
  left_join(df_hblA_final %>% group_by(acc) %>% summarize(hblA = n()), 
            by = join_by(acc)) %>%
  left_join(df_hblB_final %>% group_by(acc) %>% summarize(hblB = n()), 
            by = join_by(acc)) %>%
  left_join(df_cytK1 %>% group_by(acc) %>% summarize(cytK1 = n()), 
            by = join_by(acc)) %>%
  left_join(df_cytK2 %>% group_by(acc) %>% summarize(cytK2 = n()), 
            by = join_by(acc)) %>%
  replace(is.na(.), 0) %>%
  left_join(df, by = join_by(acc))

df_all_genes_long  = df_all_genes %>%
  pivot_longer(hblC:hblB, names_to = "gene", values_to = "count")


df_plasflow = plasflow %>%
  left_join(df_nheA %>% group_by(acc, contig_length) %>% summarize(nheA = n()), 
            by = join_by(acc, contig_length)) %>%
  left_join(df_nheB %>% group_by(acc, contig_length) %>% summarize(nheB = n()), 
            by = join_by(acc, contig_length)) %>%
  left_join(df_nheC %>% group_by(acc, contig_length) %>% summarize(nheC = n()), 
            by = join_by(acc, contig_length)) %>%
  left_join(df_hblC %>% group_by(acc, contig_length) %>% summarize(hblC = n()), 
            by = join_by(acc, contig_length))%>%
  left_join(df_hblD %>% group_by(acc, contig_length) %>% summarize(hblD = n()), 
            by = join_by(acc, contig_length)) %>%
  left_join(df_hblA_final %>% group_by(acc, contig_length) %>% summarize(hblA = n()), 
            by = join_by(acc, contig_length)) %>%
  left_join(df_hblB_final %>% group_by(acc, contig_length) %>% summarize(hblB = n()), 
            by = join_by(acc, contig_length)) %>%
  left_join(df_cytK1 %>% group_by(acc, contig_length) %>% summarize(cytK1 = n()), 
            by = join_by(acc, contig_length)) %>%
  left_join(df_cytK2 %>% group_by(acc, contig_length) %>% summarize(cytK2 = n()), 
            by = join_by(acc, contig_length)) %>%
  replace(is.na(.), 0)


######


df_plasflow %>%
  group_by(label) %>%
  summarize(n = n(), nheA = sum(nheA), nheB = sum(nheB), nheC = sum(nheC), hblC = sum(hblC), hblD = sum(hblD), hblA = sum(hblA), hblB = sum(hblB), cytK1 = sum(cytK1), cytK2 = sum(cytK2))

df_plasflow %>% 
  distinct(acc) %>%
  summarize(n())


df_all_genes %>%
  group_by(iso_source) %>%
  summarize(n = n(), nheA = sum(nheA), nheB = sum(nheB), nheC = sum(nheC), hblC = sum(hblC), hblD = sum(hblD), hblA = sum(hblA), hblB = sum(hblB), cytK1 = sum(cytK1), cytK2 = sum(cytK2))

df_all_genes %>%
  filter(hblD>1) %>%
  group_by(iso_source) %>%
  summarize(n = n())

sumD = df_hblD %>%
  group_by(acc) %>%
  summarize(total_hblD_hits = n()) 

df_hblD %>%
  group_by(acc) %>%
  mutate(hblD_copy_number = as.character(n())) %>%
  ggplot(aes(x = perc_id, y = contig_length, color = hblD_copy_number)) +
  geom_point(alpha = 0.5) +
  scale_y_log10(n.breaks = 25, labels = scales::label_number()) +
  scale_x_break(c(85, 95), scales = 0.7) +
  labs(x = "% identity to B. cereus s.s. hblD", y = "contig length (bp)")

#ggsave("hblD_copy_contig_031826.png", height = 5, width = 8, units = "in")


#write.csv(sumD, file = "df_hblD_031826.csv", row.names = FALSE)


test = lengths %>% 
  full_join(plasflow, by = join_by(acc, contig_length)) %>%
  distinct(contig, .keep_all = TRUE) %>%
  mutate(full_name = paste0(acc, "_", contig)) %>%
  filter(grepl("plasmid", label))
  
#write.table(test$full_name, file = "plasmid_contigs.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
