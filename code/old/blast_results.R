library(tidyverse)
library(readxl)
library(ivs)


#### new files are in genes_112525. work with those instead! Double check that I fixed the hblD/A/B query mess before that though... I DIDNT! RERUN hblD with correct query!!!!

### also run hblA and B mult one more time, then extract seqs 

columns = c("query", "acc", "perc_id", "ali_len", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "eval", "score")

df_hblB = read_tsv("C:/Users/cassp/OneDrive/Documents/Kovac Lab/Biomarkers paper/genes_112125/hblB_blast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_")

df_hblA = read_tsv("C:/Users/cassp/OneDrive/Documents/Kovac Lab/Biomarkers paper/genes_112125/hblA_blast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_")

df_hblC = read_tsv("C:/Users/cassp/OneDrive/Documents/Kovac Lab/Biomarkers paper/genes_112125/hblC_blast.csv", skip = 5, col_names = c("query", "acc", "perc_id", "ali_len", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "eval", "score")) %>%
  separate(acc, c("acc", "contig"), "_")

df_hblD = read_tsv("C:/Users/cassp/OneDrive/Documents/Kovac Lab/Biomarkers paper/genes_112125/hblD_blast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_")

df_cytK1 = read_tsv("C:/Users/cassp/OneDrive/Documents/Kovac Lab/Biomarkers paper/genes_112125/cytK1_blast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_") %>%
  filter(perc_id > 90)

df_cytK2 = read_tsv("C:/Users/cassp/OneDrive/Documents/Kovac Lab/Biomarkers paper/genes_112125/cytK2_blast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_") %>%
  filter(perc_id > 90)


sumD = df_hblD %>%
  group_by(acc) %>%
  summarize(n = n())





df_hblA_mega = read_tsv("C:/Users/cassp/OneDrive/Documents/Kovac Lab/Biomarkers paper/genes_112125/hblA_megablast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_") %>%
  mutate(start = ifelse(s_start < s_end, s_start,s_end)) %>%
  mutate(stop = ifelse(s_start > s_end, s_start,s_end))

df_hblB_mega = read_tsv("C:/Users/cassp/OneDrive/Documents/Kovac Lab/Biomarkers paper/genes_112125/hblB_megablast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_") %>%
  mutate(start = ifelse(s_start < s_end, s_start,s_end)) %>%
  mutate(stop = ifelse(s_start > s_end, s_start,s_end))


sumBMega = df_hblB_mega %>%
  group_by(acc) %>%
  summarize(n = n())

sumAMega = df_hblA_mega %>%
  group_by(acc) %>%
  summarize(n = n())


df = data.frame(read_excel("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Biomarkers paper\\SNP hits\\Covars 8_14_23\\SNP_hits_sheet.xlsx")) %>%
  select(Isolate, Average.Cell.Viability, Adjusted.panC.Group)

sumAMega = left_join(sumAMega, df, by = join_by(acc == Isolate))
sumBMega = left_join(sumBMega, df, by = join_by(acc == Isolate))

####################


df_hblA_mult_mega = read_tsv("C:/Users/cassp/OneDrive/Documents/Kovac Lab/Biomarkers paper/genes_112125/hblA_mult_megablast.csv", skip = 5, col_names = columns) %>%
  drop_na() %>%
  mutate(full_name = acc)%>%
  separate(query, c("gene", "spp"), "_") %>%
  separate(acc, c("acc", "contig"), "_") %>%
  mutate(start = ifelse(s_start < s_end, s_start,s_end)) %>%
  mutate(stop = ifelse(s_start > s_end, s_start,s_end))

df_hblA_mega = df_hblA_mult_mega %>% 
  filter(perc_id > 85) %>% 
  distinct(acc, contig, s_start, .keep_all = TRUE)


df_hblB_mult_mega = read_tsv("C:/Users/cassp/OneDrive/Documents/Kovac Lab/Biomarkers paper/genes_112125/hblB_mult_megablast.csv", skip = 5, col_names = columns) %>%
  drop_na() %>%
  mutate(full_name = acc)%>%
  separate(query, c("gene", "spp"), "_") %>%
  separate(acc, c("acc", "contig"), "_") %>%
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
  select(acc, ends_with("hblA")) %>%
  rename_with(~str_remove(., '_hblA'))


df_hblB_final = check %>%
  filter(best_score_gene == "hblB" | is.na(perc_id_hblA) | (overlap == FALSE & !is.na(perc_id_hblB))) %>% 
  arrange(desc(perc_id_hblB)) %>%
  distinct(acc, contig, start_hblB, .keep_all = TRUE)%>% ### Filtered out any genes with the same genomic start location as another hit. Helped limit the number same hits from different 
  distinct(acc, contig, stop_hblB, .keep_all = TRUE) %>% ### Filtered out any genes with the same genomic stop location
  select(acc, ends_with("hblB")) %>%
  rename_with(~str_remove(., '_hblB')) 



