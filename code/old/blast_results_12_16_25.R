library(tidyverse)
library(readxl)
library(ivs)
library(ggridges)
library(Biostrings)
library(msa)

setwd("C:/Users/cassp/OneDrive/Documents/Kovac Lab/Biomarkers paper/genes_030326")

columns = c("query", "acc", "perc_id", "ali_len", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "eval", "score")

lengths = read_tsv("contigs_lengths.txt", col_names = c("acc", "contig_length")) %>%
  separate(acc, c("acc", "contig"), "_")

df_hblB = read_tsv("hblB_mult_megablast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_")

df_hblA = read_tsv("hblA_mult_megablast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_")

df_hblC = read_tsv("hblC_blast.csv", skip = 5, col_names = c("query", "acc", "perc_id", "ali_len", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "eval", "score")) %>%
  separate(acc, c("acc", "contig"), "_") %>%
  left_join(lengths, by = join_by(acc, contig))

df_hblD = read_tsv("hblD_blast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_") %>%
  left_join(lengths, by = join_by(acc, contig)) 

df_cytK1 = read_tsv("cytK1_blast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_") %>%
  left_join(lengths, by = join_by(acc, contig))

df_cytK2 = read_tsv("cytK2_blast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_") %>%
  left_join(lengths, by = join_by(acc, contig))

df_nheA = read_tsv("nheA_blast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_") %>%
  left_join(lengths, by = join_by(acc, contig))

df_nheB = read_tsv("nheB_blast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_") %>%
  left_join(lengths, by = join_by(acc, contig))

df_nheC = read_tsv("nheC_blast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_") %>%
  left_join(lengths, by = join_by(acc, contig))

df = data.frame(read_excel("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Biomarkers paper\\Mastersheet_no_clones.xlsx")) %>%
  select(acc = Isolate, N50, cytotoxicity = Average.Cell.Viability....0.7.is.cytotoxic., panC_group = Adjusted.panC.Group..predicted.species.)

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



#######

#i feel weird that there are so many more hblD than hblC, but # of hblD and hblA is correct... why are we missing hblC?
#are hblC and D in close proximity?



join_DC = left_join(df_hblD, df_hblC, by = join_by(acc, contig)) %>%
  left_join(df, by = join_by(acc == acc)) #295 obs

antijoin_DC = anti_join(df_hblD, df_hblC, by = join_by(acc, contig)) %>%
  left_join(df, by = join_by(acc == acc)) %>% #73 obs 
  left_join(lengths, by = join_by(acc, contig))

df_total_contigs = lengths %>%
  group_by(acc) %>%
  summarize(total_contigs = n())

sumD = df_hblD %>%
  group_by(acc) %>%
  summarize(total_hblD_hits = n()) %>%
  left_join(select(df, acc, N50, panC_group), by = join_by(acc == acc)) %>%
  left_join(df_total_contigs, by = join_by(acc == acc))




df_all_genes = df %>%
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
  replace(is.na(.), 0) 

df_all_genes_long  = df_all_genes %>%
  pivot_longer(nheA:cytK2, names_to = "gene", values_to = "count")


# Check to see if N50 correlates with the number of hblD and hblA hits I get.
fit = aov(N50 ~ factor(hblD), data = df_all_genes)
summary(fit)
TukeyHSD(fit)

fit = aov(N50 ~ factor(hblA), data = df_all_genes)
summary(fit)
TukeyHSD(fit)

fit = aov(N50 ~ factor(hblB), data = df_all_genes)
summary(fit)
TukeyHSD(fit)

ggplot(df_all_genes_long, aes(x = factor(count), y = N50))+
  geom_violin(width = 0.1) +
  scale_y_log10(n.breaks = 25, labels = scales::label_number(), limits = c(10000, 8000000)) +
  facet_wrap(~gene)

### Does panC group correlate w hbl gene numbers? is this the right stats test for this... i want to know more than presence/absence...

fit = aov(hblD ~ factor(panC_group), data = df_all_genes)
summary(fit)
tukey = TukeyHSD(fit)

ggplot(df_all_genes_long, aes(x = panC_group, fill = factor(count))) +
  geom_bar(position = "fill") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_fill_manual(values = c("0" = "gray50", "1" = "#b3cf99" , "2" = "#87ab69", "3" = "#4b6043")) +
  facet_wrap(~gene)


### Are the hbl variants similar in % id? oh absolutely lol
df_hblD %>% left_join(df) %>%
 ggplot(aes(x = panC_group, y = perc_id, color = panC_group)) +
  geom_jitter(width = 0.2) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(title = "hblD")


### Is cytotoxicity significantly different for 

fit = aov(cytotoxicity ~ factor(hblD), data = df_all_genes)
summary(fit)
tukey = TukeyHSD(fit)
print(tukey)



ggplot(df_all_genes_long, aes(x = factor(count), y = cytotoxicity)) +
  geom_violin()+
  geom_jitter(width = 0.2, aes(color = panC_group, alpha = 0.4)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = "gene copy number") +
  facet_wrap(~gene)

ggplot(df_all_genes_long, aes(x = factor(count), y = cytotoxicity, fill = panC_group)) +
  geom_boxplot() +
  labs(x = "gene copy number") +
  theme_classic() +
  facet_wrap(~gene) 

ggsave("enterotoxin_copy_boxplot_030626.png", width = 16, height = 8, units = "in")

### Separate variants based on percent id

df_hblD_top_col = df_hblD %>%
  group_by(acc) %>%
  mutate(top = (score == max(score))) %>% ###may need to change back to perc_id
  mutate(query = "hblD") %>%
  mutate(start = ifelse(s_start < s_end, s_start,s_end)) %>%
  mutate(stop = ifelse(s_start > s_end, s_start,s_end))

df_hblA_top_col = df_hblA_final %>%
  group_by(acc) %>%
  mutate(top = (perc_id == max(perc_id))) %>%
  mutate(query = "hblA") 

df_hblB_top_col = df_hblB_final %>%
  group_by(acc) %>%
  mutate(top = (perc_id == max(perc_id))) %>%
  mutate(query = sub("_.*", "", query))

rbind(df_hblD_top_col, df_hblA_top_col, df_hblB_top_col) %>%
  mutate(query = sub("_.*", "", query)) %>%
  drop_na(perc_id) %>%
  ggplot(aes(x = top, y= contig_length, color = query)) +
  geom_violin() +
  geom_jitter(width = 0.2) +
  scale_y_log10(n.breaks = 25, labels = scales::label_number(), limits = c(500, 8000000)) +
  labs(title = "Are duplicate genes an artifact of short contigs/poor assembly?", x = "Variant with highest % ID?", y = "Length of contig containing hit") +
  facet_wrap(~query)
  

### Need to cluster the variants somehow. pairwise alignments? or maybe the way i have it is good. pairwise align all of the variants per strain? maybe use score instead....

hblD_seqs = readDNAStringSet("fastas/hblD_all.fasta")%>%
  as.data.frame() %>%
  rownames_to_column() %>%
  separate(rowname, c("acc", "contig"), "_") %>%
  separate(contig, c("contig", "range"), ":") %>%
  separate(range, c("start", "stop"), "-") %>%
  mutate(start = as.numeric(start), stop = as.numeric(stop)) %>%
  rename_with(~"seq", .cols = last_col()) %>%
  mutate(query = "hblD") %>%
  left_join(df_hblD_top_col, by = join_by(query, acc, contig, start, stop)) %>%
  left_join(df_all_genes, by = join_by(acc))


ggplot(hblD_seqs, aes(x = perc_id, fill = top)) +
  geom_histogram() # insane looking graph but... lol could probably bin based on variant using a %id cutoff of 87. or score??? bin at 1300 ish




### Save df for other analyses/plots:

#write.csv(df_all_genes, file = "df_all_genes_030526.csv", row.names = FALSE)

#write_tsv(df_hblA_final, file = "df_hblA_filt_030526.tsv")
#write_tsv(df_hblB_final, file = "df_hblB_filt_030526.tsv")
