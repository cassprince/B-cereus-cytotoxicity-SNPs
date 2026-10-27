library(tidyverse)
library(readxl)
library(ivs)
library(ggridges)
library(Biostrings)
library(msa)

setwd("C:/Users/cassp/OneDrive/Documents/GitHub/B-cereus-cytotoxicity-SNPs/blast_results/filtered_megablast_qc50_090126")

lengths = read_tsv("C:/Users/cassp/OneDrive/Documents/GitHub/B-cereus-cytotoxicity-SNPs/contig_lengths_082826.txt", col_names = c("full_name", "contig_length")) 


df_hblB = read_csv("hblB_parsed_090126.csv") %>%
  mutate(gene = "hblB")

df_hblA = read_csv("hblA_parsed_090126.csv") %>%
  mutate(gene = "hblA")

df_hblC = read_csv("hblC_blast_090126_filt.csv")

df_hblD = read_csv("hblD_blast_090126_filt.csv")

df_cytK1 = read_csv("cytK1_blast_090126_filt.csv")

df_cytK2 = read_csv("cytK2_blast_090126_filt.csv")

df_nheA = read_csv("nheA_blast_090126_filt.csv")

df_nheB = read_csv("nheB_blast_090126_filt.csv")

df_nheC = read_csv("nheC_blast_090126_filt.csv")

df = data.frame(read_excel("C:/Users/cassp/OneDrive - Cornell University/Biomarkers paper/Mastersheet_082026.xlsx")) %>%
  select(acc = Isolate, cytotoxicity = Cytotoxicity....0.7.is.cytotoxic., panC_group = Adjusted_panC_Group.predicted_species.) %>%
  mutate(panC_group = gsub("\\s*\\([^\\)]+\\)","",panC_group)) %>%
  mutate(panC_group = gsub("\\*", "", panC_group))

df_quast = read_tsv("C:/Users/cassp/OneDrive - Cornell University/Biomarkers paper/quast_results_09_01_26/transposed_report.tsv")

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
  left_join(df_hblA %>% group_by(acc) %>% summarize(hblA = n()), 
            by = join_by(acc)) %>%
  left_join(df_hblB %>% group_by(acc) %>% summarize(hblB = n()), 
            by = join_by(acc)) %>%
  left_join(df_cytK1 %>% group_by(acc) %>% summarize(cytK1 = n()), 
            by = join_by(acc)) %>%
  left_join(df_cytK2 %>% group_by(acc) %>% summarize(cytK2 = n()), 
            by = join_by(acc)) %>%
  replace(is.na(.), 0) %>%
  left_join(df_quast %>% select(Assembly, num_contigs = "# contigs", N50, N90), by = join_by(acc == "Assembly"))

df_all_genes_long  = df_all_genes %>%
  pivot_longer(nheA:cytK2, names_to = "gene", values_to = "count")


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


df_hblC %>% left_join(df) %>%
  ggplot(aes(x = panC_group, y = perc_id, color = panC_group)) +
  geom_jitter(width = 0.2) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(title = "hblC")


### Is cytotoxicity significantly different for 

fit = aov(cytotoxicity ~ factor(hblD), data = df_all_genes)
summary(fit)
tukey = TukeyHSD(fit)
print(tukey)

df_hblD %>% left_join(df) %>%
  ggplot(aes(x = perc_id, y = cytotoxicity, color = panC_group)) +
  geom_point() +
  stat_smooth(method = lm)

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

df_hblA_top_col = df_hblA %>%
  group_by(acc) %>%
  mutate(top = (perc_id == max(perc_id))) %>%
  mutate(query = "hblA") 

df_hblB_top_col = df_hblB %>%
  group_by(acc) %>%
  mutate(top = (perc_id == max(perc_id))) %>%
  mutate(query = sub("_.*", "", query))

df_hblC_top_col = df_hblC %>%
  group_by(acc) %>%
  mutate(top = (perc_id == max(perc_id))) %>%
  mutate(query = sub("_.*", "", query))

  
ggplot(df_hblD_top_col %>% left_join(df), aes(x = perc_id, fill = top)) +
  geom_histogram() +
  facet_wrap(~panC_group)


ggplot(df_hblC_top_col %>% left_join(df), aes(x = perc_id, fill = top)) +
  geom_histogram() +
  facet_wrap(~panC_group)

ggplot(df_hblB %>% left_join(df), aes(x = perc_id)) +
  geom_histogram() +
  facet_wrap(~panC_group)


### Save df for other analyses/plots:

write.csv(df_all_genes, file = "C:/Users/cassp/OneDrive/Documents/GitHub/B-cereus-cytotoxicity-SNPs/df_all_genes_090126.csv", row.names = FALSE)

#write_tsv(df_hblA_final, file = "df_hblA_filt_030526.tsv")
#write_tsv(df_hblB_final, file = "df_hblB_filt_030526.tsv")




