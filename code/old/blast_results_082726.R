library(tidyverse)
library(readxl)
library(ivs)
library(ggridges)
library(Biostrings)
library(msa)

setwd("C:/Users/cassp/OneDrive/Documents/GitHub/B-cereus-cytotoxicity-SNPs/blast_results/filtered/attempt_082726")

#lengths = read_tsv("contigs_lengths.txt", col_names = c("acc", "contig_length")) %>%
  #separate(acc, c("acc", "contig"), "_")

df_hblB = read_csv("hblB_blast_filt.csv") %>%
  mutate(gene = "hblB")

df_hblA = read_csv("hblA_blast_filt.csv") %>%
  mutate(gene = "hblA")

df_hblC = read_csv("hblC_blast_filt.csv")

df_hblD = read_csv("hblD_blast_filt.csv")

df_cytK1 = read_csv("cytK1_blast_filt.csv")

df_cytK2 = read_csv("cytK2_blast_filt.csv")

df_nheA = read_csv("nheA_blast_filt.csv")

df_nheB = read_csv("nheB_blast_filt.csv")

df_nheC = read_csv("nheC_blast_filt.csv")

df = data.frame(read_excel("../../Mastersheet_082026.xlsx")) 

#%>%
  #select(acc = Isolate, cytotoxicity = Cytotoxicity....0.7.is.cytotoxic., panC_group = Adjusted_panC_Group.predicted_species.)



nheA_shared = df %>% 
  filter(grepl("nheA", df$diarrheal_toxin_Nhe.genes.)) %>% 
  full_join(df_nheA, by = join_by(Isolate == acc))
nheB_shared = df %>% 
  filter(grepl("nheB", df$diarrheal_toxin_Nhe.genes.)) %>% 
  full_join(df_nheB, by = join_by(Isolate == acc))
nheC_shared = df %>% 
  filter(grepl("nheC", df$diarrheal_toxin_Nhe.genes.)) %>% 
  full_join(df_nheC, by = join_by(Isolate == acc))

sum(grepl("cytK-1", df$diarrheal_toxin_CytK.top_hit.))
sum(grepl("cytK-2", df$diarrheal_toxin_CytK.top_hit.))

files = list.files(path="C:/Users/cassp/OneDrive - Cornell University/Biomarkers paper/genomes_08_20_26", full.names=FALSE, recursive=FALSE)

files_notdf = anti_join(data.frame(gsub(".fasta", "", files)), df, by = join_by(gsub...fasta.......files. == Isolate))

df_notfiles = anti_join(df, data.frame(gsub(".fasta", "", files)), by = join_by(Isolate == gsub...fasta.......files.))

# Distinguish whether a shared gene hit is more like hblA or hblB and label it as such.

df_arranged = rbind(df_hblA, df_hblB) %>%
  group_by(full_name) %>%
  arrange(start, .by_group = TRUE) 

df_full = df_arranged %>%
  group_split()

df_top_hits = data.frame()
df_no_overlap = data.frame()

for (df_test in df_full){
  df_test = data.frame(df_test)
  n = 2
  scores = c()
  num_no_overlap = 0
  print(paste(df_test[1,'full_name'], "nrows of df_test:", nrow(df_test)))

  # While n-1 isn't the last row (aka there is still data in the test df)...
  while (n-1 < nrow(df_test)) {
    # If there is an overlap in genomic location between a hit and the hit directly below (which will be the next closest hit based on genomic location because of our sorting step above)...
    
    if (df_test$stop[n-1] - df_test$start[n] >= 0){
      scores = append(scores, df_test$perc_id[n-1]) # Record the percent identity in a list.
      
    } else { # If there isn't a genomic overlap...
        if (length(scores) == 0){
          scores = append(scores, df_test$perc_id[n-1]) # Record the first percent identity in a list.
        } else {
        print("no overlap")
        
        df_no_overlap = rbind(df_no_overlap, df_test[n,]) # Save the non overlapping hit.
        num_no_overlap = num_no_overlap + 1
        scores = c() # Reset the score list.
        }
    }
    n = n+1
    
  }
  scores = append(scores, df_test$perc_id[n-1]) # Record the percent identity in a list.
  print(scores)
  print(paste("max index:", which.max(scores), "... number of no overlaps:", num_no_overlap))
  print(paste("max perc_id:", df_test[which.max(scores)+num_no_overlap, 'perc_id']))
  df_top_hits = rbind(df_top_hits, df_test[which.max(scores)+num_no_overlap,]) #Save the hit with the highest percent identity.
}




df_AB = rbind(df_no_overlap, df_top_hits) %>%
  distinct(perc_id, full_name, start, stop, gene, .keep_all = TRUE)

df_hblA_final = df_AB %>%
  filter(gene == "hblA")

df_hblB_final = df_AB %>%
  filter(gene == "hblB")



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
