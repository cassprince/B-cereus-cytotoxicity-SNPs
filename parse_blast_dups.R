library(tidyverse)
library(ivs)

setwd("C:/Users/cassp/Documents/GitHub/B-cereus-cytotoxicity-SNPs/blast_results")

columns = c("query", "acc", "perc_id", "ali_len", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "eval", "score")
df = read_tsv("nheA_blast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_") %>%
  drop_na() %>%
  filter(perc_id > 85) %>%
  mutate(full_name = paste0(acc, "_", contig))%>%
  mutate(start = ifelse(s_start < s_end, s_start,s_end)) %>%
  mutate(stop = ifelse(s_start > s_end, s_start,s_end)) 

df_full = df %>%
  group_by(full_name) %>%
  arrange(start, .by_group = TRUE) %>%
  group_split()

df_top_hits = data.frame()
df_no_overlap = data.frame()

for (df_test in df_full){
  df_test = data.frame(df_test)
  #print(df_test[1,2])
  n = 2
  scores = c()
  num_no_overlap = 0
  
  while (n < nrow(df_test)+1) {
    #print(n)
    if (df_test$stop[n-1] - df_test$start[n] > 0){
      #print("overlap")
      scores = append(scores, df_test$perc_id[n-1])
      #print(paste("current max index in scores:", which.max(scores)))
      #print(paste("current max index in dt_test:", which.max(scores)+num_no_overlap))
      
    } else {
      #print("no")
      num_no_overlap = num_no_overlap + 1
      scores = c()
      df_no_overlap = rbind(df_no_overlap, df_test[n-1,])
    }
    n = n+1
    
  }
  scores = append(scores, df_test$perc_id[n-1])
  #print(paste("final max index in scores:", which.max(scores)))
  #print(paste("final max index in dt_test:", which.max(scores)+num_no_overlap))
  df_top_hits = rbind(df_top_hits, df_test[which.max(scores)+num_no_overlap,])
}






