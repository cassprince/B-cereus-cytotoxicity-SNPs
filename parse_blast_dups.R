library(tidyverse)

setwd("C:/Users/cassp/Documents/GitHub/B-cereus-cytotoxicity-SNPs/blast_results")

files = list.files(path="C:/Users/cassp/Documents/GitHub/B-cereus-cytotoxicity-SNPs/blast_results", pattern="*.csv", full.names=FALSE, recursive=FALSE)

lapply(files, function(x) {
  print(x)
  
  columns = c("query", "acc", "perc_id", "ali_len", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "eval", "score")
  df = read_tsv(x, skip = 5, col_names = columns) %>%
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
    if (df_test$stop[n-1] - df_test$start[n] >= 0){
      #print("overlap")
      scores = append(scores, df_test$perc_id[n-1])
      #print(paste("current max index in scores:", which.max(scores)))
      #print(paste("current max index in dt_test:", which.max(scores)+num_no_overlap))
      
    } else {
      scores = append(scores, df_test$perc_id[n-1])
      #print("else")
      #print(df_test[1,2])
      #print(scores)
      #print(paste("final max index in scores:", which.max(scores)))
      
      df_no_overlap = rbind(df_no_overlap, df_test[which.max(scores),])
      num_no_overlap = num_no_overlap + length(scores)
      scores = c()
    }
    n = n+1
    
  }
  scores = append(scores, df_test$perc_id[n-1])
  #print(df_test[1,2])
  #print(scores)
  #print(paste("number of no overlap:", num_no_overlap))
  #print(paste("final max index in scores:", which.max(scores)))
  #print(paste("final max index in dt_test:", which.max(scores)+num_no_overlap))
  df_top_hits = rbind(df_top_hits, df_test[which.max(scores)+num_no_overlap,])
  }
  
  df_final = rbind(df_no_overlap, df_top_hits)
  
  name_short = gsub(".csv", "", x)
  final_name = paste0(name_short, "_filt.csv")
  
  #write_csv(df_final, paste0("C:/Users/cassp/Documents/GitHub/B-cereus-cytotoxicity-SNPs/blast_results/filtered/", final_name))

})

dist = df_final %>%
  distinct(acc)

two = df_final %>%
  group_by(acc) %>%
  filter(n() > 1)
