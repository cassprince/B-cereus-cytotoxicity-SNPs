library(tidyverse)

# Get list of files in directory
files = list.files(path="C:/Users/cassp/Documents/GitHub/B-cereus-cytotoxicity-SNPs/blast_results", pattern="*.csv", full.names=FALSE, recursive=FALSE)

# Loop through each file in directory 
lapply(files, function(x) {
  print(x)
  
  # Read in BLAST output. Filter hits with > 85% percent identity to the queries. 
  columns = c("query", "acc", "perc_id", "ali_len", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "eval", "score")
  df = read_tsv(x, skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_") %>%
  drop_na() %>%
  filter(perc_id > 85) %>%
  mutate(full_name = paste0(acc, "_", contig))%>%
  mutate(start = ifelse(s_start < s_end, s_start,s_end)) %>%
  mutate(stop = ifelse(s_start > s_end, s_start,s_end)) 
  
  # Group the BLAST data by isolate name + contig (aka "full_name"). Arrange each group by the "start" of each hit (ascending). 
  # The "start" here is not necessarily the start codon of the gene. Just the lower number between s_start and s_end in the BLAST data.
  df_full = df %>%
  group_by(full_name) %>%
  arrange(start, .by_group = TRUE) %>%
  group_split() # Split each group (isolate+contig) into individual test dfs that can be looped through.
  
  # Prepare empty dfs to append data onto.
  df_top_hits = data.frame()
  df_no_overlap = data.frame()
  
  # For each individual isolate+contig (test df)...
  for (df_test in df_full){
  df_test = data.frame(df_test)
  n = 2
  scores = c()
  num_no_overlap = 0
  
  # While there is still data in the test df...
  while (n < nrow(df_test)+1) {
    # If there is an overlap in genomic location between a hit and the hit directly below (which will be the next closest hit based on genomic location because of our sorting step above)...
    if (df_test$stop[n-1] - df_test$start[n] >= 0){
      scores = append(scores, df_test$perc_id[n-1]) # Record the percent identity in a list.

    } else { # If there isn't a genomic overlap...
      scores = append(scores, df_test$perc_id[n-1]) # Record the percent identity in a list.
      df_no_overlap = rbind(df_no_overlap, df_test[which.max(scores),]) # Save the hit with the highest percent identity.
      num_no_overlap = num_no_overlap + length(scores)
      scores = c() # Reset the score list.
    }
    n = n+1
    
  }
  scores = append(scores, df_test$perc_id[n-1]) # Record the percent identity in a list.
  df_top_hits = rbind(df_top_hits, df_test[which.max(scores)+num_no_overlap,]) #Save the hit with the highest percent identity.
  }
  
  df_final = rbind(df_no_overlap, df_top_hits) # Concatenate all hits.
  name_short = gsub(".csv", "", x)
  final_name = paste0(name_short, "_filt.csv")
  
  write_csv(df_final, paste0("C:/Users/cassp/Documents/GitHub/B-cereus-cytotoxicity-SNPs/blast_results/filtered/", final_name)) # Write filtered data as .csv.

})
