library(tidyverse)

setwd("C:/Users/cassp/OneDrive/Documents/GitHub/B-cereus-cytotoxicity-SNPs/blast_results")

# Get list of files in directory
files = list.files(path="C:/Users/cassp/OneDrive/Documents/GitHub/B-cereus-cytotoxicity-SNPs/blast_results", pattern="*.csv", full.names=FALSE, recursive=FALSE)

x = "nheA_blast.csv"

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

  # For each individual isolate+contig (test df):
for (df_test in df_full){ 
    df_test = data.frame(df_test)
    print(df_test[1,2])
    n = 1
    scores = data.frame() # Create a df of rows to compare between the % id's. This will only be added to in the case of overlapping gene pairs.
    
    # If df_test only has one hit (no comparisons can be made):
    if (nrow(df_test) == 1){
      df_top_hits = rbind(df_top_hits, df_test)
    }
    # While there is still data in the test df:
    while (n < nrow(df_test)) {
      
      #If this is the first instance of df_test OR after a string of overlapping genes:
      if (nrow(scores) == 0){ 
        scores = df_test[n,] # Record data for the row.
      }
      
      # If there is an overlap in genomic location between a hit and the hit directly below (which will be the next closest hit based on genomic location because of our sorting step above):
      if (df_test$stop[n] - df_test$start[n+1] >= 0){ 
        scores = rbind(scores, df_test[n+1,]) 
        
        # If n+1 is the last row in df_test (last possible comparison):
        if (n+1 == nrow(df_test)){ 
          df_top_hits = rbind(df_top_hits, scores[which.max(scores$perc_id),]) # Record the row with max % id in df of scores.
        }
        
      } else { # If no overlap between n and n+1:
        # If there were overlapping genes directly preceding the non-overlapping gene pair:
        if (nrow(scores)>1){ 
            df_top_hits = rbind(df_top_hits, scores[which.max(scores$perc_id),]) # Record the row with max % id in df of scores.
            scores = data.frame() # Reset the df of scores.
          } else{ # If no overlapping gene pairs were directly preceding the non-overlapping gene pair:
          df_no_overlap = rbind(df_no_overlap, df_test[n,]) # Save the first gene in the tested gene pair
            scores = data.frame()
        }
        # If n+1 is the last row in df_test (last possible comparison):
        if (n+1 == nrow(df_test)){ 
          df_no_overlap = rbind(df_no_overlap, df_test[n+1,]) # Save the second gene in the tested gene pair too 
        }
      }
      
      n = n+1  
    }
    
}
  
df_final = rbind(df_no_overlap, df_top_hits) # Concatenate all hits.
