library(tidyverse)

setwd("C:/Users/cassp/OneDrive/Documents/GitHub/B-cereus-cytotoxicity-SNPs/blast_results/filtered")

df_hblB = read_csv("hblB_blast_filt.csv") %>%
  mutate(gene = "hblB")

df_hblA = read_csv("hblA_blast_filt.csv") %>%
  mutate(gene = "hblA")

df_arranged = rbind(df_hblA, df_hblB) %>%
  group_by(full_name) %>%
  arrange(start, .by_group = TRUE) 

df_full = df_arranged %>%
  group_split()

df_top_hits = data.frame()
df_no_overlap = data.frame()


for (df_test in df_full){
  df_test = data.frame(df_test)
  n = 1
  scores = data.frame()
  
  print(paste(df_test[1,'full_name'], "nrows of df_test:", nrow(df_test)))
  while (n < nrow(df_test)) {
    if (nrow(scores) == 0){ #if first instance
      scores = df_test[n,]
    }
    
    if (df_test$stop[n] - df_test$start[n+1] >= 0){ #if overlap
      scores = rbind(scores, df_test[n+1,]) 
      
      if (n+1 == nrow(df_test)){ #if last possible comparison
        df_top_hits = rbind(df_top_hits, scores[which.max(scores$perc_id),]) #determine max in list of scores
      }
      
    } 
    else { #no overlap between n and n+1
      if (nrow(scores)>1){ #if there were overlapping genes directly preceding the non overlap
          df_top_hits = rbind(df_top_hits, scores[which.max(scores$perc_id),])
          scores = data.frame()
        } 
      else{ #no overlapping genes were directly preceding the non overlap
          df_no_overlap = rbind(df_no_overlap, scores) #save first gene tested in non overlap
          scores = data.frame()
        }
      if (n+1 == nrow(df_test)){ #if last possible comparison
        df_no_overlap = rbind(df_no_overlap, df_test[n+1,]) #save last gene as non-overlap
      }
    }
    
    n = n+1  
  }
  
}
  
  
  
df_AB = rbind(df_no_overlap, df_top_hits)

df_hblA_final = df_AB %>%
  filter(gene == "hblA")

df_hblB_final = df_AB %>%
  filter(gene == "hblB")

df_hblB_final %>%
  group_by(acc) %>%
  filter(n() >1)
  
df_hblA_final %>%
  group_by(acc) %>%
  filter(n() >1)
