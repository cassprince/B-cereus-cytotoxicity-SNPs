library(tidyverse)

setwd("C:/Users/cassp/OneDrive/Documents/Kovac Lab/Biomarkers paper/genes_nano_031826")

columns = c("query", "acc", "perc_id", "ali_len", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "eval", "score")

lengths = read_tsv("contig_lengths_nano.txt", col_names = c("acc", "contig_length")) %>%
  separate(acc, c("acc", "contig"), "_")

df_hblD = read_tsv("hblD_blast.csv", skip = 5, col_names = columns) %>%
  separate(acc, c("acc", "contig"), "_")  %>%
  left_join(lengths, by = join_by(acc, contig)) 

sumD = df_hblD %>%
  filter(query == "hblD_NZ_JBJJSR010000001.1:c3065421-3064201") %>%
  group_by(acc) %>%
  summarize(total_hblD_hits = n()) 

df_hblD %>%
  filter(query == "hblD_NZ_JBJJSR010000001.1:c3065421-3064201") %>%
  group_by(acc) %>%
  mutate(hblD_copy_number = as.character(n())) %>%
  ggplot(aes(x = perc_id, y = contig_length, color = hblD_copy_number)) +
  geom_point(alpha = 0.5) +
  scale_y_log10(n.breaks = 25, labels = scales::label_number()) +
  scale_x_break(c(85, 95), scales = 0.7) +
  labs(x = "% identity to B. cereus s.s. hblD", y = "contig length (bp)")

ggsave("hblD_copy_contig_031826.png", height = 5, width = 8, units = "in")


write.csv(sumD, file = "df_hblD_031826.csv", row.names = FALSE)
