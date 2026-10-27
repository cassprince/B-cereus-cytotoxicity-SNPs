library(tidyverse)
library(Biostrings)
library(seqinr)


########## shorten hblC

hblC_up_seqs = readDNAStringSet("upstream_seqs/hblC_up.fasta")%>%
  as.data.frame() %>%
  rownames_to_column() %>%
  #separate(rowname, c("contig", "range"), ":") %>%
  #separate(range, c("start", "stop"), "-") %>%
  mutate(rowname = gsub(" <unknown description>", "", rowname)) %>%
  #mutate(start = as.numeric(start), stop = as.numeric(stop)) %>%
  rename_with(~"seq", .cols = last_col()) 

hblC_up_short = hblC_up_seqs %>%
  mutate(seq_length = str_length(seq)) %>%
  mutate(subseq = substr(seq, str_length(seq)-600+1, str_length(seq))) %>%
  mutate(subseq_length = str_length(subseq)) %>%
  filter(subseq_length > 50)


write.fasta(as.list(hblC_up_short$subseq), names = as.list(hblC_up_short$rowname), file.out = "upstream_seqs/hblC_up_600bp.fasta", nbchar = 60)



