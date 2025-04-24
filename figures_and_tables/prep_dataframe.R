library(tidyverse)
setwd("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Biomarkers paper\\SNP hits\\Covars 8_14_23")


hblA = data.frame(read_csv("hblA_logreg_8_14_23.csv")) %>%
  select("Gene", "Position", "Nucleotide", "True.Negatives","False.Positives", "False.Negatives", "True.Positives", "Avg.Cytotoxicity.With", "Avg.Cytotoxicity.Without", "Stdev.Cytotoxicity.With", "Stdev.Cytotoxicity.Without", "Change.the.aa.residue.")
hblB = data.frame(read_csv("hblB_logreg_8_14_23.csv")) %>%
  select("Gene", "Position", "Nucleotide", "True.Negatives","False.Positives", "False.Negatives", "True.Positives", "Avg.Cytotoxicity.With", "Avg.Cytotoxicity.Without", "Stdev.Cytotoxicity.With", "Stdev.Cytotoxicity.Without", "Change.the.aa.residue.")
hblC = data.frame(read_csv("hblC_logreg_8_14_23.csv")) %>%
  select("Gene", "Position", "Nucleotide", "True.Negatives","False.Positives", "False.Negatives", "True.Positives", "Avg.Cytotoxicity.With", "Avg.Cytotoxicity.Without", "Stdev.Cytotoxicity.With", "Stdev.Cytotoxicity.Without", "Change.the.aa.residue.")
hblD = data.frame(read_csv("hblD_logreg_8_14_23.csv")) %>%
  select("Gene", "Position", "Nucleotide", "True.Negatives","False.Positives", "False.Negatives", "True.Positives", "Avg.Cytotoxicity.With", "Avg.Cytotoxicity.Without", "Stdev.Cytotoxicity.With", "Stdev.Cytotoxicity.Without", "Change.the.aa.residue.")

nheA = data.frame(read_csv("nheA_logreg_8_14_23.csv"))%>%
  select("Gene", "Position", "Nucleotide", "True.Negatives","False.Positives", "False.Negatives", "True.Positives", "Avg.Cytotoxicity.With", "Avg.Cytotoxicity.Without", "Stdev.Cytotoxicity.With", "Stdev.Cytotoxicity.Without", "Change.the.aa.residue.")
nheB = data.frame(read_csv("nheB_logreg_8_14_23.csv")) %>%
  select("Gene", "Position", "Nucleotide", "True.Negatives","False.Positives", "False.Negatives", "True.Positives", "Avg.Cytotoxicity.With", "Avg.Cytotoxicity.Without", "Stdev.Cytotoxicity.With", "Stdev.Cytotoxicity.Without", "Change.the.aa.residue.")
nheC = data.frame(read_csv("nheC_logreg_8_14_23.csv")) %>%
  select("Gene", "Position", "Nucleotide", "True.Negatives","False.Positives", "False.Negatives", "True.Positives", "Avg.Cytotoxicity.With", "Avg.Cytotoxicity.Without", "Stdev.Cytotoxicity.With", "Stdev.Cytotoxicity.Without", "Change.the.aa.residue.")

nhe_up = data.frame(read_csv("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Biomarkers paper\\upstream_4_21_25\\nhe_promoter_filt_logreg_4_21_25.csv")) %>%
  select("Gene", "Position", "Nucleotide", "True.Negatives","False.Positives", "False.Negatives", "True.Positives", "Avg.Cytotoxicity.With", "Avg.Cytotoxicity.Without", "Stdev.Cytotoxicity.With", "Stdev.Cytotoxicity.Without") %>%
  mutate("Change.the.aa.residue." = "FALSE")

hbl_up = data.frame(read_csv("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Biomarkers paper\\upstream_4_21_25\\hbl_promoter_filt_logreg_4_21_25.csv")) %>%
  select("Gene", "Position", "Nucleotide", "True.Negatives","False.Positives", "False.Negatives", "True.Positives", "Avg.Cytotoxicity.With", "Avg.Cytotoxicity.Without", "Stdev.Cytotoxicity.With", "Stdev.Cytotoxicity.Without")%>%
  mutate("Change.the.aa.residue." = "FALSE")

all = rbind(hblA, hblB, hblC, hblD, nheA, nheB, nheC, nhe_up, hbl_up) %>%
  mutate(Nucleotide = sub("(.)", "\\U\\1", Nucleotide, perl=TRUE)) %>%
  mutate(Gene = sub("_filt", "", Gene)) %>%
  mutate(label_nt = paste0(Gene, "_", Position, Nucleotide))
  

test = all %>%
  mutate(acc = (True.Positives + True.Negatives)/(True.Positives + False.Positives + False.Negatives + True.Negatives),
         prec = True.Positives/(True.Positives + False.Positives),
         sens = True.Positives/(True.Positives+ False.Negatives),
         spec = True.Negatives/(True.Negatives+ False.Positives), 
         score = acc+prec+sens+spec)

filt = test %>%
  filter(spec > 0.6 & sens > 0.6)
  
filt %>%
  summarize(total = n(),
            nonsynonymous = sum(as.logical(Change.the.aa.residue.)), 
            synonymous = n() - nonsynonymous)

write_csv(filt, "C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Biomarkers paper\\SNP hits\\4_21_25\\SNP_mastersheet_4_22_25.csv")



##### SNP_hit_mastersheet

snp_list = paste0(filt$Gene, "_", "SNP", "_", filt$Position)

hblA_snp = data.frame(read_csv("hblA_snps_8_14_23.csv")) %>% 
  select(-Adjusted.panC.Group..predicted.species., -Average.Cell.Viability....0.7.is.cytotoxic., -...1)
hblB_snp = data.frame(read_csv("hblB_snps_8_14_23.csv")) %>% 
  select(-Adjusted.panC.Group..predicted.species., -Average.Cell.Viability....0.7.is.cytotoxic., -...1)
hblC_snp = data.frame(read_csv("hblC_snps_8_14_23.csv")) %>% 
  select(-Adjusted.panC.Group..predicted.species., -Average.Cell.Viability....0.7.is.cytotoxic., -...1)
hblD_snp = data.frame(read_csv("hblD_snps_8_14_23.csv")) %>% 
  select(-Adjusted.panC.Group..predicted.species., -Average.Cell.Viability....0.7.is.cytotoxic., -...1)
nheA_snp = data.frame(read_csv("nheA_snps_8_14_23.csv")) %>% 
  select(-Adjusted.panC.Group..predicted.species., -Average.Cell.Viability....0.7.is.cytotoxic., -...1)
nheB_snp = data.frame(read_csv("nheB_snps_8_14_23.csv")) %>% 
  select(-Adjusted.panC.Group..predicted.species., -Average.Cell.Viability....0.7.is.cytotoxic., -...1)
nhe_up_snp = data.frame(read_csv("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Biomarkers paper\\upstream_4_21_25\\nhe_promoter_filt_snps_4_21_25.csv")) 
colnames(nhe_up_snp) = sub("_filt", "", colnames(nhe_up_snp))

df_master = data.frame(read_excel("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Biomarkers paper\\Mastersheet_no_clones.xlsx")) %>% 
  select(Isolate, Adjusted.panC.Group = Adjusted.panC.Group..predicted.species., Average.Cell.Viability = Average.Cell.Viability....0.7.is.cytotoxic.) %>%
  mutate(Adjusted.panC.Group = gsub("[*].*$", "", Adjusted.panC.Group))

df_genes = data.frame(read_csv("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Biomarkers paper\\SNP_gene_hits.csv")) %>%
  select(Isolate = "...1", nhe, hbl, cytk1, cytk2) 

var_names = read.csv("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Biomarkers paper\\SNP hits\\4_21_25\\SNP_mastersheet_labels_4_22_25.csv") %>% 
  mutate(label_snp = paste0(Gene, "_", "SNP", "_", Position)) %>%
  select(label_final, label_snp) %>%
  mutate(label_snp = sub("_filt", "", label_snp)) %>% 
  deframe()

df_snp_final = left_join(df_master, hblA_snp, by = join_by("Isolate")) %>%
  left_join(hblB_snp, by = join_by("Isolate")) %>%
  left_join(hblC_snp, by = join_by("Isolate")) %>%
  left_join(hblD_snp, by = join_by("Isolate")) %>%
  left_join(nheA_snp, by = join_by("Isolate")) %>%
  left_join(nheB_snp, by = join_by("Isolate")) %>%
  left_join(nhe_up_snp, by = join_by("Isolate")) %>% 
  left_join(df_genes, by = join_by("Isolate")) %>% 
  select(Isolate, Adjusted.panC.Group, Average.Cell.Viability, nhe, hbl, cytk1, cytk2, snp_list) %>%
  mutate(Adjusted.panC.Group, Adjusted.panC.Group = recode(Adjusted.panC.Group, "Group_clarus" = "0", "Group_I(pseudomycoides)" = "1", "Group_II(mosaicus/luti)" = "2", "Group_III(mosaicus)" = "3", "Group_IV(cereus_sensu_stricto)" = "4", "Group_V(toyonensis)" = "5", "Group_VI(mycoides/paramycoides)" = "6", "Group_VII(cytotoxicus)" = "7", "Group_VIII(mycoides)" = "8")) %>%
  rename(!!!var_names)


write_csv(df_snp_final, "C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Biomarkers paper\\SNP hits\\4_21_25\\SNP_hits_4_22_25.csv")
