library(tidyverse)
library(ggdist)
library(ggsignif)

setwd("C:/Users/princecr/Documents/GitHub/B-cereus-cytotoxicity-SNPs/data")
df_genes = read_csv("df_all_genes_090126.csv") %>%
  mutate(nheABC = ifelse(nheA>0 & nheB>0 & nheC >0, "1", "0")) %>%
  mutate(hblCDAB = ifelse(hblC>0 & hblD>0 & hblA>0 & hblB >0, "1", "0")) %>%
  mutate(across(c(nheABC, hblCDAB, cytK1, cytK2), ~ as.factor(.)))

df_plot = df_genes %>%
  pivot_longer(cols = c(nheABC, hblCDAB, cytK1, cytK2), names_to = "genes", values_to = "presence") %>%
  mutate(presence = gsub("1", "present", presence)) %>%
  mutate(presence = gsub("0", "absent", presence)) %>%
  mutate(genes = gsub("cytK1", "cytK-1", genes)) %>%
  mutate(genes = gsub("cytK2", "cytK-2", genes)) %>%
  mutate(genes = factor(genes, levels = c("nheABC", "hblCDAB", "cytK-1", "cytK-2"))) 
df_summ = df_plot %>%
  group_by(genes, presence) %>%
  summarize(n = n())


ggplot(df_plot, aes(x = presence, y = cytotoxicity)) +
  geom_hline(yintercept = 0.51, color = "gray", linetype = 1) +
  stat_halfeye(justification = -0.9, width = 0.3, .width = 0) +
  geom_boxplot(width = 0.4, aes(fill = presence)) +
  geom_signif(comparisons = list(c("present", "absent")), y_position = 1.35) +
  geom_text(data = df_summ, aes(x = presence, y = 1.65, label = paste0("n=", n)), size = 4) +
  theme_classic() +
  theme(text = element_text(color = "black", size = 14), 
        axis.text = element_text(color = "black", size = 14), 
        axis.title.y = element_text(color = "black", size = 15), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        strip.text = element_text(face = "italic")) +
  labs(x = "", y = "Normalized Cytotoxicity") +
  scale_fill_manual(values = c("#FCC983", "#B5DCFF")) +
  facet_wrap(~genes, nrow = 1)

ggsave("../fig3_raincloud_genes_111126.png", width = 7.4, height = 4.3, units = "in")



#scale_color_manual(values = c("gray", "#FC8383", "#FCC983", "#FCF083", "#B3FC83", "#B5DCFF", "#A0A5F7", "#C1A0F7", "#F7A0E7")) +
