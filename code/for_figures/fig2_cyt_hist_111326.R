library(ggprism)
library(diptest)
library(ggpubr)
library(cowplot)
library(colorspace)
library(tidyverse)

setwd("C:/Users/princecr/Documents/GitHub/B-cereus-cytotoxicity-SNPs/data")
data= read_csv("Mastersheet_082026.csv") %>%
  select(Isolate, cytotoxicity = "Cytotoxicity ( >0.7 is cytotoxic)", panC_group = "Adjusted_panC_Group(predicted_species)") %>%
  mutate(panC_group = gsub("\\s*\\([^\\)]+\\)","",panC_group)) %>%
  mutate(panC_group = gsub("\\*", "", panC_group)) %>%
  mutate(panC_group = gsub("_", " ", panC_group)) %>%
  drop_na()

ticks=c(-0.5, 0, 0.5, 1, 1.5)

# Test for unimodality/normal distribution (p<0.05 is not unimodal)
dip.test(data$cytotoxicity)

ggplot(data, aes(x=cytotoxicity))+ 
  geom_histogram(binwidth=0.1, color="black", fill="gray", aes(y=..density..))+
  #geom_vline(aes(xintercept=0.7), color="black", linetype="dashed", linewidth=1)+
  geom_density(color="#56B4E9", linewidth = 0.75)+
  theme_prism(base_size = 17)+
  scale_y_continuous(expand= c(0,0))+
  xlab("Normalized Cytotoxicity")+
  ylab("Kernel Density Estimate")+
  scale_x_continuous(limits = c(-0.5, 1.5), breaks= ticks) +
  theme(axis.title = element_text(size = 18))

ggsave("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Figures\\Biomarkers\\Cytotoxicity histogram 2_15_24.png", plot = plt, device = "png", width = 8, height = 4.5, dpi = 300)




ggplot(data, aes(x=cytotoxicity))+ 
  geom_histogram(binwidth=0.1, color="black", aes(fill=panC_group))+
  scale_y_continuous(expand= c(0,0))+
  labs(x = "Normalized Cytotoxicity", y = "Frequency")+
  scale_x_continuous(limits = c(-0.5, 1.5), breaks= ticks) +
  theme_classic()+
  theme(text = element_text(color = "black", size = 14), 
        axis.text = element_text(color = "black", size = 14), 
        axis.title.y = element_text(color = "black", size = 15), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "italic")) +
  scale_fill_manual(name = "", values = c("gray", "#FC8383", "#FCC983", "#FCF083", "#B3FC83", "#B5DCFF", "#A0A5F7", "#C1A0F7", "#F7A0E7"))

ggsave("../figures/fig2_histogram_111326.png",  width = 7, height = 3.5, units = "in")
