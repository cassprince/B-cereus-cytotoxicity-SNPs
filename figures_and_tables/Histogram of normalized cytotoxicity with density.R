library(ggplot2)
library(readxl)
library(ggprism)
library(diptest)
library(ggpubr)
library(cowplot)
library(colorspace)

setwd("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Biomarkers paper")
data<-data.frame(read_excel("Mastersheet_no_clones.xlsx"))

cytotoxicity<-data[,"Average.Cell.Viability....0.7.is.cytotoxic."]
ticks<-c(-0.5, 0, 0.5, 1, 1.5)
group<-data[,"Adjusted.panC.Group..predicted.species."]
group<-gsub("\\s*\\([^\\)]+\\)", "", group)
data$Adjusted.panC.Group..predicted.species.<-gsub("_", " ", group)



dip.test(cytotoxicity)

plt <- ggplot(data, aes(x=cytotoxicity))+ 
  geom_histogram(binwidth=0.1, color="black", fill="gray", aes(y=..density..))+
  #geom_vline(aes(xintercept=0.7), color="black", linetype="dashed", linewidth=1)+
  geom_density(color="#56B4E9", linewidth = 0.75)+
  theme_prism(base_size = 17)+
  scale_y_continuous(expand= c(0,0))+
  xlab("Normalized Cytotoxicity")+
  ylab("Kernel Density Estimate")+
  scale_x_continuous(limits = c(-0.5, 1.5), breaks= ticks) +
  theme(axis.title = element_text(size = 18))

plt

ggsave("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Figures\\Biomarkers\\Cytotoxicity histogram 2_15_24.png", plot = plt, device = "png", width = 8, height = 4.5, dpi = 300)




stacked = ggplot(data, aes(x=Average.Cell.Viability....0.7.is.cytotoxic.))+ 
  geom_histogram(binwidth=0.1, color="black", aes(fill=Adjusted.panC.Group..predicted.species.))+
  theme_prism(base_size = 17)+
  scale_y_continuous(expand= c(0,0))+
  labs(x = "Normalized Cytotoxicity", y = "Frequency")+
  scale_x_continuous(limits = c(-0.5, 1.5), breaks= ticks) +
  theme(axis.title = element_text(size = 18)) +
  scale_fill_discrete_qualitative(palette = "Set2")

ggsave("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Figures\\Biomarkers\\Cytotoxicity histogram 9_2_24.png", plot = stacked, device = "png", width = 9, height = 4.5, dpi = 300)
