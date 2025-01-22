library(ggplot2)
library(ggprism)
library(ggbreak)
library(runner)
library(dplyr)
library(stringr)
library(viridis)
setwd("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Biomarkers paper\\Protein alignments 8_21_23")


############################################
#-----------AA IDENTITY PLOT------------#

df <- data.frame(read.csv("cytK2_identity_8_21_23.csv"))

running_mean <- runner(df[,3], k = 1, f = function(x) mean(x)) * 100

df <- cbind(df, Mean = running_mean)
df <- mutate(df, Score = running_mean * Coverage)

xticks <- c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450)
#xticks <- c(240, 243, 246, 250)
yticks <- c(0, 20, 40, 60, 80, 100)

p = ggplot(df, aes(x = Position)) +
  geom_area(mapping = aes(y = ifelse(Mean > 90, Mean, 0)),  fill = "gray") +  
  geom_line(aes(y = Mean), color = "black", linewidth = 0.75) +
  theme_prism() +
  scale_y_continuous(expand= c(0,0), limits = c(0, 105), breaks = yticks) +
  scale_x_continuous(expand= c(0,0), breaks = xticks) +
  xlab("Position (aa)") +
  ylab("Average percent identity (%)") +
  theme(text = element_text(size = 20))

p

ggsave("cytK2_aa_identity_bigtext.png", p, width = 12, height = 4.5)

