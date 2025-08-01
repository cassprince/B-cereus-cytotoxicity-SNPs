library(ggplot2)
library(ggprism)
library(ggbreak)
library(runner)
library(dplyr)
library(stringr)
library(viridis)
setwd("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Biomarkers paper\\upstream_4_21_25")


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

#-----------NT IDENTITY PLOT------------#
nheUp = data.frame(read.csv("nheUp_identity_4_25_25.csv"))
running_mean = runner(nheUp['Identity'], k = 1, f = function(x) mean(x)) * 100
nheUp = cbind(nheUp, Mean = running_mean, gene = "nheUp") 

nheA <- data.frame(read.csv("nheA_identity_4_25_25.csv"))
running_mean = runner(nheA['Identity'], k = 1, f = function(x) mean(x)) * 100
nheA = cbind(nheA, Mean = running_mean, gene = "nheA") 

nheB <- data.frame(read.csv("nheB_identity_4_25_25.csv"))
running_mean = runner(nheB['Identity'], k = 1, f = function(x) mean(x)) * 100
nheB = cbind(nheB, Mean = running_mean, gene = "nheB") 

nheC <- data.frame(read.csv("nheC_identity_4_25_25.csv"))
running_mean = runner(nheC['Identity'], k = 1, f = function(x) mean(x)) * 100
nheC = cbind(nheC, Mean = running_mean, gene = "nheC") 

nhe_operon = rbind(nheUp, nheA, nheB, nheC)

xticks <- seq(from = 0, to = 1500, by = 250)
yticks <- c(0, 20, 40, 60, 80, 100)

nhe_operon$gene = factor(nhe_operon$gene, levels = c("nheUp", "nheA", "nheB", "nheC"))

p = ggplot(nhe_operon, aes(x = Position)) +
  geom_area(mapping = aes(y = ifelse(Mean > 90, Mean, 0)),  fill = "gray") +
  geom_line(aes(y = Mean), color = "black", linewidth = 0.75) +
  theme_prism() +
  scale_y_continuous(expand= c(0,0), limits = c(0, 105), breaks = yticks) +
  scale_x_continuous(expand= c(0,0), breaks = xticks) +
  labs(x = "", y = "") +
  #xlab("Position (nt)") +
  #ylab("Average percent identity (%)") +
  theme(text = element_text(size = 36)) +
  facet_wrap(vars(gene), ncol = 4)
p

ggsave("nhe_nt_identity_4_25_25.png", p, width = 27, height = 5)




#### hbl

hblUp = data.frame(read.csv("hblUp_identity_4_25_25.csv"))
running_mean = runner(hblUp['Identity'], k = 1, f = function(x) mean(x)) * 100
hblUp = cbind(hblUp, Mean = running_mean, gene = "hblUp") %>%
  filter(Position > max(Position)-max(nheUp$Position))

hblA = data.frame(read.csv("hblA_identity_4_25_25.csv"))
running_mean = runner(hblA['Identity'], k = 1, f = function(x) mean(x)) * 100
hblA = cbind(hblA, Mean = running_mean, gene = "hblA") 

hblB = data.frame(read.csv("hblB_identity_4_25_25.csv"))
running_mean = runner(hblB['Identity'], k = 1, f = function(x) mean(x)) * 100
hblB = cbind(hblB, Mean = running_mean, gene = "hblB") 

hblC = data.frame(read.csv("hblC_identity_4_25_25.csv"))
running_mean = runner(hblC['Identity'], k = 1, f = function(x) mean(x)) * 100
hblC = cbind(hblC, Mean = running_mean, gene = "hblC") 

hblD = data.frame(read.csv("hblD_identity_4_25_25.csv"))
running_mean = runner(hblD['Identity'], k = 1, f = function(x) mean(x)) * 100
hblD = cbind(hblD, Mean = running_mean, gene = "hblD") 

hbl_operon = rbind(hblUp, hblA, hblB, hblC, hblD)

xticks <- seq(from = 0, to = 1500, by = 300)
yticks <- c(0, 20, 40, 60, 80, 100)

hbl_operon$gene = factor(hbl_operon$gene, levels = c("hblUp", "hblC", "hblD", "hblA", "hblB"))

p = ggplot(hbl_operon, aes(x = Position)) +
  geom_area(mapping = aes(y = ifelse(Mean > 90, Mean, 0)),  fill = "gray") +
  geom_line(aes(y = Mean), color = "black", linewidth = 0.75) +
  theme_prism() +
  scale_y_continuous(expand= c(0,0), limits = c(0, 105), breaks = yticks) +
  scale_x_continuous(expand= c(0,0), breaks = xticks) +
  labs(x = "", y = "") +
  #xlab("Position (nt)") +
  #ylab("Average percent identity (%)") +
  theme(text = element_text(size = 38)) +
  facet_wrap(vars(gene), ncol = 5)
p

ggsave("hbl_nt_identity_4_25_25.png", p, width = 32, height = 5)


