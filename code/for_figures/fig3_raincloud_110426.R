library(tidyverse)
library(ggdist)
library(ggprism)

setwd("C:/Users/princecr/Documents/GitHub/B-cereus-cytotoxicity-SNPs/data")

df = read_csv("df_all_genes_090126.csv") %>%
  mutate(panC_group = gsub("_", " ", panC_group))

df_summ = df %>%
  group_by(panC_group) %>%
  summarize(n = n())

df_ast = data.frame(panC_group = df_summ$panC_group, asterisk = c())

df_ast = data.frame(panC_group =c("Group I", "Group IV", "Group V"))

### Raincloud plot 

ggplot(df, aes(x = panC_group, y = cytotoxicity)) +
  stat_halfeye(adjust = 1.5, width = 0.5, .width = 0, justification = -.5, point_colour = NA) +
  geom_boxplot(width = 0.3, outlier.color = NA, aes(fill = panC_group)) +
  geom_jitter(width = 0.1, alpha = 0.2, size = 1)+
  geom_text(data = df_summ, aes(x = panC_group, y = 1.5, label = paste0("n = ", n)), size = 4.5) +
  geom_text(data = df_ast, aes(x = panC_group, y = 1.37, label = "*"), size = 5) +
  theme_classic() +
  labs(y = "Normalized Cytotoxicity", x = "") +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme(axis.text = element_text(color = "black", size = 13),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(size = 14),
        panel.grid.major.y = element_line(color  = "gray75"),
        panel.grid.minor.y = element_line(color  = "gray90")) +
  scale_fill_manual(values = c("white", "#fc8383", "#fcc983","#fcf083", "#b3fc83", "#b5dcff", "#a0a5f7", "#c1a0f7", "#f7a0e7"))

ggsave("../figures/fig3_110926.png", width = 7, height = 5, units = "in", dpi = 350)
  
  



