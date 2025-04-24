library(ggplot2)
library(ggprism)
library(tidyverse)
library(readxl)
library(grid)

df = data.frame(read_csv("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Biomarkers paper\\SNP hits\\4_21_25\\SNP_mastersheet_4_22_25.csv"))
df$diff = abs(df$Avg.Cytotoxicity.With - df$Avg.Cytotoxicity.Without)

df_SNP = data.frame(read_excel("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Biomarkers paper\\SNP hits\\Covars 8_14_23\\SNP_hits_sheet.xlsx"))
df_gene = data.frame(read.csv("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Biomarkers paper\\gene_presence_stats_statsmodels_split.csv"))
df_gene$Gene = gsub("vir_", "", df_gene$Gene)
df_gene$diff = df_gene$Avg.Cytotoxicity.With - df_gene$Avg.Cytotoxicity.Without
df_gene$Sensitivity = df_gene$True.Positives/(df_gene$True.Positives + df_gene$False.Negatives)
df_gene$Specificity = df_gene$True.Negatives/(df_gene$True.Negatives + df_gene$False.Positives)
df_gene$NPV = df_gene$True.Negatives/(df_gene$True.Negatives + df_gene$False.Negatives)
df_gene = df_gene %>% filter(row_number() <= n()-1) %>% filter(Gene != "sph")


cytotoxicity<-df_SNP[,"Average.Cell.Viability"]

ticks<-c(-1,-0.5,0,0.5,1,1.5,2,2.5,3)

stat_box_data <- function(y, upper_limit = max(cytotoxicity) * 1.5) {
  return( 
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('n =', length(y), '\n',
                    'mean =', round(mean(y), 1), '\n')
    )
  )
}

nval <- function(y, upper_limit = max(cytotoxicity) * 1.5) {
  return( 
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('n =', length(y), '\n')
    )
  )
}



####----violin plot... uggo

df_SNP_long = pivot_longer(df_SNP, cols = colnames(df_SNP)[4:length(colnames(df_SNP))])

df_SNP_long$SNP = paste(df_SNP_long$name, df_SNP_long$value, sep = "_")
df_SNP_long = filter(df_SNP_long, value != "NA")


test = left_join(df_SNP_long, select(df, c(SNP, diff)), by = c('name' = 'SNP'))

ggplot(df_SNP_long, aes(x=SNP, y=Average.Cell.Viability, fill=value))+
  geom_violin(trim=FALSE, width = 1)+
  geom_boxplot(width=0.1,outlier.shape = NA)+
  stat_summary(fun=mean, geom="point", size=3, color="black")+
  theme_prism()+
  scale_y_continuous(breaks = ticks)+
  xlab("Group")+
  ylab("Normalized Cytotoxicity")+
  theme(text = element_text(size = 17))+
  theme(axis.text.x = element_text(size=12, angle=45, hjust = 1, vjust = 1))+
  stat_summary(fun.data = nval, geom = "text", fun.y = median, 
               position = position_dodge(width = 0.75)) 

####----accuracy and precision

df2 = pivot_longer(df, cols = c('acc', 'prec', 'sens', 'spec'))
df2$label_tile = with(df2, factor(label_tile, levels = rev(unique(df2$label_tile))))

tile = ggplot(df2, aes(x = name, y = label_tile, fill= value)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="cornflowerblue") +
  geom_text(aes(label = round(value, digits = 2)), color = "black", size = 4) +
  theme_prism()+
  theme(text = element_text(size = 16),
        axis.text.x = element_text(angle=45, hjust = 1, vjust = 1),
        legend.position = "none") +
  labs(x = "", y = "") +
  scale_x_discrete(labels = c("Accuracy", "Precision", "Sensitivity", "Specificity"))

tile

ggsave("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Figures\\Biomarkers\\sens_spec_heatmap_SNP_4_22_25.png", tile, width = 4.2, height = 9.2, units = "in")


df2_gene = pivot_longer(df_gene, cols = c('Accuracy', 'Precision', 'Sensitivity', 'Specificity'))
df2_gene$Gene = with(df2_gene, factor(Gene, levels = rev(unique(df2_gene$Gene))))

tile2 = ggplot(df2_gene, aes(x = name, y = Gene, fill= value)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="cornflowerblue") +
  geom_text(aes(label = round(value, digits = 2)), color = "black", size = 4) +
  theme_prism()+
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1), axis.text.y=element_text(face="bold.italic"))

tile2

ggsave("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Figures\\Biomarkers\\sen_spec_heatmap_gene_2_11_25.png", tile2, width = 5, height = 5, units = "in")




####----cytotoxicity

df3 = pivot_longer(df, cols = c('Avg.Cytotoxicity.With', 'Avg.Cytotoxicity.Without'))



df$label_tile = with(df,factor(label_tile,levels = unique(df$label_tile)))

bar = ggplot(df, aes(x = label_tile, y = diff, fill = diff)) + 
  geom_col(color = "black") +
  scale_fill_gradient(low="white", high="cornflowerblue") +
  scale_y_continuous(expand= c(0,0), limits = c(0,0.4)) +
  ylab("Difference in Cytotoxicity") +
  theme_prism(axis_text_angle = 45)+
  theme(axis.text.y = element_text(angle = 90, vjust = 2, hjust=0.5),
        text = element_text(size = 18),
        legend.position = "none")
bar

ggsave("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Figures\\Biomarkers\\bar_SNP_4_22_25.png", bar, width = 8, height = 4, units = "in")

df_gene$Gene = with(df_gene,factor(Gene,levels = unique(df_gene$Gene)))

bar2 = ggplot(df_gene, aes(x = Gene, y = diff, fill = diff)) + 
  geom_col(color = "black") +
  scale_fill_gradient(low="white", high="cornflowerblue") +
  scale_y_continuous(expand= c(0,0), limits = c(-0.5,0.5)) +
  ylab("Difference in Cytotoxicity") +
  theme_prism(axis_text_angle = 45)+
  theme(axis.text.y = element_text(angle = 90, vjust = 2, hjust=0.5))+
  theme(text = element_text(size = 14))

bar2

ggsave("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Figures\\Biomarkers\\bar_gene_all_negs_1_11_24.png", bar2, width = 5, height = 4, units = "in")



png(file = "C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Figures\\Biomarkers\\table_gene.png", units = "in", width = 2, height = 4, res = 600)

grid.draw(tableGrob(df_gene[,6:7], rows = NULL, cols = NULL))

dev.off()

######Hydrophobicity


df$label_aa = with(df,factor(label_aa,levels = unique(df$label_aa)))
yticks = c(-60, -30, 0, 30, 60)

p = ggplot(df, aes(x = label_aa, y = diff_hp, fill = diff_hp)) + 
  geom_col(color = "black") + 
  scale_fill_gradient(low="white", high="cornflowerblue") +
  scale_y_continuous(breaks = yticks) +
  labs(y = "Difference in Hydrophobicity", x= "") +
  theme_prism() +
  theme_prism(axis_text_angle = 45) +
  theme(axis.text.y = element_text(angle = 90, vjust = 2, hjust=0.5),
        text = element_text(size = 18),
        legend.position = "none") 

p

ggsave("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Figures\\Biomarkers\\snp_hydrophobicity_4_23_25.png", p, dpi = 600, width = 4, height = 4, units = "in")


