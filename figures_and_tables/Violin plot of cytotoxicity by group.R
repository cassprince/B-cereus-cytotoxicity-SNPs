library(ggplot2)
library(readxl)
library(ggprism)
library(tidyverse)
library(stringr)
library(paletteer)
library(ggsignif)

setwd("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Biomarkers paper")
data<-data.frame(read_excel("Mastersheet_no_clones.xlsx"))

##Clean species names.
species<-data[,"final.taxon.names"]
species <- gsub(";.*", "", species)
species <- recode(species, "B. mosaicus subsp. cereus biovar Emeticus*" = "B. mosaicus subsp. cereus biovar Emeticus")
biovar <- str_extract(species, "biovar .*")
biovar <- gsub("biovar ", "", biovar)
spp <- str_extract(data[,"BTyper3_Species"], "^[^\\(]+")
spp = gsub("\\*", "", data$spp)

group<-data[,"Adjusted.panC.Group..predicted.species."]
group<-gsub("\\s*\\([^\\)]+\\)", "", group)
group<-gsub("_", " ", group)
cytotoxicity<-data[,"Average.Cell.Viability....0.7.is.cytotoxic."]

df <- data.frame(cbind(group = group, species = species, biovar = biovar, spp = spp, cytotoxicity = cytotoxicity))
df_short <- df[!(df$species %in% c("B. clarus", "B. luti", "B. mosaicus biovar Emeticus", "B. paramycoides", "B. toyonensis biovar Thuringiensis")),]


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

##Violin plot with jitter points overlayed.
plt_jitter <- ggplot(data, aes(x=group, y=cytotoxicity))+
  geom_violin(trim=FALSE, fill="gray", width = 1.7)+
  geom_jitter(width=0.05)+
  theme_prism()+
  scale_y_continuous(breaks = ticks)+
  xlab("Group")+
  ylab("Normalized Cytotoxicity")+
  theme(text = element_text(size = 17))+
  theme(axis.text.x = element_text(size=12, angle=45))+ 
  stat_summary(fun.data = stat_box_data, geom = "text", fun.y = median, 
               position = position_dodge(width = 0.75))

ggsave("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Figures\\Biomarkers\\Cytotoxicity by group violin plot jitter reversed.png", plot = plt_jitter, device = "png", width = 11, height = 9, dpi = 300)



#test = data.matrix(crossing(df$group, df$group))

#test2 = split(test, seq(nrow(test)))

##Violin plot with box plot overlayed.
plt_box <- ggplot(data, aes(x=group, y=cytotoxicity))+
  geom_violin(trim=FALSE, width = 1.7, aes(fill = group))+
  geom_boxplot(width=0.1, outlier.shape = NA)+
  stat_summary(fun=mean, geom="point", size=3, color="black")+
  theme_prism()+
  scale_y_continuous(breaks = ticks)+
  scale_fill_manual(values = c("#e5c9ff", "gray", "gray", "#e5c9ff", "#e5c9ff","gray", "gray", "gray" ))+
  xlab("")+
  ylab("Normalized Cytotoxicity")+
  theme(text = element_text(size = 15))+
  theme(axis.text.x = element_text(size=15, angle=45, hjust = 1, vjust = 1))+
  stat_summary(fun.data = nval, geom = "text", fun.y = median, 
               position = position_dodge(width = 0.75)) 

plt_box

ggsave("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Figures\\Biomarkers\\Cytotoxicity group violin plot 9_6_24.png", plot = plt_box, device = "png", width = 11, height = 6, dpi = 300)


##Violin plot by species. Need to clean up species names to make it prettier.
plt_spp <- ggplot(df_short, aes(x=as.character(species), y=as.numeric(cytotoxicity),  fill=spp))+
  geom_violin(trim=FALSE, width = 2)+
  geom_boxplot(width=0.1)+
  stat_summary(fun=mean, geom="point", size=3, color="black")+
  theme_prism()+
  scale_y_continuous(breaks = ticks)+
  xlab("Group")+
  ylab("Normalized Cytotoxicity")+
  theme(text = element_text(size = 17))+
  theme(axis.text.x = element_text(size=12, angle=45, hjust = 1, vjust = 1))+
  stat_summary(fun.data = nval, geom = "text", fun.y = median)+
  scale_fill_brewer(palette = "Pastel1")
plt_spp
ggsave("C:\\Users\\cassp\\OneDrive\\Documents\\Kovac Lab\\Figures\\Biomarkers\\Cytotoxicity by species violin plot box reversed.png", plot = plt_spp, device = "png", width = 14, height = 10, dpi = 300)






#########################################


mean(data$Average.Cell.Viability.Reversed[data$Adjusted.panC.Group..predicted.species. == "Group_I(pseudomycoides)"])
sd(data$Average.Cell.Viability.Reversed[data$Adjusted.panC.Group..predicted.species. == "Group_I(pseudomycoides)"])
nrow(data[data$Adjusted.panC.Group..predicted.species. == "Group_I(pseudomycoides)",])


mean(data$Average.Cell.Viability.Reversed[data$Adjusted.panC.Group..predicted.species. == "Group_II(mosaicus/luti)"])
sd(data$Average.Cell.Viability.Reversed[data$Adjusted.panC.Group..predicted.species. == "Group_II(mosaicus/luti)"])

mean(data$Average.Cell.Viability.Reversed[data$Adjusted.panC.Group..predicted.species. == "Group_III(mosaicus)"])
sd(data$Average.Cell.Viability.Reversed[data$Adjusted.panC.Group..predicted.species. == "Group_III(mosaicus)"])

mean(data$Average.Cell.Viability.Reversed[data$Adjusted.panC.Group..predicted.species. == "Group_IV(cereus_sensu_stricto)"])
sd(data$Average.Cell.Viability.Reversed[data$Adjusted.panC.Group..predicted.species. == "Group_IV(cereus_sensu_stricto)"])

mean(data$Average.Cell.Viability.Reversed[data$Adjusted.panC.Group..predicted.species. == "Group_V(toyonensis)"])
sd(data$Average.Cell.Viability.Reversed[data$Adjusted.panC.Group..predicted.species. == "Group_V(toyonensis)"])

mean(data$Average.Cell.Viability.Reversed[data$Adjusted.panC.Group..predicted.species. == "Group_VI(mycoides/paramycoides)"])
sd(data$Average.Cell.Viability.Reversed[data$Adjusted.panC.Group..predicted.species. == "Group_VI(mycoides/paramycoides)"])

mean(data$Average.Cell.Viability.Reversed[data$Adjusted.panC.Group..predicted.species. == "Group_VII(cytotoxicus)"])
sd(data$Average.Cell.Viability.Reversed[data$Adjusted.panC.Group..predicted.species. == "Group_VII(cytotoxicus)"])

mean(data$Average.Cell.Viability.Reversed[data$Adjusted.panC.Group..predicted.species. == "Group_VIII(mycoides)"])
sd(data$Average.Cell.Viability.Reversed[data$Adjusted.panC.Group..predicted.species. == "Group_VIII(mycoides)"])

