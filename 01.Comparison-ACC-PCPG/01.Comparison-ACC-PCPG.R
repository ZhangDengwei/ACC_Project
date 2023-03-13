# load libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(vegan)
library(UpSetR)
library(ggsci)
library(doBy)


# load metadata 
df_metadata <- read.csv("01.Data/Metadata-TCGA-Kraken-17625-Samples.csv")
df_metadata_ACC <- df_metadata %>% filter(investigation=="TCGA-ACC")
df_metadata_PCPG <- df_metadata %>% filter(investigation=="TCGA-PCPG")

df_metadata_ACC <- df_metadata %>% filter(investigation=="TCGA-ACC")
# remove two patients missing pathologic_stage_label
df_metadata_ACC_f <- df_metadata_ACC %>% filter(pathologic_stage_label != "Not available")
# meta data from TCGA web
df_meta_ACC_2 <- read.delim("01.Data/ACC_77_meta.txt")
df_metadata_ACC_final <- merge(df_metadata_ACC_f, df_meta_ACC_2[,c("id", setdiff(names(df_meta_ACC_2), names(df_metadata_ACC_f)))], by="id", all.x=TRUE)
df_metadata_ACC_final$CENSOR <- ifelse(df_metadata_ACC_final$vital_status=="Alive", 0, 1)
df_metadata_ACC_final$race_ethnicity <- ifelse(df_metadata_ACC_final$race=="WHITE", "WHITE", "Not_WHITE")

write.csv(df_metadata_ACC_final, file = "Tables/metadata_77_ACC.csv")

# load microbial data from Poore's study
df_raw <- read.csv("01.Data/Kraken-TCGA-Raw-Data-17625-Samples.csv")
df_normalization <- read.csv("01.Data/Kraken-TCGA-Voom-SNM-Full-Data.csv")
df_normalization_filter_likely <- read.csv("01.Data/Kraken-TCGA-Voom-SNM-Likely-Contaminants-Removed-Data.csv")
df_normalization_filter_Plate_Center <- read.csv("01.Data/Kraken-TCGA-Voom-SNM-Plate-Center-Filtering-Data.csv")
df_normalization_filter_putative <- read.csv("01.Data/Kraken-TCGA-Voom-SNM-All-Putative-Contaminants-Removed-Data.csv")
df_normalization_filter_stringent <- read.csv("01.Data/Kraken-TCGA-Voom-SNM-Most-Stringent-Filtering-Data.csv")

# microbial for ACC and PCPG, where only include "Primary Tumor"
meta_ACC_PCPG <- df_metadata %>% filter(investigation %in% c("TCGA-ACC", "TCGA-PCPG")) %>% filter(sample_type=="Primary Tumor")
samples_ACC_PCPG <- meta_ACC_PCPG$id

df_raw_ACC_PCPG <- df_raw %>% filter(X %in% samples_ACC_PCPG)%>%
                                remove_rownames %>% column_to_rownames(var="X")
df_normalization_ACC_PCPG  <- df_normalization %>% filter(X %in% samples_ACC_PCPG) %>%
                                remove_rownames %>% column_to_rownames(var="X")
df_normalization_filter_likely_ACC_PCPG <- df_normalization_filter_likely %>% filter(X %in% samples_ACC_PCPG)%>%
                                              remove_rownames %>% column_to_rownames(var="X")
df_normalization_filter_Plate_Center_ACC_PCPG <- df_normalization_filter_Plate_Center %>% filter(X %in% samples_ACC_PCPG)%>%
                                                  remove_rownames %>% column_to_rownames(var="X")
df_normalization_filter_putative_ACC_PCPG <- df_normalization_filter_putative %>% filter(X %in% samples_ACC_PCPG)%>%
                                                remove_rownames %>% column_to_rownames(var="X")
df_normalization_filter_stringent_ACC_PCPG <- df_normalization_filter_stringent %>% filter(X %in% samples_ACC_PCPG)%>%
                                                remove_rownames %>% column_to_rownames(var="X")

# virus for ACC and PCPG
df_raw_ACC_PCPG_Vir <- df_raw_ACC_PCPG %>% select(grep("Virus", names(df_raw_ACC_PCPG)))
df_normalization_ACC_PCPG_Vir <- df_normalization_ACC_PCPG %>% select(grep("Virus", names(df_normalization_ACC_PCPG)))
df_normalization_filter_likely_ACC_PCPG_Vir <- df_normalization_filter_likely_ACC_PCPG %>% select(grep("Virus", names(df_normalization_filter_likely_ACC_PCPG)))
df_normalization_filter_Plate_Center_ACC_PCPG_Vir <- df_normalization_filter_Plate_Center_ACC_PCPG %>% select(grep("Virus", names(df_normalization_filter_Plate_Center_ACC_PCPG)))
df_normalization_filter_putative_ACC_PCPG_Vir <- df_normalization_filter_putative_ACC_PCPG %>% select(grep("Virus", names(df_normalization_filter_putative_ACC_PCPG)))
df_normalization_filter_stringent_ACC_PCPG_Vir <- df_normalization_filter_stringent_ACC_PCPG %>% select(grep("Virus", names(df_normalization_filter_stringent_ACC_PCPG)))

# archaea for ACC and PCPG
df_raw_ACC_PCPG_Arc <- df_raw_ACC_PCPG %>% select(grep("Archaea", names(df_raw_ACC_PCPG)))
df_normalization_ACC_PCPG_Arc <- df_normalization_ACC_PCPG %>% select(grep("Archaea", names(df_normalization_ACC_PCPG)))
df_normalization_filter_likely_ACC_PCPG_Arc <- df_normalization_filter_likely_ACC_PCPG %>% select(grep("Archaea", names(df_normalization_filter_likely_ACC_PCPG)))
df_normalization_filter_Plate_Center_ACC_PCPG_Arc <- df_normalization_filter_Plate_Center_ACC_PCPG %>% select(grep("Archaea", names(df_normalization_filter_Plate_Center_ACC_PCPG)))
df_normalization_filter_putative_ACC_PCPG_Arc <- df_normalization_filter_putative_ACC_PCPG %>% select(grep("Archaea", names(df_normalization_filter_putative_ACC_PCPG)))
df_normalization_filter_stringent_ACC_PCPG_Arc <- df_normalization_filter_stringent_ACC_PCPG %>% select(grep("Archaea", names(df_normalization_filter_stringent_ACC_PCPG)))

# bacteria for ACC and PCPG
df_raw_ACC_PCPG_Bac <- df_raw_ACC_PCPG %>% select(grep("Bacteria", names(df_raw_ACC_PCPG)))
df_normalization_ACC_PCPG_Bac <- df_normalization_ACC_PCPG %>% select(grep("Bacteria", names(df_normalization_ACC_PCPG)))
df_normalization_filter_likely_ACC_PCPG_Bac <- df_normalization_filter_likely_ACC_PCPG %>% select(grep("Bacteria", names(df_normalization_filter_likely_ACC_PCPG)))
df_normalization_filter_Plate_Center_ACC_PCPG_Bac <- df_normalization_filter_Plate_Center_ACC_PCPG %>% select(grep("Bacteria", names(df_normalization_filter_Plate_Center_ACC_PCPG)))
df_normalization_filter_putative_ACC_PCPG_Bac <- df_normalization_filter_putative_ACC_PCPG %>% select(grep("Bacteria", names(df_normalization_filter_putative_ACC_PCPG)))
df_normalization_filter_stringent_ACC_PCPG_Bac <- df_normalization_filter_stringent_ACC_PCPG %>% select(grep("Bacteria", names(df_normalization_filter_stringent_ACC_PCPG)))

write.csv(df_raw_ACC_PCPG, file = "Tables/Raw-ACC_PCPG.csv", row.names = T)
write.csv(df_normalization_ACC_PCPG, file = "Tables/Voom-SNM-ACC_PCPG.csv", row.names = T)
write.csv(df_normalization_filter_likely_ACC_PCPG, file = "Tables/Voom-SNM-Filter-Likely-ACC_PCPG.csv", row.names = T)
write.csv(df_normalization_filter_Plate_Center_ACC_PCPG, file = "Tables/Voom-SNM-Filter-Plate_Center-ACC_PCPG.csv", row.names = T)
write.csv(df_normalization_filter_putative_ACC_PCPG, file = "Tables/Voom-SNM-Filter-Putative-ACC_PCPG.csv", row.names = T)
write.csv(df_normalization_filter_stringent_ACC_PCPG, file = "Tables/Voom-SNM-Filter-Stringent-ACC_PCPG.csv", row.names = T)

df_raw_ACC <- df_raw_ACC_PCPG[df_metadata_ACC_final$id, ]
df_normalization_ACC <- df_normalization_ACC_PCPG[df_metadata_ACC_final$id, ]
df_normalization_filter_likely_ACC <- df_normalization_filter_likely_ACC_PCPG[df_metadata_ACC_final$id, ]
df_normalization_filter_Plate_Center_ACC <- df_normalization_filter_Plate_Center_ACC_PCPG[df_metadata_ACC_final$id, ]
df_normalization_filter_putative_ACC <- df_normalization_filter_putative_ACC_PCPG[df_metadata_ACC_final$id, ]
df_normalization_filter_stringent_ACC <- df_normalization_filter_stringent_ACC_PCPG[df_metadata_ACC_final$id, ]

write.csv(df_raw_ACC, file = "Tables/Raw-ACC.csv", row.names = T)
write.csv(df_normalization_ACC, file = "Tables/Voom-SNM-ACC.csv", row.names = T)
write.csv(df_normalization_filter_likely_ACC, file = "Tables/Voom-SNM-Filter-Likely-ACC.csv", row.names = T)
write.csv(df_normalization_filter_Plate_Center_ACC, file = "Tables/Voom-SNM-Filter-Plate_Center-ACC.csv", row.names = T)
write.csv(df_normalization_filter_putative_ACC, file = "Tables/Voom-SNM-Filter-Putative-ACC.csv", row.names = T)
write.csv(df_normalization_filter_stringent_ACC, file = "Tables/Voom-SNM-Filter-Stringent-ACC.csv", row.names = T)



######################################################################
# Examining difference between ACC and PCPG using PCoA
######################################################################
function_PcoA <- function(in_matrix, title){
  in_matrix[in_matrix<0] <- 0 
  distance <- vegdist(in_matrix, method = 'bray')
  pcoa <- cmdscale(distance, k = (nrow(in_matrix) - 1), eig = TRUE)
  point <- data.frame(pcoa$point)
  pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
  # extract first two coordinate value
  sample_eig <- data.frame({pcoa$point})[1:2]
  sample_eig$id <- rownames(sample_eig)
  names(sample_eig)[1:2] <- c('PCoA1', 'PCoA2')
  group_m <- merge(sample_eig, df_metadata[,c("id", "investigation")], by = "id", all.x = TRUE)
  # PERMANOVA
  dt <- in_matrix[group_m$id,]
  print(identical(rownames(dt), group_m$id))
  
  group_m[,"investigation"] <- as.factor(group_m[,"investigation"])
  adonis_result <- adonis2(dt~investigation, group_m, permutations = 999, distance = 'bray')
  
  R2 <- round(adonis_result$R2[1], 2)
  pvalue <- round(adonis_result$`Pr(>F)`[1], 2)
  plot <- ggscatter(group_m, x= "PCoA1", y = "PCoA2",color="investigation",
                    ellipse = TRUE,
                    mean.point = TRUE, star.plot = TRUE,
                    ellipse.level = 0.95,
                    ggtheme = theme_minimal()) +
    labs(x = paste('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'),
         y = paste('PCoA2: ', round(100 * pcoa_eig[2], 2), '%'))+
    theme_classic()+
    scale_color_manual(values = c("#0072B5CC","#E18727CC"))+
    geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
    geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
    theme(panel.grid = element_line(color = 'black', linetype = 2, size = 0.1), 
          panel.background = element_rect(color = 'black', fill = 'transparent'), 
          legend.title=element_blank())+
    theme(axis.title = element_text(size = 18, colour="black"),
          axis.text = element_text(size = 16, colour = "black"),
          legend.text = element_text(size = 16))+
    geom_text(x=0, y=max(group_m$PCoA2)*0.8, label=paste("PERMANOVA R2 = ",R2 ,"\n","P = ", pvalue, sep = ""))+
    ggtitle(title)
  return(plot)
}


# PCA examination, only for primary tumor
pca_nor <- function_PcoA(df_normalization_ACC_PCPG, "Overall-Norm")
pca_nor_likely <- function_PcoA(df_normalization_filter_likely_ACC_PCPG, "Overall-Norm-Likely")
pca_nor_Plate_Center <- function_PcoA(df_normalization_filter_Plate_Center_ACC_PCPG, "Overall-Norm-Plate_Center")
pca_nor_putative <- function_PcoA(df_normalization_filter_putative_ACC_PCPG, "Overall-Norm-Putative")
pca_nor_stringent <- function_PcoA(df_normalization_filter_stringent_ACC_PCPG, "Overall-Norm-Stringent")

pca_nor_Vir <- function_PcoA(df_normalization_ACC_PCPG_Vir, "Virus-Norm")
pca_nor_likely_Vir  <- function_PcoA(df_normalization_filter_likely_ACC_PCPG_Vir, "Virus-Likely")
pca_nor_Plate_Center_Vir  <- function_PcoA(df_normalization_filter_Plate_Center_ACC_PCPG_Vir, "Virus-Plate_Center")
pca_nor_putative_Vir  <- function_PcoA(df_normalization_filter_putative_ACC_PCPG_Vir, "Virus-Putative")
pca_nor_stringent_Vir  <- function_PcoA(df_normalization_filter_stringent_ACC_PCPG_Vir, "Virus-Stringent")

pca_nor_Arc <- function_PcoA(df_normalization_ACC_PCPG_Arc, "Archaea-Norm")
pca_nor_likely_Arc  <- function_PcoA(df_normalization_filter_likely_ACC_PCPG_Arc, "Archaea-Likely")
pca_nor_Plate_Center_Arc  <- function_PcoA(df_normalization_filter_Plate_Center_ACC_PCPG_Arc, "Archaea-Plate_Center")
pca_nor_putative_Arc  <- function_PcoA(df_normalization_filter_putative_ACC_PCPG_Arc, "Archaea-Putative")
pca_nor_stringent_Arc  <- function_PcoA(df_normalization_filter_stringent_ACC_PCPG_Arc, "Archaea-Stringent")

pca_nor_Bac <- function_PcoA(df_normalization_ACC_PCPG_Bac, "Bacteria-Norm")
pca_nor_likely_Bac  <- function_PcoA(df_normalization_filter_likely_ACC_PCPG_Bac, "Bacteria-Likely")
pca_nor_Plate_Center_Bac  <- function_PcoA(df_normalization_filter_Plate_Center_ACC_PCPG_Bac, "Bacteria-Norm-Plate_Center")
pca_nor_putative_Bac  <- function_PcoA(df_normalization_filter_putative_ACC_PCPG_Bac, "Bacteria-Norm-Putative")
pca_nor_stringent_Bac  <- function_PcoA(df_normalization_filter_stringent_ACC_PCPG_Bac, "Bacteria-Norm-Stringent")

pdf(file = "Figures/Comparison_PCA_ACC-PCPG.pdf", width = 7, height = 5)
pca_nor
pca_nor_likely
pca_nor_Plate_Center
pca_nor_putative
pca_nor_stringent
pca_nor_Vir
pca_nor_likely_Vir
pca_nor_Plate_Center_Vir
pca_nor_putative_Vir
pca_nor_stringent_Vir
pca_nor_Arc
pca_nor_likely_Arc
pca_nor_Plate_Center_Arc
pca_nor_putative_Arc
pca_nor_stringent_Arc
pca_nor_Bac
pca_nor_likely_Bac
pca_nor_Plate_Center_Bac
pca_nor_putative_Bac
pca_nor_stringent_Bac
dev.off()




######################################################################
# Intersection of microbes in ACC
######################################################################
microbe_without <- names(df_normalization_ACC)[grep("Virus|Bacteria|Archae", names(df_normalization_ACC))]
microbe_likely <- names(df_normalization_filter_likely_ACC)[grep("Virus|Bacteria|Archae", names(df_normalization_filter_likely_ACC))]
microbe_PC <- names(df_normalization_filter_Plate_Center_ACC)[grep("Virus|Bacteria|Archae", names(df_normalization_filter_Plate_Center_ACC))]
microbe_putative <- names(df_normalization_filter_putative_ACC)[grep("Virus|Bacteria|Archae", names(df_normalization_filter_putative_ACC))]
microbe_stringent <- names(df_normalization_filter_stringent_ACC)[grep("Virus|Bacteria|Archae", names(df_normalization_filter_stringent_ACC))]

df_without_all <- data.frame(microbe=microbe_without,without=1)
df_without_Vir <- data.frame(virus=microbe_without[grep("Virus",microbe_without)],
                             without=1)
df_without_Arc <- data.frame(Archae=microbe_without[grep("Archae",microbe_without)],
                             without=1)
df_without_Bac <- data.frame(Bacteria=microbe_without[grep("Bacteria",microbe_without)],
                             without=1)
df_likely_all <- data.frame(microbe=microbe_likely,likely=1)
df_likely_Vir <- data.frame(virus=microbe_likely[grep("Virus",microbe_likely)],
                             likely=1)
df_likely_Arc <- data.frame(Archae=microbe_likely[grep("Archae",microbe_likely)],
                             likely=1)
df_likely_Bac <- data.frame(Bacteria=microbe_likely[grep("Bacteria",microbe_likely)],
                             likely=1)
df_PC_all <- data.frame(microbe=microbe_PC,PC=1)
df_PC_Vir <- data.frame(virus=microbe_PC[grep("Virus",microbe_PC)],
                            PC=1)
df_PC_Arc <- data.frame(Archae=microbe_PC[grep("Archae",microbe_PC)],
                            PC=1)
df_PC_Bac <- data.frame(Bacteria=microbe_PC[grep("Bacteria",microbe_PC)],
                            PC=1)
df_putative_all <- data.frame(microbe=microbe_putative,putative=1)
df_putative_Vir <- data.frame(virus=microbe_putative[grep("Virus",microbe_putative)],
                            putative=1)
df_putative_Arc <- data.frame(Archae=microbe_putative[grep("Archae",microbe_putative)],
                            putative=1)
df_putative_Bac <- data.frame(Bacteria=microbe_putative[grep("Bacteria",microbe_putative)],
                            putative=1)
df_stringent_all <- data.frame(microbe=microbe_stringent,stringent=1)
df_stringent_Vir <- data.frame(virus=microbe_stringent[grep("Virus",microbe_stringent)],
                            stringent=1)
df_stringent_Arc <- data.frame(Archae=microbe_stringent[grep("Archae",microbe_stringent)],
                            stringent=1)
df_stringent_Bac <- data.frame(Bacteria=microbe_stringent[grep("Bacteria",microbe_stringent)],
                            stringent=1)

df_combine_all <- Reduce(function(x,y){merge(x,y,by="microbe", all=TRUE)},
                         list(df_without_all, df_likely_all, df_PC_all, df_putative_all, df_stringent_all))
df_combine_Vir <- Reduce(function(x,y){merge(x,y,by="virus", all=TRUE)},
                         list(df_without_Vir, df_likely_Vir, df_PC_Vir, df_putative_Vir, df_stringent_Vir))
df_combine_Arc <- Reduce(function(x,y){merge(x,y,by="Archae", all=TRUE)},
                         list(df_without_Arc, df_likely_Arc, df_PC_Arc, df_putative_Arc, df_stringent_Arc))
df_combine_Bac <- Reduce(function(x,y){merge(x,y,by="Bacteria", all=TRUE)},
                         list(df_without_Bac, df_likely_Bac, df_PC_Bac, df_putative_Bac, df_stringent_Bac))

df_combine_all[is.na(df_combine_all)] <- 0
df_combine_Vir[is.na(df_combine_Vir)] <- 0
df_combine_Arc[is.na(df_combine_Arc)] <- 0
df_combine_Bac[is.na(df_combine_Bac)] <- 0
df_combine_all[df_combine_all!=0] <- 1
df_combine_Vir[df_combine_Vir!=0] <- 1
df_combine_Arc[df_combine_Arc!=0] <- 1
df_combine_Bac[df_combine_Bac!=0] <- 1

plot_upset_All <- upset(df_combine_all[,2:6], 
                        keep.order=TRUE, sets.x.label = "Number of microbes")
plot_upset_Vir <- upset(df_combine_Vir[,2:6], 
                    keep.order=TRUE, sets.x.label = "Number of virus")
plot_upset_Arc <- upset(df_combine_Arc[,2:6], 
                        keep.order=TRUE, sets.x.label = "Number of archae")
plot_upset_Bac <- upset(df_combine_Bac[,2:6], 
                        keep.order=TRUE, sets.x.label = "Number of bacteria")

pdf(file = "Figures/Intersection_of_microbes_at_different_levels_of_filter.pdf",height = 4,width = 5)
plot_upset_All
plot_upset_Vir
plot_upset_Arc
plot_upset_Bac
dev.off()



######################################################################
# Relative abundance of microbes in ACC
######################################################################
df_raw_ACC <- df_raw_ACC_PCPG[df_metadata_ACC_final$id, ]
df_normalization_ACC <- df_normalization_ACC_PCPG[df_metadata_ACC_final$id, ]
df_normalization_filter_likely_ACC <- df_normalization_filter_likely_ACC_PCPG[df_metadata_ACC_final$id, ]
df_normalization_filter_Plate_Center_ACC <- df_normalization_filter_Plate_Center_ACC_PCPG[df_metadata_ACC_final$id, ]
df_normalization_filter_putative_ACC <- df_normalization_filter_putative_ACC_PCPG[df_metadata_ACC_final$id, ]
df_normalization_filter_stringent_ACC <- df_normalization_filter_stringent_ACC_PCPG[df_metadata_ACC_final$id, ]

df_without_raw <- df_raw_ACC[,intersect(names(df_raw_ACC), names(df_normalization_ACC))]
df_likely_raw <- df_raw_ACC[,intersect(names(df_raw_ACC), names(df_normalization_filter_likely_ACC))]
df_PC_raw <- df_raw_ACC[,intersect(names(df_raw_ACC), names(df_normalization_filter_Plate_Center_ACC))]
df_putative_raw <- df_raw_ACC[,intersect(names(df_raw_ACC), names(df_normalization_filter_putative_ACC))]
df_stringent_raw <- df_raw_ACC[,intersect(names(df_raw_ACC), names(df_normalization_filter_stringent_ACC))]

df_without_raw <- df_without_raw[,colSums(df_without_raw)>0]
df_likely_raw <- df_likely_raw[,colSums(df_likely_raw)>0]
df_PC_raw <- df_PC_raw[,colSums(df_PC_raw)>0]
df_putative_raw <- df_putative_raw[,colSums(df_putative_raw)>0]
df_stringent_raw <- df_stringent_raw[,colSums(df_stringent_raw)>0]

df_without_raw_rb <- df_without_raw/rowSums(df_without_raw)
df_likely_raw_rb <- df_likely_raw/rowSums(df_likely_raw)
df_PC_raw_rb <- df_PC_raw/rowSums(df_PC_raw)
df_putative_raw_rb <- df_putative_raw/rowSums(df_putative_raw)
df_stringent_raw_rb <- df_stringent_raw/rowSums(df_stringent_raw)

names(df_without_raw_rb) <- str_extract(names(df_without_raw_rb), "g__.+")
names(df_likely_raw_rb) <- str_extract(names(df_likely_raw_rb), "g__.+")
names(df_PC_raw_rb) <- str_extract(names(df_PC_raw_rb), "g__.+")
names(df_putative_raw_rb) <- str_extract(names(df_putative_raw_rb), "g__.+")
names(df_stringent_raw_rb) <- str_extract(names(df_stringent_raw_rb), "g__.+")

function_top <- function(in_matrix, identifier, top_no){
                  colsum <- colSums(in_matrix)
                  colsum <- colsum[order(colsum, decreasing = TRUE)]
                  top10 <- head(names(colsum), top_no)
                  df_top10 <- in_matrix[,top10]
                  df_others <- rowSums(in_matrix[,! names(in_matrix) %in% top10])
                  df_top10$Others <- df_others
                  df_top10$id <- rownames(df_top10)
                  df_top10_melt <- reshape2::melt(df_top10, id="id")
                  df_top10_melt[,"identifier"] <- identifier
                  return(df_top10_melt)
}

df_without_raw_rb_top <- function_top(df_without_raw_rb, "without", 5)
df_likely_raw_rb_top <- function_top(df_likely_raw_rb, "likely", 5)
df_PC_raw_rb_top <- function_top(df_PC_raw_rb, "PC", 5)
df_putative_raw_rb_top <- function_top(df_putative_raw_rb, "putative", 5)
df_stringent_raw_rb_top <- function_top(df_stringent_raw_rb, "stringent", 5)

df_combine <- rbind(df_without_raw_rb_top, df_likely_raw_rb_top,
                    df_PC_raw_rb_top, df_putative_raw_rb_top, df_stringent_raw_rb_top)
df_combine$identifier <- factor(df_combine$identifier, levels = c("without", "likely", "PC", "putative", "stringent"))

df_order <- doBy::summaryBy(value~variable, data=df_combine[,c("variable", "value")], FUN=sum)
df_order <- df_order[order(df_order$value.sum),]
df_combine$variable <- factor(df_combine$variable, levels = df_order$variable)

p_top <- ggplot(df_combine, aes(x=id,  y=value)) + 
          geom_area(aes(fill = variable, group = variable))+
          facet_wrap(vars(identifier), nrow = 5) +
          scale_fill_d3("category20")+
          theme_bw()+
          xlab("")+
          ylab("Relative abundance")+
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank())
ggsave(filename = "Figures/Relative_abundance_different_levels_filter.pdf",plot = p_top, useDingbats = FALSE, width = 12,height = 5)
