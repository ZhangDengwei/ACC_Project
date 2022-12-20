# load libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(vegan)
library(phyloseq)
library(ggpubr)
library(survival)
library(survminer)
library(ggthemes)


# load metadata
df_meta_ACC <- read.csv("Tables/metadata_77_ACC.csv")

# load data
df_raw <- read.csv("01.Data/Kraken-TCGA-Raw-Data-17625-Samples.csv")
df_normalization <- read.csv("01.Data/Kraken-TCGA-Voom-SNM-Full-Data.csv")
df_normalization_filter_likely <- read.csv("01.Data/Kraken-TCGA-Voom-SNM-Likely-Contaminants-Removed-Data.csv")
df_normalization_filter_Plate_Center <- read.csv("01.Data/Kraken-TCGA-Voom-SNM-Plate-Center-Filtering-Data.csv")
df_normalization_filter_putative <- read.csv("01.Data/Kraken-TCGA-Voom-SNM-All-Putative-Contaminants-Removed-Data.csv")
df_normalization_filter_stringent <- read.csv("01.Data/Kraken-TCGA-Voom-SNM-Most-Stringent-Filtering-Data.csv")

df_raw_ACC <- df_raw %>% filter(X %in% df_meta_ACC$id) %>% remove_rownames %>% column_to_rownames(var="X")
df_normalization_ACC <- df_normalization %>% filter(X %in% df_meta_ACC$id) %>% remove_rownames %>% column_to_rownames(var="X")
df_normalization_filter_likely_ACC <- df_normalization_filter_likely %>% filter(X %in% df_meta_ACC$id) %>% remove_rownames %>% column_to_rownames(var="X")
df_normalization_filter_Plate_Center_ACC <- df_normalization_filter_Plate_Center %>% filter(X %in% df_meta_ACC$id) %>% remove_rownames %>% column_to_rownames(var="X")
df_normalization_filter_putative_ACC <- df_normalization_filter_putative %>% filter(X %in% df_meta_ACC$id) %>% remove_rownames %>% column_to_rownames(var="X")
df_normalization_filter_stringent_ACC <- df_normalization_filter_stringent %>% filter(X %in% df_meta_ACC$id) %>% remove_rownames %>% column_to_rownames(var="X")

# raw counts after contamination removal
df_raw_ACC_filter_likely <- df_raw_ACC[,intersect(names(df_normalization_filter_likely_ACC), names(df_raw_ACC))]
df_raw_ACC_filter_Plate_Center <- df_raw_ACC[,intersect(names(df_normalization_filter_Plate_Center_ACC), names(df_raw_ACC))]
df_raw_ACC_filter_putative <- df_raw_ACC[,intersect(names(df_normalization_filter_putative_ACC), names(df_raw_ACC))]
df_raw_ACC_filter_stringent <- df_raw_ACC[,intersect(names(df_normalization_filter_stringent_ACC), names(df_raw_ACC))]

# remove empty cols and two patients missing stage information
df_raw_ACC <- df_raw_ACC[df_meta_ACC$id,colSums(df_raw_ACC)>0]
df_raw_ACC_filter_likely <- df_raw_ACC_filter_likely[df_meta_ACC$id, colSums(df_raw_ACC_filter_likely)>0]
df_raw_ACC_filter_Plate_Center <- df_raw_ACC_filter_Plate_Center[df_meta_ACC$id, colSums(df_raw_ACC_filter_Plate_Center)>0]
df_raw_ACC_filter_putative <- df_raw_ACC_filter_putative[df_meta_ACC$id, colSums(df_raw_ACC_filter_putative)>0]
df_raw_ACC_filter_stringent <- df_raw_ACC_filter_stringent[df_meta_ACC$id, colSums(df_raw_ACC_filter_stringent)>0]

# Virus
df_raw_ACC_Vir <- df_raw_ACC %>% select(grep("Virus", names(df_raw_ACC)))
df_raw_ACC_filter_likely_Vir <- df_raw_ACC_filter_likely %>% select(grep("Virus", names(df_raw_ACC_filter_likely)))
df_raw_ACC_filter_Plate_Center_Vir <- df_raw_ACC_filter_Plate_Center %>% select(grep("Virus", names(df_raw_ACC_filter_Plate_Center)))
df_raw_ACC_filter_putative_Vir <- df_raw_ACC_filter_putative %>% select(grep("Virus", names(df_raw_ACC_filter_putative)))
df_raw_ACC_filter_stringent_Vir <- df_raw_ACC_filter_stringent %>% select(grep("Virus", names(df_raw_ACC_filter_stringent)))

# Archaea
df_raw_ACC_Arc <- df_raw_ACC %>% select(grep("Archaea", names(df_raw_ACC)))
df_raw_ACC_filter_likely_Arc <- df_raw_ACC_filter_likely %>% select(grep("Archaea", names(df_raw_ACC_filter_likely)))
df_raw_ACC_filter_Plate_Center_Arc <- df_raw_ACC_filter_Plate_Center %>% select(grep("Archaea", names(df_raw_ACC_filter_Plate_Center)))
df_raw_ACC_filter_putative_Arc <- df_raw_ACC_filter_putative %>% select(grep("Archaea", names(df_raw_ACC_filter_putative)))
df_raw_ACC_filter_stringent_Arc <- df_raw_ACC_filter_stringent %>% select(grep("Archaea", names(df_raw_ACC_filter_stringent)))

# Bacteria
df_raw_ACC_Bac <- df_raw_ACC %>% select(grep("Bacteria", names(df_raw_ACC)))
df_raw_ACC_filter_likely_Bac <- df_raw_ACC_filter_likely %>% select(grep("Bacteria", names(df_raw_ACC_filter_likely)))
df_raw_ACC_filter_Plate_Center_Bac <- df_raw_ACC_filter_Plate_Center %>% select(grep("Bacteria", names(df_raw_ACC_filter_Plate_Center)))
df_raw_ACC_filter_putative_Bac <- df_raw_ACC_filter_putative %>% select(grep("Bacteria", names(df_raw_ACC_filter_putative)))
df_raw_ACC_filter_stringent_Bac <- df_raw_ACC_filter_stringent %>% select(grep("Bacteria", names(df_raw_ACC_filter_stringent)))

## convert raw counts to relative abundance
df_raw_ACC_Vir_rb <- df_raw_ACC_Vir / rowSums(df_raw_ACC_Vir)
df_raw_ACC_filter_likely_Vir_rb <- df_raw_ACC_filter_likely_Vir / rowSums(df_raw_ACC_filter_likely_Vir)
df_raw_ACC_filter_Plate_Center_Vir_rb <- df_raw_ACC_filter_Plate_Center_Vir / rowSums(df_raw_ACC_filter_Plate_Center_Vir)
df_raw_ACC_filter_putative_Vir_rb <- df_raw_ACC_filter_putative_Vir / rowSums(df_raw_ACC_filter_putative_Vir)
df_raw_ACC_filter_stringent_Vir_rb <- df_raw_ACC_filter_stringent_Vir / rowSums(df_raw_ACC_filter_stringent_Vir)

df_raw_ACC_Arc_rb <- df_raw_ACC_Arc / rowSums(df_raw_ACC_Arc)
df_raw_ACC_filter_likely_Arc_rb <- df_raw_ACC_filter_likely_Arc / rowSums(df_raw_ACC_filter_likely_Arc)
df_raw_ACC_filter_Plate_Center_Arc_rb <- df_raw_ACC_filter_Plate_Center_Arc / rowSums(df_raw_ACC_filter_Plate_Center_Arc)
df_raw_ACC_filter_putative_Arc_rb <- df_raw_ACC_filter_putative_Arc / rowSums(df_raw_ACC_filter_putative_Arc)
df_raw_ACC_filter_stringent_Arc_rb <- df_raw_ACC_filter_stringent_Arc / rowSums(df_raw_ACC_filter_stringent_Arc)

df_raw_ACC_Bac_rb <- df_raw_ACC_Bac / rowSums(df_raw_ACC_Bac)
df_raw_ACC_filter_likely_Bac_rb <- df_raw_ACC_filter_likely_Bac / rowSums(df_raw_ACC_filter_likely_Bac)
df_raw_ACC_filter_Plate_Center_Bac_rb <- df_raw_ACC_filter_Plate_Center_Bac / rowSums(df_raw_ACC_filter_Plate_Center_Bac)
df_raw_ACC_filter_putative_Bac_rb <- df_raw_ACC_filter_putative_Bac / rowSums(df_raw_ACC_filter_putative_Bac)
df_raw_ACC_filter_stringent_Bac_rb <- df_raw_ACC_filter_stringent_Bac / rowSums(df_raw_ACC_filter_stringent_Bac)

# Voom-SNM normalized data
df_Voom_SNM_ACC_PCPG <- read.csv("Tables/Voom-SNM-ACC_PCPG.csv")
df_Voom_SNM_f_likely_ACC_PCPG <- read.csv("Tables/Voom-SNM-Filter-Likely-ACC_PCPG.csv")
df_Voom_SNM_f_PC_ACC_PCPG <- read.csv("Tables/Voom-SNM-Filter-Plate_Center-ACC_PCPG.csv")
df_Voom_SNM_f_putative_ACC_PCPG <- read.csv("Tables/Voom-SNM-Filter-Putative-ACC_PCPG.csv")
df_Voom_SNM_f_strigent_ACC_PCPG <- read.csv("Tables/Voom-SNM-Filter-Stringent-ACC_PCPG.csv")

df_Voom_SNM_ACC <- df_Voom_SNM_ACC_PCPG %>% filter(X %in% df_meta_ACC$id)
df_Voom_SNM_f_likely_ACC <- df_Voom_SNM_f_likely_ACC_PCPG %>% filter(X %in% df_meta_ACC$id)
df_Voom_SNM_f_PC_ACC <- df_Voom_SNM_f_PC_ACC_PCPG %>% filter(X %in% df_meta_ACC$id)
df_Voom_SNM_f_putative_ACC <- df_Voom_SNM_f_putative_ACC_PCPG %>% filter(X %in% df_meta_ACC$id)
df_Voom_SNM_f_strigent_ACC <- df_Voom_SNM_f_strigent_ACC_PCPG %>% filter(X %in% df_meta_ACC$id)

df_Voom_SNM_ACC <- df_Voom_SNM_ACC %>% remove_rownames %>% column_to_rownames(var="X")
df_Voom_SNM_f_likely_ACC <- df_Voom_SNM_f_likely_ACC %>% remove_rownames %>% column_to_rownames(var="X")
df_Voom_SNM_f_PC_ACC <- df_Voom_SNM_f_PC_ACC %>% remove_rownames %>% column_to_rownames(var="X")
df_Voom_SNM_f_putative_ACC <- df_Voom_SNM_f_putative_ACC %>% remove_rownames %>% column_to_rownames(var="X")
df_Voom_SNM_f_strigent_ACC <- df_Voom_SNM_f_strigent_ACC %>% column_to_rownames(var="X")

df_Voom_SNM_ACC_Vir <- df_Voom_SNM_ACC  %>% select(grep("Virus", names(df_Voom_SNM_ACC)))
df_Voom_SNM_f_likely_ACC_Vir <- df_Voom_SNM_f_likely_ACC  %>% select(grep("Virus", names(df_Voom_SNM_f_likely_ACC)))
df_Voom_SNM_f_PC_ACC_Vir <- df_Voom_SNM_f_PC_ACC  %>% select(grep("Virus", names(df_Voom_SNM_f_PC_ACC)))
df_Voom_SNM_f_putative_ACC_Vir <- df_Voom_SNM_f_putative_ACC  %>% select(grep("Virus", names(df_Voom_SNM_f_putative_ACC)))
df_Voom_SNM_f_strigent_ACC_Vir <- df_Voom_SNM_f_strigent_ACC  %>% select(grep("Virus", names(df_Voom_SNM_f_strigent_ACC)))

df_Voom_SNM_ACC_Arc <- df_Voom_SNM_ACC  %>% select(grep("Archaea", names(df_Voom_SNM_ACC)))
df_Voom_SNM_f_likely_ACC_Arc <- df_Voom_SNM_f_likely_ACC  %>% select(grep("Archaea", names(df_Voom_SNM_f_likely_ACC)))
df_Voom_SNM_f_PC_ACC_Arc <- df_Voom_SNM_f_PC_ACC  %>% select(grep("Archaea", names(df_Voom_SNM_f_PC_ACC)))
df_Voom_SNM_f_putative_ACC_Arc <- df_Voom_SNM_f_putative_ACC  %>% select(grep("Archaea", names(df_Voom_SNM_f_putative_ACC)))
df_Voom_SNM_f_strigent_ACC_Arc <- df_Voom_SNM_f_strigent_ACC  %>% select(grep("Archaea", names(df_Voom_SNM_f_strigent_ACC)))

df_Voom_SNM_ACC_Bac <- df_Voom_SNM_ACC  %>% select(grep("Bacteria", names(df_Voom_SNM_ACC)))
df_Voom_SNM_f_likely_ACC_Bac <- df_Voom_SNM_f_likely_ACC  %>% select(grep("Bacteria", names(df_Voom_SNM_f_likely_ACC)))
df_Voom_SNM_f_PC_ACC_Bac <- df_Voom_SNM_f_PC_ACC  %>% select(grep("Bacteria", names(df_Voom_SNM_f_PC_ACC)))
df_Voom_SNM_f_putative_ACC_Bac <- df_Voom_SNM_f_putative_ACC  %>% select(grep("Bacteria", names(df_Voom_SNM_f_putative_ACC)))
df_Voom_SNM_f_strigent_ACC_Bac <- df_Voom_SNM_f_strigent_ACC  %>% select(grep("Bacteria", names(df_Voom_SNM_f_strigent_ACC)))




# #################################################################################
# # examine sequencing depth
# #################################################################################
function_rearecurve <- function(in_matrix, title){
  rarecurve(in_matrix, step=1000, lwd=2,
            ylab="Number of genera",
            xlab="Sequencing reads per sample",
            main=title,
            label=F, se = TRUE)
  abline(v=min(rowSums(in_matrix)), col="red", lty=2)
}

pdf(file = "Figures/Rarecurve_ACC.pdf", width = 5, height = 5)
function_rearecurve(df_raw_ACC_Vir, "Raw_ACC_Virus")
function_rearecurve(df_raw_ACC_filter_likely_Vir, "Raw_ACC_Virus_filter_likely")
function_rearecurve(df_raw_ACC_filter_Plate_Center_Vir, "Raw_ACC_Virus_filter_Plate_Center")
function_rearecurve(df_raw_ACC_filter_putative_Vir, "Raw_ACC_Virus_filter_putative")
function_rearecurve(df_raw_ACC_filter_stringent_Vir, "Raw_ACC_Virus_filter_stringent")
function_rearecurve(df_raw_ACC_Arc, "Raw_ACC_Archaea")
function_rearecurve(df_raw_ACC_filter_likely_Arc, "Raw_ACC_Archaea_filter_likely")
function_rearecurve(df_raw_ACC_filter_Plate_Center_Arc, "Raw_ACC_Archaea_filter_Plate_Center")
function_rearecurve(df_raw_ACC_filter_putative_Arc, "Raw_ACC_Archaea_filter_putative")
function_rearecurve(df_raw_ACC_filter_stringent_Arc, "Raw_ACC_Archaea_filter_stringent")
function_rearecurve(df_raw_ACC_Bac, "Raw_ACC_Bacteria")
function_rearecurve(df_raw_ACC_filter_likely_Bac, "Raw_ACC_Bacteria_filter_likely")
function_rearecurve(df_raw_ACC_filter_Plate_Center_Bac, "Raw_ACC_Bacteria_filter_Plate_Center")
function_rearecurve(df_raw_ACC_filter_putative_Bac, "Raw_ACC_Bacteria_filter_putative")
function_rearecurve(df_raw_ACC_filter_stringent_Bac, "Raw_ACC_Bacteria_filter_stringent")
dev.off()



#################################################################################
# examine alpha-diversity
#################################################################################
ps_raw <- phyloseq(otu_table(as.matrix(df_raw_ACC), taxa_are_rows=FALSE))
ps_raw_filter_likely <- phyloseq(otu_table(as.matrix(df_raw_ACC_filter_likely), taxa_are_rows=FALSE))
ps_raw_filter_Plate_Center <- phyloseq(otu_table(as.matrix(df_raw_ACC_filter_Plate_Center), taxa_are_rows=FALSE))
ps_raw_filter_putative <- phyloseq(otu_table(as.matrix(df_raw_ACC_filter_putative), taxa_are_rows=FALSE))
ps_raw_filter_stringent <- phyloseq(otu_table(as.matrix(df_raw_ACC_filter_stringent), taxa_are_rows=FALSE))
ps_raw_Vir <- phyloseq(otu_table(as.matrix(df_raw_ACC_Vir), taxa_are_rows=FALSE))
ps_raw_filter_likely_Vir <- phyloseq(otu_table(as.matrix(df_raw_ACC_filter_likely_Vir), taxa_are_rows=FALSE))
ps_raw_filter_Plate_Center_Vir <- phyloseq(otu_table(as.matrix(df_raw_ACC_filter_Plate_Center_Vir), taxa_are_rows=FALSE))
ps_raw_filter_putative_Vir <- phyloseq(otu_table(as.matrix(df_raw_ACC_filter_putative_Vir), taxa_are_rows=FALSE))
ps_raw_filter_stringent_Vir <- phyloseq(otu_table(as.matrix(df_raw_ACC_filter_stringent_Vir), taxa_are_rows=FALSE))
ps_raw_Arc <- phyloseq(otu_table(as.matrix(df_raw_ACC_Arc), taxa_are_rows=FALSE))
ps_raw_filter_likely_Arc <- phyloseq(otu_table(as.matrix(df_raw_ACC_filter_likely_Arc), taxa_are_rows=FALSE))
ps_raw_filter_Plate_Center_Arc <- phyloseq(otu_table(as.matrix(df_raw_ACC_filter_Plate_Center_Arc), taxa_are_rows=FALSE))
ps_raw_filter_putative_Arc <- phyloseq(otu_table(as.matrix(df_raw_ACC_filter_putative_Arc), taxa_are_rows=FALSE))
ps_raw_filter_stringent_Arc <- phyloseq(otu_table(as.matrix(df_raw_ACC_filter_stringent_Arc), taxa_are_rows=FALSE))
ps_raw_Bac <- phyloseq(otu_table(as.matrix(df_raw_ACC_Bac), taxa_are_rows=FALSE))
ps_raw_filter_likely_Bac <- phyloseq(otu_table(as.matrix(df_raw_ACC_filter_likely_Bac), taxa_are_rows=FALSE))
ps_raw_filter_Plate_Center_Bac <- phyloseq(otu_table(as.matrix(df_raw_ACC_filter_Plate_Center_Bac), taxa_are_rows=FALSE))
ps_raw_filter_putative_Bac <- phyloseq(otu_table(as.matrix(df_raw_ACC_filter_putative_Bac), taxa_are_rows=FALSE))
ps_raw_filter_stringent_Bac <- phyloseq(otu_table(as.matrix(df_raw_ACC_filter_stringent_Bac), taxa_are_rows=FALSE))
# rarefy
set.seed(123)
ps_raw_raref <- rarefy_even_depth(ps_raw, sample.size=min(sample_sums(ps_raw)))
ps_raw_filter_likely_raref <- rarefy_even_depth(ps_raw_filter_likely, sample.size=min(sample_sums(ps_raw_filter_likely)))
ps_raw_filter_Plate_Center_raref <- rarefy_even_depth(ps_raw_filter_Plate_Center, sample.size=min(sample_sums(ps_raw_filter_Plate_Center)))
ps_raw_filter_putative_raref <- rarefy_even_depth(ps_raw_filter_putative, sample.size=min(sample_sums(ps_raw_filter_putative)))
ps_raw_filter_stringent_raref <- rarefy_even_depth(ps_raw_filter_stringent, sample.size=min(sample_sums(ps_raw_filter_stringent)))
ps_raw_Vir_raref <- rarefy_even_depth(ps_raw_Vir, sample.size=min(sample_sums(ps_raw_Vir)))
ps_raw_filter_likely_Vir_raref <- rarefy_even_depth(ps_raw_filter_likely_Vir, sample.size=min(sample_sums(ps_raw_filter_likely_Vir)))
ps_raw_filter_Plate_Center_Vir_raref <- rarefy_even_depth(ps_raw_filter_Plate_Center_Vir, sample.size=min(sample_sums(ps_raw_filter_Plate_Center_Vir)))
ps_raw_filter_putative_Vir_raref <- rarefy_even_depth(ps_raw_filter_putative_Vir, sample.size=min(sample_sums(ps_raw_filter_putative_Vir)))
# ps_raw_filter_stringent_Vir_raref <- rarefy_even_depth(ps_raw_filter_stringent_Vir, sample.size=min(sample_sums(ps_raw_filter_stringent_Vir)))
ps_raw_Arc_raref <- rarefy_even_depth(ps_raw_Arc, sample.size=min(sample_sums(ps_raw_Arc)))
ps_raw_filter_likely_Arc_raref <- rarefy_even_depth(ps_raw_filter_likely_Arc, sample.size=min(sample_sums(ps_raw_filter_likely_Arc)))
ps_raw_filter_Plate_Center_Arc_raref <- rarefy_even_depth(ps_raw_filter_Plate_Center_Arc, sample.size=min(sample_sums(ps_raw_filter_Plate_Center_Arc)))
ps_raw_filter_putative_Arc_raref <- rarefy_even_depth(ps_raw_filter_putative_Arc, sample.size=min(sample_sums(ps_raw_filter_putative_Arc)))
ps_raw_filter_stringent_Arc_raref <- rarefy_even_depth(ps_raw_filter_stringent_Arc, sample.size=min(sample_sums(ps_raw_filter_stringent_Arc)))
ps_raw_Bac_raref <- rarefy_even_depth(ps_raw_Bac, sample.size=min(sample_sums(ps_raw_Bac)))
ps_raw_filter_likely_Bac_raref <- rarefy_even_depth(ps_raw_filter_likely_Bac, sample.size=min(sample_sums(ps_raw_filter_likely_Bac)))
ps_raw_filter_Plate_Center_Bac_raref <- rarefy_even_depth(ps_raw_filter_Plate_Center_Bac, sample.size=min(sample_sums(ps_raw_filter_Plate_Center_Bac)))
ps_raw_filter_putative_Bac_raref <- rarefy_even_depth(ps_raw_filter_putative_Bac, sample.size=min(sample_sums(ps_raw_filter_putative_Bac)))
ps_raw_filter_stringent_Bac_raref <- rarefy_even_depth(ps_raw_filter_stringent_Bac, sample.size=min(sample_sums(ps_raw_filter_stringent_Bac)))

# alpha diversity
alpha_raw <- estimate_richness(ps_raw, split = TRUE)
alpha_raw_filter_likely <- estimate_richness(ps_raw_filter_likely, split = TRUE)
alpha_raw_filter_Plate_Center <- estimate_richness(ps_raw_filter_Plate_Center, split = TRUE)
alpha_raw_filter_putative <- estimate_richness(ps_raw_filter_putative, split = TRUE)
alpha_raw_filter_stringent <- estimate_richness(ps_raw_filter_stringent, split = TRUE)
alpha_raw_Vir <- estimate_richness(ps_raw_Vir, split = TRUE)
alpha_raw_filter_likely_Vir <- estimate_richness(ps_raw_filter_likely_Vir, split = TRUE)
alpha_raw_filter_Plate_Center_Vir <- estimate_richness(ps_raw_filter_Plate_Center_Vir, split = TRUE)
alpha_raw_filter_putative_Vir <- estimate_richness(ps_raw_filter_putative_Vir, split = TRUE)
# alpha_raw_filter_stringent_Vir <- estimate_richness(ps_raw_filter_stringent_Vir, split = TRUE)
alpha_raw_Arc <- estimate_richness(ps_raw_Arc, split = TRUE)
alpha_raw_filter_likely_Arc <- estimate_richness(ps_raw_filter_likely_Arc, split = TRUE)
alpha_raw_filter_Plate_Center_Arc <- estimate_richness(ps_raw_filter_Plate_Center_Arc, split = TRUE)
alpha_raw_filter_putative_Arc <- estimate_richness(ps_raw_filter_putative_Arc, split = TRUE)
# alpha_raw_filter_stringent_Arc <- estimate_richness(ps_raw_filter_stringent_Arc, split = TRUE)
alpha_raw_Bac <- estimate_richness(ps_raw_Bac, split = TRUE)
alpha_raw_filter_likely_Bac <- estimate_richness(ps_raw_filter_likely_Bac, split = TRUE)
alpha_raw_filter_Plate_Center_Bac <- estimate_richness(ps_raw_filter_Plate_Center_Bac, split = TRUE)
alpha_raw_filter_putative_Bac <- estimate_richness(ps_raw_filter_putative_Bac, split = TRUE)
alpha_raw_filter_stringent_Bac <- estimate_richness(ps_raw_filter_stringent_Bac, split = TRUE)

alpha_raw_raref <- estimate_richness(ps_raw_raref, split = TRUE)
alpha_raw_filter_likely_raref <- estimate_richness(ps_raw_filter_likely_raref, split = TRUE)
alpha_raw_filter_Plate_Center_raref <- estimate_richness(ps_raw_filter_Plate_Center_raref, split = TRUE)
alpha_raw_filter_putative_raref <- estimate_richness(ps_raw_filter_putative_raref, split = TRUE)
alpha_raw_filter_stringent_raref <- estimate_richness(ps_raw_filter_stringent_raref, split = TRUE)

alpha_raw_Vir_raref <- estimate_richness(ps_raw_Vir_raref, split = TRUE)
alpha_raw_filter_likely_Vir_raref <- estimate_richness(ps_raw_filter_likely_Vir_raref, split = TRUE)
alpha_raw_filter_Plate_Center_Vir_raref <- estimate_richness(ps_raw_filter_Plate_Center_Vir_raref, split = TRUE)
alpha_raw_filter_putative_Vir_raref <- estimate_richness(ps_raw_filter_putative_Vir_raref, split = TRUE)
# alpha_raw_filter_stringent_Vir_raref <- estimate_richness(ps_raw_filter_stringent_Vir_raref, split = TRUE)
alpha_raw_Arc_raref <- estimate_richness(ps_raw_Arc_raref, split = TRUE)
alpha_raw_filter_likely_Arc_raref <- estimate_richness(ps_raw_filter_likely_Arc_raref, split = TRUE)
alpha_raw_filter_Plate_Center_Arc_raref <- estimate_richness(ps_raw_filter_Plate_Center_Arc_raref, split = TRUE)
alpha_raw_filter_putative_Arc_raref <- estimate_richness(ps_raw_filter_putative_Arc_raref, split = TRUE)
# alpha_raw_filter_stringent_Arc_raref <- estimate_richness(ps_raw_filter_stringent_Arc_raref, split = TRUE)
alpha_raw_Bac_raref <- estimate_richness(ps_raw_Bac_raref, split = TRUE)
alpha_raw_filter_likely_Bac_raref <- estimate_richness(ps_raw_filter_likely_Bac_raref, split = TRUE)
alpha_raw_filter_Plate_Center_Bac_raref <- estimate_richness(ps_raw_filter_Plate_Center_Bac_raref, split = TRUE)
alpha_raw_filter_putative_Bac_raref <- estimate_richness(ps_raw_filter_putative_Bac_raref, split = TRUE)
alpha_raw_filter_stringent_Bac_raref <- estimate_richness(ps_raw_filter_stringent_Bac_raref, split = TRUE)

# compare the diversity with and without rarefaction
cor.test(alpha_raw$Shannon, alpha_raw_raref$Shannon, method="spearman")
cor.test(alpha_raw_filter_likely$Shannon, alpha_raw_filter_likely_raref$Shannon, method="spearman")
cor.test(alpha_raw_filter_Plate_Center$Shannon, alpha_raw_filter_Plate_Center_raref$Shannon, method="spearman")
cor.test(alpha_raw_filter_putative$Shannon, alpha_raw_filter_putative_raref$Shannon, method="spearman")
cor.test(alpha_raw_filter_stringent$Shannon, alpha_raw_filter_stringent_raref$Shannon, method="spearman")
cor.test(alpha_raw_Vir$Shannon, alpha_raw_Vir_raref$Shannon, method="spearman")
cor.test(alpha_raw_filter_likely_Vir$Shannon, alpha_raw_filter_likely_Vir_raref$Shannon, method="spearman")
cor.test(alpha_raw_filter_Plate_Center_Vir$Shannon, alpha_raw_filter_Plate_Center_Vir_raref$Shannon, method="spearman")
cor.test(alpha_raw_filter_putative_Vir$Shannon, alpha_raw_filter_putative_Vir_raref$Shannon, method="spearman")
# cor.test(alpha_raw_filter_stringent_Vir$Shannon, alpha_raw_filter_stringent_Vir_raref$Shannon, method="spearman")
cor.test(alpha_raw_Arc$Shannon, alpha_raw_Arc_raref$Shannon, method="spearman")
cor.test(alpha_raw_filter_likely_Arc$Shannon, alpha_raw_filter_likely_Arc_raref$Shannon, method="spearman")
cor.test(alpha_raw_filter_Plate_Center_Arc$Shannon, alpha_raw_filter_Plate_Center_Arc_raref$Shannon, method="spearman")
cor.test(alpha_raw_filter_putative_Arc$Shannon, alpha_raw_filter_putative_Arc_raref$Shannon, method="spearman")
# cor.test(alpha_raw_filter_stringent_Arc$Shannon, alpha_raw_filter_stringent_Arc_raref$Shannon, method="spearman")
cor.test(alpha_raw_Bac$Shannon, alpha_raw_Bac_raref$Shannon, method="spearman")
cor.test(alpha_raw_filter_likely_Bac$Shannon, alpha_raw_filter_likely_Bac_raref$Shannon, method="spearman")
cor.test(alpha_raw_filter_Plate_Center_Bac$Shannon, alpha_raw_filter_Plate_Center_Bac_raref$Shannon, method="spearman")
cor.test(alpha_raw_filter_putative_Bac$Shannon, alpha_raw_filter_putative_Bac_raref$Shannon, method="spearman")
cor.test(alpha_raw_filter_stringent_Bac$Shannon, alpha_raw_filter_stringent_Bac_raref$Shannon, method="spearman")

df_alpha <- data.frame(id=rownames(alpha_raw_Vir),
                       alpha=alpha_raw$Shannon, alpha_rf=alpha_raw_raref$Shannon,
                       alpha_f_likely=alpha_raw_filter_likely$Shannon, alpha_f_likely_rf=alpha_raw_filter_likely_raref$Shannon,
                       alpha_f_PC=alpha_raw_filter_Plate_Center$Shannon, alpha_f_PC_rf=alpha_raw_filter_Plate_Center_raref$Shannon,
                       alpha_f_Putative=alpha_raw_filter_putative$Shannon, alpha_f_Putative_rf=alpha_raw_filter_Plate_Center_raref$Shannon,
                       alpha_f_Strigent=alpha_raw_filter_stringent$Shannon, alpha_f_Strigent_rf=alpha_raw_filter_stringent_raref$Shannon,
                       alpha_Vir=alpha_raw_Vir$Shannon, alpha_Vir_rf=alpha_raw_Vir_raref$Shannon,
                       alpha_Vir_f_Likely=alpha_raw_filter_likely_Vir$Shannon, alpha_Vir_f_Likely_rf=alpha_raw_filter_likely_Vir_raref$Shannon,
                       alpha_Vir_f_PC=alpha_raw_filter_Plate_Center_Vir$Shannon, alpha_Vir_f_PC_rf=alpha_raw_filter_Plate_Center_Vir_raref$Shannon,
                       alpha_Vir_f_Putative=alpha_raw_filter_putative_Vir$Shannon, alpha_Vir_f_Putative_rf=alpha_raw_filter_putative_Vir_raref$Shannon,
                       alpha_Arc=alpha_raw_Arc$Shannon, alpha_Arc_rf=alpha_raw_Arc_raref$Shannon,
                       alpha_Arc_f_Likely=alpha_raw_filter_likely_Arc$Shannon, alpha_Arc_f_Likely_rf=alpha_raw_filter_likely_Arc_raref$Shannon,
                       alpha_Arc_f_PC=alpha_raw_filter_Plate_Center_Arc$Shannon, alpha_Arc_f_PC_rf=alpha_raw_filter_Plate_Center_Arc_raref$Shannon,
                       alpha_Arc_f_Putative=alpha_raw_filter_putative_Arc$Shannon, alpha_Arc_f_Putative_rf=alpha_raw_filter_putative_Arc_raref$Shannon,
                       alpha_Bac=alpha_raw_Bac$Shannon, alpha_Bac_rf=alpha_raw_Bac_raref$Shannon,
                       alpha_Bac_f_Likely=alpha_raw_filter_likely_Bac$Shannon, alpha_Bac_f_Likely_rf=alpha_raw_filter_likely_Bac_raref$Shannon,
                       alpha_Bac_f_PC=alpha_raw_filter_Plate_Center_Bac$Shannon, alpha_Bac_f_PC_rf=alpha_raw_filter_Plate_Center_Bac_raref$Shannon,
                       alpha_Bac_f_Putative=alpha_raw_filter_putative_Bac$Shannon, alpha_Bac_f_Putative_rf=alpha_raw_filter_putative_Bac_raref$Shannon,
                       alpha_Bac_f_Strigent=alpha_raw_filter_stringent_Bac$Shannon, alpha_Bac_f_Strigent_rf=alpha_raw_filter_stringent_Bac_raref$Shannon)
write.csv(df_alpha, file = "Tables/Alpha_diversity_ACC.csv", row.names = FALSE)

df_meta_alpha <- merge(df_meta_ACC, df_alpha, by="id", all.x = TRUE)
write.csv(df_meta_alpha, file = "Tables/Alpha_diversity_ACC_metadata.csv", row.names = FALSE)



#################################################################################
# divide patients into two groups according to overall survival
#################################################################################
# Long-term survivors and short-term survivors
median_ST <- median(df_meta_ACC$OS_MONTHS)
df_meta_ACC$Group_ST <- ifelse(df_meta_ACC$OS_MONTHS>median_ST, "LTS", "STS")
# check the microbial difference between two groups
function_PcoA <- function(in_matrix, title){
  in_matrix[in_matrix<0] <- 0 # assign zero to negative values
  distance <- vegdist(in_matrix, method = 'bray')
  pcoa <- cmdscale(distance, k = (nrow(in_matrix) - 1), eig = TRUE)
  point <- data.frame(pcoa$point)
  pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
  # extract first two coordinate value
  sample_eig <- data.frame({pcoa$point})[1:2]
  sample_eig$id <- rownames(sample_eig)
  names(sample_eig)[1:2] <- c('PCoA1', 'PCoA2')
  group_m <- merge(sample_eig, df_meta_ACC, by = "id", all.x = TRUE)
  # PERMANOVA
  dt <- in_matrix[group_m$id,]
  identical(rownames(dt), group_m$id)
  adonis_result <- adonis2(dt~Group_ST, group_m, permutations = 999, distance = 'bray')
  R2 <- round(adonis_result$R2[1], 2)
  pvalue <- round(adonis_result$`Pr(>F)`[1], 2)
  plot <- ggscatter(group_m, x= "PCoA1", y = "PCoA2",color="Group_ST",
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

pdf(file = "Figures/PcoA_ACC_LTS_STS_Relative_Abundance.pdf", width = 6, height = 5)
function_PcoA(df_raw_ACC_Vir_rb, "Raw_virus")
function_PcoA(df_raw_ACC_filter_likely_Vir_rb, "Raw_filter_likely_virus")
function_PcoA(df_raw_ACC_filter_Plate_Center_Vir_rb, "Raw_filter_Plate_Center_virus")
function_PcoA(df_raw_ACC_filter_putative_Vir_rb, "Raw_filter_putative_virus")
# function_PcoA(df_raw_ACC_filter_stringent_Vir_rb, "Raw_filter_stringent_virus")
function_PcoA(df_raw_ACC_Arc_rb, "Raw_archaea")
function_PcoA(df_raw_ACC_filter_likely_Arc_rb, "Raw_filter_likely_archaea")
function_PcoA(df_raw_ACC_filter_Plate_Center_Arc_rb, "Raw_filter_Plate_Center_archaea")
function_PcoA(df_raw_ACC_filter_putative_Arc_rb, "Raw_filter_putative_archaea")
# function_PcoA(df_raw_ACC_filter_stringent_Arc_rb, "Raw_filter_stringent_archaea")
function_PcoA(df_raw_ACC_Bac_rb, "Raw_bacteria")
function_PcoA(df_raw_ACC_filter_likely_Bac_rb, "Raw_filter_likely_bacteria")
function_PcoA(df_raw_ACC_filter_Plate_Center_Bac_rb, "Raw_filter_Plate_Center_bacteria")
function_PcoA(df_raw_ACC_filter_putative_Bac_rb, "Raw_filter_putative_bacteria")
function_PcoA(df_raw_ACC_filter_stringent_Bac_rb, "Raw_filter_stringent_bacteria")
dev.off()

pdf(file = "Figures/PcoA_ACC_LTS_STS_Voom-SNM.pdf", width = 6, height = 5)
function_PcoA(df_Voom_SNM_ACC, "Voom_SNM_Overall")
function_PcoA(df_Voom_SNM_f_likely_ACC, "Voom_SNM_filter_likely_Overall")
function_PcoA(df_Voom_SNM_f_PC_ACC, "Voom_SNM_filter_Plate_Center_Overall")
function_PcoA(df_Voom_SNM_f_putative_ACC, "Voom_SNM_filter_putative_Overall")
function_PcoA(df_Voom_SNM_f_strigent_ACC, "Voom_SNM_filter_stringent_Overall")
function_PcoA(df_Voom_SNM_ACC_Vir, "Voom_SNM_virus")
function_PcoA(df_Voom_SNM_f_likely_ACC_Vir, "Voom_SNM_filter_likely_virus")
function_PcoA(df_Voom_SNM_f_PC_ACC_Vir, "Voom_SNM_filter_Plate_Center_virus")
function_PcoA(df_Voom_SNM_f_putative_ACC_Vir, "Voom_SNM_filter_putative_virus")
function_PcoA(df_Voom_SNM_ACC_Arc, "Voom_SNM_archaea")
function_PcoA(df_Voom_SNM_f_likely_ACC_Arc, "Voom_SNM_filter_likely_archaea")
function_PcoA(df_Voom_SNM_f_PC_ACC_Arc, "Voom_SNM_filter_Plate_Center_archaea")
function_PcoA(df_Voom_SNM_f_putative_ACC_Arc, "Voom_SNM_filter_putative_archaea")
function_PcoA(df_Voom_SNM_ACC_Bac, "Voom_SNM_bacteria")
function_PcoA(df_Voom_SNM_f_likely_ACC_Bac, "Voom_SNM_filter_likely_bacteria")
function_PcoA(df_Voom_SNM_f_PC_ACC_Bac, "Voom_SNM_filter_Plate_Center_bacteria")
function_PcoA(df_Voom_SNM_f_putative_ACC_Bac, "Voom_SNM_filter_putative_bacteria")
function_PcoA(df_Voom_SNM_f_strigent_ACC_Bac, "Voom_SNM_filter_stringent_bacteria")
dev.off()


#################################################################################
# divide patients into two groups according to shannon diversity
#################################################################################
df_meta_alpha$Group_bac <- ifelse(df_meta_alpha$alpha_Bac>median(df_meta_alpha$alpha_Bac), "high", "low")
df_meta_alpha$Group_bac_fL <- ifelse(df_meta_alpha$alpha_Bac_f_Likely>median(df_meta_alpha$alpha_Bac_f_Likely), "high", "low")
df_meta_alpha$Group_bac_fPC <- ifelse(df_meta_alpha$alpha_Bac_f_PC>median(df_meta_alpha$alpha_Bac_f_PC), "high", "low")
df_meta_alpha$Group_bac_fP <- ifelse(df_meta_alpha$alpha_Bac_f_Putative>median(df_meta_alpha$alpha_Bac_f_Putative), "high", "low")
df_meta_alpha$Group_bac_fS <- ifelse(df_meta_alpha$alpha_Bac_f_Strigent>median(df_meta_alpha$alpha_Bac_f_Strigent), "high", "low")

function_sur <- function(group, title){
  commands <- paste("fit <- survfit(Surv(OS_MONTHS, CENSOR) ~ ", group, ", data = df_meta_alpha)", sep = "")
  eval(parse(text = commands))
  
#fit <- survfit(Surv(OS_MONTHS, CENSOR) ~ group, data = df_meta_alpha)
  p <- ggsurvplot(fit, # 创建的拟合对象
                             data = df_meta_alpha,  # 指定变量数据来源
                             conf.int = TRUE, # 显示置信区间
                             pval = TRUE, # 添加P值
                             palette = c("#0072B5CC","#E18727CC"),
                             title = title,
                             xlab = "Time (Monthes)", 
                             legend = "right",
                             legend.title = "",
                             ggtheme = theme_base(),
                             break.x.by = 40)
  return(p)
}

pdf(file = "Figures/Survivor_BAC_alpha_diversity.pdf", width = 7, height = 5)
function_sur("Group_bac", "Raw_bacteria")
function_sur("Group_bac_fL", "Raw_bacteria_filter_likely")
function_sur("Group_bac_fPC", "Raw_bacteria_filter_Plate_Center")
function_sur("Group_bac_fP", "Raw_bacteria_filter_putative")
function_sur("Group_bac_fS", "Raw_bacteria_filter_strigent")
dev.off()
