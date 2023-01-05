# load libraries
library(ggplot2)
library(dplyr)
library(Maaslin2)
library(tidyverse)
library(UpSetR)
library(stringr)
library(circlize)
library(ComplexHeatmap)
library(survival)
library(survminer)
library(ggthemes)
library(ggsci)



df_meta_ACC <- read.csv("Tables/metadata_77_ACC.csv")
df_clinical <- read.csv("01.Data/ACC_clinic.csv")

# load microbial data from Poore's study
df_ACC_Voom_SNM <- read.csv("Tables/Voom-SNM-ACC.csv", row.names = 1)
df_ACC_filter_likely_Voom_SNM  <- read.csv("Tables/Voom-SNM-Filter-Likely-ACC.csv", row.names = 1)
df_ACC_filter_Plate_Center_Voom_SNM <- read.csv("Tables/Voom-SNM-Filter-Plate_Center-ACC.csv", row.names = 1)
df_ACC_filter_putative_Voom_SNM <- read.csv("Tables/Voom-SNM-Filter-Putative-ACC.csv", row.names = 1)
df_ACC_filter_stringent_Voom_SNM <- read.csv("Tables/Voom-SNM-Filter-Stringent-ACC.csv", row.names = 1)

# load clustering results
df_cluster <- read.csv("Tables/Unspervised_clustering.csv")
df_cluster <- df_cluster %>% dplyr::select(c("id", "Without_combined", "Likely_combined", "PC_combined", "Putative_combined"))
df_cluster[df_cluster==1] <- "MS1"
df_cluster[df_cluster==2] <- "MS2"
df_cluster_clinical <- merge(df_cluster, df_clinical, by="id", all.x = TRUE)

df_cluster_clinical <- df_cluster_clinical %>% remove_rownames() %>% column_to_rownames(var="id")
df_cluster_clinical$Without_combined <- factor(df_cluster_clinical$Without_combined, levels = c("MS1", "MS2"))
df_cluster_clinical$Likely_combined <- factor(df_cluster_clinical$Likely_combined, levels = c("MS1", "MS2"))
df_cluster_clinical$PC_combined <- factor(df_cluster_clinical$PC_combined, levels = c("MS1", "MS2"))
df_cluster_clinical$Putative_combined <- factor(df_cluster_clinical$Putative_combined, levels = c("MS1", "MS2"))



#######################################################################################################
# Differential genus between two groups using "Maaslin2"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
#######################################################################################################
function_Maaslin2 <- function(input_data, groups, out){
  input_data <- input_data[rownames(df_cluster_clinical), ]
  
  fit_data = Maaslin2(
    input_data = input_data, 
    input_metadata = df_cluster_clinical,
    min_prevalence = 0,
    normalization="NONE",
    transform="NONE",
    analysis_method = "LM",
    output = out, 
    fixed_effects = groups,
    cores = 6,
    random_effects = c("SEX", "Ethnicity_Race", "AGE", "Neoplasm_Status", "Clinical_Status_3_Mo_Post.Op", "Surgical_margin", "ATYPICAL_MITOTIC_FIGURES"))
  return(fit_data)
}

maaslin_without_group_without <- function_Maaslin2(df_ACC_Voom_SNM, "Without_combined", "Temp/without_group_without")
maaslin_without_group_likely <- function_Maaslin2(df_ACC_Voom_SNM, "Likely_combined", "Temp/without_group_likely")
maaslin_without_group_PC <- function_Maaslin2(df_ACC_Voom_SNM, "PC_combined", "Temp/without_group_PC")
maaslin_without_group_putative <- function_Maaslin2(df_ACC_Voom_SNM, "Putative_combined", "Temp/without_group_putative")

maaslin_likely_group_without <- function_Maaslin2(df_ACC_filter_likely_Voom_SNM, "Without_combined", "Temp/likely_group_without")
maaslin_likely_group_likely <- function_Maaslin2(df_ACC_filter_likely_Voom_SNM, "Likely_combined", "Temp/likely_group_likely")
maaslin_likely_group_PC <- function_Maaslin2(df_ACC_filter_likely_Voom_SNM, "PC_combined", "Temp/likely_group_PC")
maaslin_likely_group_putative <- function_Maaslin2(df_ACC_filter_likely_Voom_SNM, "Putative_combined", "Temp/likely_group_putative")

maaslin_PC_group_without <- function_Maaslin2(df_ACC_filter_Plate_Center_Voom_SNM, "Without_combined", "Temp/PC_group_without")
maaslin_PC_group_likely <- function_Maaslin2(df_ACC_filter_Plate_Center_Voom_SNM, "Likely_combined", "Temp/PC_group_likely")
maaslin_PC_group_PC <- function_Maaslin2(df_ACC_filter_Plate_Center_Voom_SNM, "PC_combined", "Temp/PC_group_PC")
maaslin_PC_group_putative <- function_Maaslin2(df_ACC_filter_Plate_Center_Voom_SNM, "Putative_combined", "Temp/PC_group_putative")

maaslin_putative_group_without <- function_Maaslin2(df_ACC_filter_putative_Voom_SNM, "Without_combined", "Temp/putative_group_without")
maaslin_putative_group_likely <- function_Maaslin2(df_ACC_filter_putative_Voom_SNM, "Likely_combined", "Temp/putative_group_likely")
maaslin_putative_group_PC <- function_Maaslin2(df_ACC_filter_putative_Voom_SNM, "PC_combined", "Temp/putative_group_PC")
maaslin_putative_group_putative <- function_Maaslin2(df_ACC_filter_putative_Voom_SNM, "Putative_combined", "Temp/putative_group_putative")

maaslin_strigent_group_without <- function_Maaslin2(df_ACC_filter_stringent_Voom_SNM, "Without_combined", "Temp/strigent_group_without")
maaslin_strigent_group_likely <- function_Maaslin2(df_ACC_filter_stringent_Voom_SNM, "Likely_combined", "Temp/strigent_group_likely")
maaslin_strigent_group_PC <- function_Maaslin2(df_ACC_filter_stringent_Voom_SNM, "PC_combined", "Temp/strigent_group_PC")
maaslin_strigent_group_putative <- function_Maaslin2(df_ACC_filter_stringent_Voom_SNM, "Putative_combined", "Temp/strigent_group_putative")

# extract significantly differential features
fuction_feature <- function(df1, df2, df3, df4){
  maaslin_without <- df1$results %>% filter(qval<0.05) %>% dplyr::select(c("feature", "coef", "qval")) %>%
    setNames(c("feature", "coef_without", "qval_without"))
  maaslin_likely<- df2$results %>% filter(qval<0.05) %>% dplyr::select(c("feature", "coef", "qval")) %>%
    setNames(c("feature", "coef_likely", "qval_likely"))
  maaslin_PC <- df3$results %>% filter(qval<0.05) %>% dplyr::select(c("feature", "coef", "qval")) %>%
    setNames(c("feature", "coef_PC", "qval_PC"))
  maaslin_putative <- df4$results %>% filter(qval<0.05) %>% dplyr::select(c("feature", "coef", "qval")) %>%
    setNames(c("feature", "coef_putative", "qval_putative"))
  
  df_combine <- Reduce(function(x,y){merge(x,y,by="feature",all=TRUE)},
                       list(maaslin_without, maaslin_likely, maaslin_PC, maaslin_putative))
  df_combine_feature <- df_combine %>% dplyr::select(c("coef_without", "coef_likely", "coef_PC", "coef_putative"))
  df_combine_feature[is.na(df_combine_feature)] <- 0
  df_combine_feature[df_combine_feature!=0] <- 1
  
  plot_upset <- upset(df_combine_feature, 
                      keep.order=TRUE, sets.x.label = "Number of differential microbes",
                      queries=list(list(query = intersects, params = list("coef_without","coef_likely","coef_PC","coef_putative"), color = "red",active = T)))
  return(list(df_combine, plot_upset))
}

list_without <- fuction_feature(maaslin_without_group_without, maaslin_without_group_likely, maaslin_without_group_PC, maaslin_without_group_putative)
list_likely <- fuction_feature(maaslin_likely_group_without, maaslin_likely_group_likely, maaslin_likely_group_PC, maaslin_likely_group_putative)
list_PC <- fuction_feature(maaslin_PC_group_without, maaslin_PC_group_likely, maaslin_PC_group_PC, maaslin_PC_group_putative)
list_putative <- fuction_feature(maaslin_putative_group_without, maaslin_putative_group_likely, maaslin_putative_group_PC, maaslin_putative_group_putative)
# list_strigent <- fuction_feature(maaslin_strigent_group_without, maaslin_strigent_group_likely, maaslin_strigent_group_PC, maaslin_strigent_group_putative)
# No siginificantly differential features detected in stringent group under four classifications

df_feature_without <- list_without[[1]]
df_feature_likely <- list_likely[[1]]
df_feature_PC <- list_PC[[1]]
df_feature_putative <- list_putative[[1]]

df_feature_without_core <- na.omit(df_feature_without)
df_feature_likely_core <- na.omit(df_feature_likely)
df_feature_PC_core <- na.omit(df_feature_PC)
df_feature_putative_core <- na.omit(df_feature_putative)

pdf(file = "Figures/Significantly_differential_microbes_Upset.pdf",height = 3,width = 6)
list_without[[2]]
list_likely[[2]]
list_PC[[2]]
list_putative[[2]]
dev.off()


# core features that are invariably detected under different classifications
core_features <- Reduce(intersect, list(df_feature_without_core$feature,df_feature_likely_core$feature,
                                        df_feature_PC_core$feature,df_feature_putative_core$feature))

# Visulazation

df_core_features <- Reduce(function(x,y){merge(x,y,by="feature",all=TRUE)},
                           list(df_feature_without_core, df_feature_likely_core, df_feature_PC_core, df_feature_putative_core))
df_core_features <- na.omit(df_core_features)

df_core_features$avg.coef <- rowMeans(df_core_features[ , grep("coef", names(df_core_features))])
df_core_features$avg.qval <- rowMeans(df_core_features[ , grep("qval", names(df_core_features))])
df_core_features$log.avg.qval <- -log(df_core_features$avg.qval)
df_core_features <- df_core_features[order(df_core_features$avg.coef), ]
df_core_features$class <- ifelse(df_core_features$avg.coef < -1, "Enrich_MS1", ifelse(df_core_features$avg.coef > 1, "Enrich_MS2", "Not"))
df_core_features$class <- factor(df_core_features$class, levels = c("Enrich_MS1", "Enrich_MS2", "Not"))

# p <- ggplot(data = df_core_features[,c("avg.coef", "log.avg.qval", "class")], aes(x=avg.coef, y=log.avg.qval, color=class))+
#         geom_point()+
#         scale_color_manual(values = c("#0072B5CC","#E18727CC", "#A0A0A0"))+
#         theme_bw()+
#         xlab("Average Coefficients")+
#         ylab("-Log adjusted p-value")+
#         geom_vline(xintercept = c(-1, 1), linetype="dotted", color = "red", size=1)

dt <- df_core_features[,c("feature", "avg.coef", "log.avg.qval")]
dt$kingdom <- str_extract(df_core_features$feature, "k__[a-zA-Z]+")
dt$kingdom <- factor(dt$kingdom, levels = c("k__Viruses", "k__Archaea", "k__Bacteria"))
p<- ggplot(data = dt, aes(x=avg.coef, y=log.avg.qval, fill=kingdom))+
        geom_point(aes(colour=kingdom)) +
        #scale_color_gradient2(low="#0072B5CC", high="#E18727CC", midpoint=0)+
        scale_color_aaas()+
        theme_bw()+
        xlab("Average Coefficients")+
        ylab("-Log adjusted p-value")+
        geom_vline(xintercept = c(-1, 1), linetype="dotted", color = "red", size=1)

ggsave(filename = "Figures/Core_features_coefficients.pdf", p, useDingbats = FALSE, width = 5,height = 3.8)

write_csv(df_core_features, file = "Tables/Core_features_coefficients.csv")



#######################################################################################################
# Stratify patients based on the abundance of significant microbes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#######################################################################################################
sig_microbes <- df_core_features[abs(df_core_features$avg.coef)>1, ]$feature

function_ex_sur <- function(in_matrix){
  df_cox <- data.frame()
  list_re <- list()
  for (i in 1:length(sig_microbes)){
    genus <- sig_microbes[i]
    if (genus %in% names(in_matrix)){
      df_s <- in_matrix[,sig_microbes[i], drop=FALSE]
      df_s$class <- ifelse(df_s[,1]>median(df_s[,1]), "high", "low")
      df_s$id <- rownames(df_s)
      
      df_meta_ACC_m <- merge(df_meta_ACC[, c("id", "OS_MONTHS", "CENSOR")], df_s, by="id", all.x = TRUE)
      
      # survival plot
      fit <- survfit(Surv(OS_MONTHS, CENSOR) ~ class, data = df_meta_ACC_m)
      
      p <- ggsurvplot(fit,
                      data = df_meta_ACC_m,
                      conf.int = TRUE,
                      pval = TRUE,
                      palette = c("#0072B5CC","#E18727CC"),
                      xlab = "Time (Monthes)", 
                      legend = "right",
                      legend.title = "",
                      ggtheme = theme_base(),
                      break.x.by = 40)+
        ggtitle(str_extract(genus, "g__.+"))
      # cox analysis
      df_meta_ACC_m$class <- factor(df_meta_ACC_m$class, levels = c("low", "high"))
      res.cox <- coxph(Surv(OS_MONTHS, CENSOR) ~ class, data = df_meta_ACC_m)
      cox <- summary(res.cox)
      df_hr <- cox$coefficients %>% as.data.frame()
      rownames(df_hr) <- genus
      df_cox <- rbind(df_cox, df_hr)
      
      list_re[[i]] <- p
    }
  }
  return(list(list_re, df_cox))
}

p_without <- function_ex_sur(df_ACC_Voom_SNM)
p_likely <- function_ex_sur(df_ACC_filter_likely_Voom_SNM)
p_PC <- function_ex_sur(df_ACC_filter_Plate_Center_Voom_SNM)
p_putative <- function_ex_sur(df_ACC_filter_putative_Voom_SNM)
p_stringent <- function_ex_sur(df_ACC_filter_stringent_Voom_SNM)

pdf(file = "Figures/Biomarker_survial_examination_Voom-SNM_Without.pdf",height = 5,width = 6)
p_without[[1]]
dev.off()
pdf(file = "Figures/Biomarker_survial_examination_Voom-SNM_Filter_Likely.pdf",height = 5,width = 6)
p_likely[[1]]
dev.off()
pdf(file = "Figures/Biomarker_survial_examination_Voom-SNM_Filter_PC.pdf",height = 5,width = 6)
p_PC[[1]]
dev.off()
pdf(file = "Figures/Biomarker_survial_examination_Voom-SNM_Filter_Putative.pdf",height = 5,width = 6)
p_putative[[1]]
dev.off()
pdf(file = "Figures/Biomarker_survial_examination_Voom-SNM_Filter_Stringent.pdf",height = 5,width = 6)
p_stringent[[1]]
dev.off()


# combine the cox results
df_cox_without <- p_without[[2]]
df_cox_likely <- p_likely[[2]]
df_cox_PC <- p_PC[[2]]
df_cox_putative <- p_putative[[2]]

df_cox_without$genus <- rownames(df_cox_without)
df_cox_likely$genus <- rownames(df_cox_likely)
df_cox_PC$genus <- rownames(df_cox_PC)
df_cox_putative$genus <- rownames(df_cox_putative)

df_cox_without_s <- df_cox_without %>% dplyr::select(c("genus" , "exp(coef)"), "Pr(>|z|)") %>% setNames(c("genus", "HR_without", "pvalue_without"))
df_cox_likely_s <- df_cox_likely %>% dplyr::select(c("genus" , "exp(coef)"), "Pr(>|z|)") %>% setNames(c("genus", "HR_likely", "pvalue_likely"))
df_cox_PC_s <- df_cox_PC %>% dplyr::select(c("genus" , "exp(coef)"), "Pr(>|z|)") %>% setNames(c("genus", "HR_PC", "pvalue_PC"))
df_cox_putative_s <- df_cox_putative %>% dplyr::select(c("genus" , "exp(coef)"), "Pr(>|z|)") %>% setNames(c("genus", "HR_putative", "pvalue_putative"))

df_cox_combine <- Reduce(function(x,y){merge(x,y,by="genus",all=TRUE)},
                         list(df_cox_without_s, df_cox_likely_s, df_cox_PC_s, df_cox_putative_s))
df_cox_combine$names <- str_extract(df_cox_combine$genus, "g__.+")
  
df_cox_combine <- df_cox_combine %>% remove_rownames() %>% column_to_rownames(var = "names")
  
df_HR_values <- df_cox_combine[,grep("HR",names(df_cox_combine))]
df_pvalues <- df_cox_combine[,grep("pvalue",names(df_cox_combine))]

pdf(file = "Figures/Heatmap_Hazard_Ration_35_biomarkers.pdf",height = 10,width = 5)
Heatmap(as.matrix(df_HR_values), name = "Hazard Ratio", 
        col = colorRamp2(c(0, 1, max(df_HR_values)),c("#0072B5CC", "#FFFFFF", "#E18727CC")), column_dend_height = unit(10, "mm"),
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "complete",
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        row_names_gp = gpar(fontsize=10),
        column_names_gp = gpar(fontsize=12),
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        row_title_side = "left",row_dend_side="left",
        row_names_side = "left", column_names_rot = 45,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(df_pvalues[i, j] < 0.001) {
            grid.text("***", x, y)
          } else if(df_pvalues[i, j] < 0.01) {
            grid.text("**", x, y)
          } else if(df_pvalues[i, j] < 0.05) {
            grid.text("*", x, y)
          }
        })
dev.off()

write.csv(df_cox_combine, file = "Tables/Cox_35_microbial_sigantures.csv", row.names = TRUE)
save(list = ls(all.names = TRUE), file = "RData/06.ACC_Microbial_Feature_TWO_Clusters.RData", envir = .GlobalEnv)


