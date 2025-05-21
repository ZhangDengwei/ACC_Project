rm(list = ls())
getwd()

{
  library(AnnotationDbi)
  library(KEGGREST) #用于提取通路及基因信息
  library(msigdbr)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
}

# hsa_path<-keggLink("pathway","hsa")
# meta=unique(hsa_path)[grepl('hsa00',unique(hsa_path))]
# hsa_info<-lapply(meta,keggGet)
# 
# nm=unlist(lapply(hsa_info,function(x)x[[1]]$NAME))
# 
# genes=unlist(lapply(hsa_info,function(x){
#   g=x[[1]]$GENE
#   paste(str_split(g[seq(2,length(g),by=2)],';',simplify=T)[,1],collapse=';')
# }))
# df=data.frame(
#   hsa=meta,
#   nm=nm,
#   genes=genes
# )

gene_tpm <- read.xlsx( './00_data_preparation/gene_table.xlsx',sheet = 1,rowNames = TRUE) %>% 
  # .[rowMeans(.>0)>0.1,] %>% 
  mutate_all(as.numeric)
dim(gene_tpm)

# clinical data
df_meta_ACC_tcga <- read.csv("./00_data_preparation/metadata.csv",na.strings = "") %>% 
  filter(COHORT=='TCGA')

df_cluster_clinical <- read.csv("./clustering/whole/pam/Unsupervised_clustering.csv") %>% 
  dplyr::select(c("PATIENT_ID", "genus_tcga_s_euclidean")) %>% 
  setNames(c("PATIENT_ID", "MS")) %>% 
  dplyr::mutate(MS=recode(.$MS,'1'="MS1",'2'="MS2")) %>% 
  dplyr::mutate(MS = as.factor(MS)) %>% 
  merge(., df_meta_ACC_tcga, by="PATIENT_ID", all = TRUE) %>% 
  set_rownames(.$PATIENT_ID)



### microbial cluster and immune score
fun_metabolism_score <- function(exp){
  
  # exp: gene*sample, data.frame

  # gsva
  metabolic_genesets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% 
    dplyr::filter(grepl('hsa00|hsa01',gs_exact_source))
  metabolic_genes_by_pathway <- metabolic_genesets %>%
    dplyr::select(gene_symbol,gs_name)
  
  gene_set<-data.frame(gene=metabolic_genes_by_pathway$gene_symbol,type=metabolic_genes_by_pathway$gs_name)
  gene_set_list<- split(as.matrix(gene_set)[,1], gene_set[,2])
  gsvaPar <- gsvaParam(as.matrix(exp), gene_set_list,maxDiff = TRUE)
  metabolism_res <- gsva(gsvaPar) %>% 
    t() %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column('PATIENT_ID')
  
  return(metabolism_res)
}

res_metabolism <- fun_metabolism_score(gene_tpm)
metabolic_category <- read.csv(file = './source_data/kegg_metabolism_category.csv',header = TRUE)
df_cluster_clinical_metabolism <- merge(df_cluster_clinical,res_metabolism,by='PATIENT_ID') 


# heatmap plot
data <- df_cluster_clinical_metabolism %>% 
  dplyr::select("PATIENT_ID", contains(c('KEGG'))) %>%
  tibble::column_to_rownames('PATIENT_ID') %>% 
  t() 
range(data)
pd1 <- df_cluster_clinical_metabolism %>% 
  dplyr::select("PATIENT_ID", 'MS')%>% 
  dplyr::mutate(MS=factor(.$MS,levels=c('MS1','MS2')))
metabolic_category <- read.csv(file = './source_data/kegg_metabolism_category.csv',header = TRUE)
pd2 <- data.frame(gsva_score=rownames(data),
                  category=metabolic_category$category[match(rownames(data),metabolic_category$pathway_new)])

setequal(pd1$PATIENT_ID,colnames(data))

data_nor <- scale(t(data)) %>%  #scale默认对列标准化
  t() %>% 
  as.data.frame() %>% 
  dplyr::mutate(across(everything(), ~ ifelse(. > 1, 1, .))) %>% 
  dplyr::mutate(across(everything(), ~ ifelse(. < -1, -1, .)))
group <- pd1$MS
test_results <- apply(data_nor, 1, function(row) {
  group1_values <- row[group == "MS1"]
  group2_values <- row[group == "MS2"]
  test <-wilcox.test(group1_values, group2_values)
  test$p.value
})

significance_group <- ifelse(test_results < 0.05, "Significant", "Not Significant")
row_annotation <- rowAnnotation(p_value = anno_barplot(-log10(test_results)),
                                signif = significance_group,
                                col=list(signif=c("Significant" = "red", "Not Significant" = "#eaebea")))

color_fun <- colorRamp2(c(1,0,-1), c("#FFFF00", "white", "#2D004B"))
# color_fun <- colorRamp2(c(1,0.25,0,-0.25,-1), c(c('#08306C','#4778A9','#68AAD1','#9FC7E0','#C9DBEA')))
column_annotation <- HeatmapAnnotation(MS = pd1$MS,
                                       col = list(MS = c("MS1" = "#0072B5CC", "MS2" = "#E18727CC")))

# pdf('./cluster_and_immune/heatmap_cluster_and_metabolism.pdf',width = 15,height = 15)
ComplexHeatmap::Heatmap(data_nor, 
                        na_col = "white",
                        col = color_fun,  # 添加颜色映射函数
                        show_column_names = F,
                        row_names_side = "right",
                        name = "fraction",
                        column_split = pd1$MS, 
                        column_gap = unit(5, "mm"),
                        row_split = pd2$category,
                        row_gap = unit(5, "mm"),
                        column_title = NULL,
                        cluster_columns = TRUE,
                        cluster_rows = F, 
                        top_annotation = column_annotation,
                        right_annotation = row_annotation,
                        width = ncol(data_nor)*unit(3, "mm"), 
                        height = nrow(data_nor)*unit(4, "mm"),
                        row_names_gp = gpar(fontsize = 8)  # 调整行名字体大小
)
# dev.off()

