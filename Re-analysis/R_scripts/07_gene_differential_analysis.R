rm(list = ls())
getwd()

{
  library(GSVA)
  library(ggpubr)
  library(tidyr)
  library(dplyr)
  library(magrittr)
  library(Biobase)
  library(openxlsx)
  library(stringr)
  library(IOBR)
  library(clusterProfiler)
  library(tidyverse)
  library(ggplot2)
  library(cowplot)
  library(ggrepel)
  library(ggstatsplot)
  library(limma)
  library(edgeR)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(UpSetR)
  library(GseaVis)
  library(ggthemes)
  library(igraph)
  library(Hmisc)
  library(ggraph)
}


### load data

# gene expression
gene_tpm <- read.xlsx( './00_data_preparation/gene_table.xlsx',sheet = 1,rowNames = TRUE) %>%   mutate_all(as.numeric);dim(gene_tpm)
# [1] 25531    79

# clinical info
df_meta_acc_tcga <- read.csv("./00_data_preparation/metadata.csv",na.strings = "") %>%  filter(COHORT=='TCGA')

df_clinical_cluster <- read.csv("./clustering/whole/pam/Unsupervised_clustering.csv") %>% 
  dplyr::select(c("PATIENT_ID", "genus_tcga_s_euclidean")) %>% 
  setNames(c("PATIENT_ID", "MS")) %>% 
  dplyr::mutate(MS=recode(.$MS,'1'="MS1",'2'="MS2")) %>% 
  dplyr::mutate(MS = as.factor(MS)) %>% 
  merge(., df_meta_acc_tcga, by="PATIENT_ID", all = TRUE) %>% 
  set_rownames(.$PATIENT_ID)

### DEG analysis 
fun_DEG_analysis <- function(clinic,exp){
  exp <- exp[,clinic$PATIENT_ID]
  if(identical(clinic$PATIENT_ID,colnames(exp))){
    design <- model.matrix(~0+factor(clinic$MS)) %>%  # group matrix 
      set_colnames(c(levels(factor(clinic$MS)))) %>% 
      set_rownames(colnames(exp))
    y = voom(exp, design, plot = T) # important！
    contrast.matrix<-makeContrasts(paste0(c('MS1','MS2'),collapse = "-"),levels = design) 
    contrast.matrix
    fit <- lmFit(y,design) %>%  #step1 
      contrasts.fit(., contrast.matrix) %>% # step2 important!
      eBayes(.)  # default no trend !!!
    nrDEG = topTable(fit, coef=1, n=Inf) %>% # step3
      na.omit(tempOutput) 
    head(nrDEG)
    
    logFC_cutoff <- 0.5
    logFC_cutoff <- with(nrDEG,mean(abs( logFC)) + 2*sd(abs( logFC)))
    DEG=nrDEG %>% 
      dplyr::mutate(change=as.factor(ifelse(.$P.Value < 0.05 & abs(.$logFC) > logFC_cutoff,
                                            ifelse(.$logFC > logFC_cutoff ,'UP','DOWN'),'NOT'))) %>% 
      dplyr::mutate(GENE_NAME=rownames(.))
    table(DEG$change)
    
    #KEGG
    DEG <- tibble::rownames_to_column(DEG,var = 'Gene_symbol')
    dff <- bitr(unique(DEG$Gene_symbol), fromType = "SYMBOL",
                toType = c( "ENTREZID"),
                OrgDb = org.Hs.eg.db)
    head(dff)
    DEG=merge(DEG,dff,by.x='Gene_symbol',by.y='SYMBOL')
    head(DEG)
    gene_change= DEG[DEG$change %in% c('UP','DOWN'),'ENTREZID']
    kk.change <- enrichKEGG(gene= gene_change,
                            organism = 'hsa',
                            pvalueCutoff = 0.05)
    kegg_result <- kk.change@result 
    
    # GSEA
    geneList<-DEG$logFC
    names(geneList)=DEG$ENTREZID #使用转换好的ID
    geneList=sort(geneList,decreasing = T) #从高到低排序
    set.seed(123)
    KEGG_gsearesult <- gseKEGG(geneList, keyType='kegg',organism = "hsa",nPerm = 1000, minGSSize = 25, maxGSSize = 300,pvalueCutoff = 0.1, pAdjustMethod = "BH",use_internal_data=FALSE,seed = T)
    KEGG_gsearesult2 <- setReadable(KEGG_gsearesult,     #将结果中的ENTREZID转化为symbol
                                   OrgDb = "org.Hs.eg.db",
                                   keyType = "ENTREZID")
    gsea_result <- KEGG_gsearesult2 %>% .@result 
    
    # output
    return(list(logFC_cutoff,
                DEG,
                kk.change,
                kegg_result,
                KEGG_gsearesult,
                gsea_result))
  }
}
res_DEG <- fun_DEG_analysis(df_clinical_cluster,gene_tpm)

# volcano plot
logFC_cutoff <- res_DEG[[1]]
df_volcano <- res_DEG[[2]]
key_genes <- df_volcano[df_volcano$change != "NOT", ]
# pdf('./gene_differential_analysis/volcano_deg_analysis.pdf',width = 4,height = 5,onefile = FALSE)
ggplot(df_volcano, aes(x = logFC, y = -log10(P.Value), color = change)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("NOT" = "grey", 
                                "UP" = "#E18727CC", 
                                "DOWN" = "#0072B5CC")) +
  labs(title = "Volcano Plot",x = "Log2 Fold Change",y = "-Log10(P-value)") +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = key_genes, aes(label = GENE_NAME), size = 3) +
  theme_few() +
  theme(
    legend.title = element_blank(),
    legend.position = "top"
  )
# dev.off()


# gsea plot
gsea_result <- res_DEG[[5]]
gsea_terms <- gsea_result@result[gsea_result@result$p.adjust < 0.05, ]
terms <- c('hsa04110','hsa04610','hsa05150','hsa00350','hsa00830','hsa00140')
gseaList <- lapply(terms, function(x){
  gseaNb(object = gsea_result,
         geneSetID = x,
         addPval = T,
         pvalX = 0.5,
         pvalY = 0.75,
         pCol = 'black',
         pHjust = 0)
}) 

# pdf('./gene_differential_analysis/signaling_pathways.pdf',width = 12,height = 8)
cowplot::plot_grid(plotlist = gseaList,ncol = 3,align = 'hv')
# dev.off()

# netplot top10
geneList<-res_DEG[[2]]$logFC
names(geneList)=res_DEG[[2]]$ENTREZID #使用转换好的ID
geneList=sort(geneList,decreasing = T) #从高到低排序
# pdf('./gene_differential_analysis/netplot_signaling_pathways.pdf',width = 6,height = 6)
cnetplot(res_DEG[[5]], 
         layout = igraph::layout_with_fr,
         node_label='category',
         colorEdge=TRUE,
         colorCategory = "p.adjust",
         showCategory = 10, foldChange = geneList, circular = TRUE)
# dev.off()

# interplay plot
otu_data <- read.xlsx('./00_data_preparation/otu_genus_cpm_table.xlsx',sheet = 2,rowNames = TRUE)  %>%  mutate(across(where(is.character), as.numeric)) %>%log1p() %>% t()
gene_data <- read.xlsx( './00_data_preparation/gene_table.xlsx',sheet = 1,rowNames = TRUE) %>% mutate_all(as.numeric) %>% t()

otu_feature <- read.csv(file = './microbial_feature/Core_features_coefficients.csv',header = TRUE) %>% dplyr::filter(pval<0.05) %>% .$feature %>% unique()
# gene_feature <- res_DEG[[6]]$Description[ res_DEG[[6]]$ID%in%c('hsa04110','hsa04610','hsa05150','hsa00350','hsa00830','hsa00140')]
gene_feature_hsa04110 <- unlist(strsplit(res_DEG[[6]]$core_enrichment[res_DEG[[6]]$ID=='hsa04110'], split = "/"))
gene_feature_hsa04610 <- unlist(strsplit(res_DEG[[6]]$core_enrichment[res_DEG[[6]]$ID=='hsa04610'], split = "/"))
gene_feature_hsa05150 <- unlist(strsplit(res_DEG[[6]]$core_enrichment[res_DEG[[6]]$ID=='hsa05150'], split = "/"))
gene_feature_hsa00350 <- unlist(strsplit(res_DEG[[6]]$core_enrichment[res_DEG[[6]]$ID=='hsa00350'], split = "/"))
gene_feature_hsa00830 <- unlist(strsplit(res_DEG[[6]]$core_enrichment[res_DEG[[6]]$ID=='hsa00830'], split = "/"))
gene_feature_hsa00140 <- unlist(strsplit(res_DEG[[6]]$core_enrichment[res_DEG[[6]]$ID=='hsa00140'], split = "/"))

fun_interplay_plot <- function(otu_data,otu_feature,gene_data,gene_feature){
  
  otu_matrix <- as.matrix(otu_data[,otu_feature])
  gene_matrix <- as.matrix(gene_data[,gene_feature])
  
  # 计算相关性矩阵和p值（同之前的示例）
  correlation_matrix <- Hmisc::rcorr(otu_matrix, gene_matrix, type = "spearman")
  
  # 提取相关系数和p值
  r_values <- correlation_matrix$r[1:ncol(otu_matrix), (ncol(otu_matrix)+1):(ncol(otu_matrix)+ncol(gene_matrix))]
  p_values <- correlation_matrix$P[1:ncol(otu_matrix), (ncol(otu_matrix)+1):(ncol(otu_matrix)+ncol(gene_matrix))]
  
  # correlation_values <- correlation_matrix$r
  # p_values <- correlation_matrix$P
  
  # 4. 多重检验校正（FDR）
  p_adjusted <- matrix(p.adjust(p_values, method = "fdr"), 
                       nrow = nrow(p_values), ncol = ncol(p_values))
  rownames(p_adjusted) <- rownames(r_values)
  colnames(p_adjusted) <- colnames(r_values)
  
  # 提取显著相关性
  edges <- which(abs(r_values)>0.3 & p_values < 0.05, arr.ind = TRUE)
  edge_list <- data.frame(
    from = rownames(r_values)[edges[,1]],
    to = colnames(r_values)[edges[,2]],
    cor = r_values[edges],
    p = p_adjusted[edges]
  )
  
  # significant_results <- data.frame(
  #   OTU = colnames(otu_matrix)[significant_correlations[, "row"]],
  #   Gene = colnames(gene_data)[significant_correlations[, "col"]],
  #   Correlation = correlation_values[significant_correlations],
  #   P_value = p_values[significant_correlations]
  # ) 
  
  # 创建边列表
  # edges <- data.frame(
  #   from = paste0("OTU:", significant_results$OTU),
  #   to = paste0("Gene:", significant_results$Gene),
  #   weight = abs(significant_results$Correlation)
  # )
  
  # 创建网络图对象
  network_graph <- graph_from_data_frame(edge_list, directed = FALSE)
  
  # 设置节点颜色和大小
  V(network_graph)$shape <- ifelse(grepl("g__", V(network_graph)$name), "circle", "square")
  V(network_graph)$color <- ifelse(grepl("g__", V(network_graph)$name), "#00A1D5FF", "#DF8F44FF")
  V(network_graph)$size <- 8
  
  # 设置边的宽度和颜色
  E(network_graph)$width <- E(network_graph)$cor * 5
  
  get_edge_colors <- function(cor_values) {
    # 创建从蓝(负)到白(中性)到红(正)的渐变色
    # col_pal <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
    #                               "#F7F7F7",  # 中性白色
    #                               "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))
    col_pal <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

    # 将相关性值(-1到1)映射到颜色索引(1到n)
    col_index <- round((cor_values + 1)/2 * (100-1)) + 1
    colors <- col_pal[col_index]
    return(colors)
  }
  edge_colors <- get_edge_colors(edge_list$cor)
  
  
  E(network_graph)$color <- edge_colors
  
  # E(network_graph)$color <- ifelse(significant_results$Correlation > 0, "#B22222", "#4682B4")
  
  # 绘制网络图
  plot(network_graph,
       layout = layout_in_circle(network_graph), # 使用Fruchterman-Reingold布局
       vertex.label.cex = 0.85,
       edge.curved = 0.2,
       vertex.frame.color = NA,
       vertex.label.color = "black",
       # vertex.label.color= ifelse(grepl("OTU", V(network_graph)$name), "#00A1D5FF", "#DF8F44FF"),
       # vertex.label.dist = 1.5,
       main = "OTU-Gene Correlation Network"
  )
  
  # 添加图例
  legend("bottomright",
         legend = c("OTU", "Gene"),
         fill = c("#00A1D5FF", "#DF8F44FF"),
         cex = 0.8
  )
  
  # 添加渐变色图例
  color_gradient <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
  legend_image <- as.raster(matrix(rev(color_gradient), ncol=1))

  # 图例位置和参数
  # plot(c(0,2), c(0,1), type = 'n', axes = F, xlab = '', ylab = '')
  rasterImage(legend_image, 0.9, 0.6, 1.0, 0.9)
  # text(x = 1.0, y = seq(0.2, 0.8, l = 5),
  #      labels = round(seq(-1, 1, length.out = 5), 1),
  #      cex = 0.8)
  text(x = 0.9, y = 1.0, labels = "Correlation", cex = 0.9)
  
  
  # p <- ggraph(network_graph, layout = 'fr') + 
  #   geom_edge_link(aes(color = weight, width = abs(weight)), 
  #                  alpha = 0.7, show.legend = TRUE) +
  #   geom_node_point(aes(shape = shape, color = color), size = 8) +
  #   geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  #   scale_edge_color_gradient2(low = "blue", mid = "gray", high = "red") +
  #   scale_shape_manual(values = c("circle" = 16, "square" = 15)) +
  #   scale_color_identity() +
  #   theme_graph() +
  #   labs(title = "OTU-Gene Network with ggraph")
  # print(p)
  
  
  
}

# pdf('./gene_differential_analysis/microbe_gene_interplay_plot.pdf',width = 6,height = 6)
fun_interplay_plot(otu_data,otu_feature,gene_data,gene_feature_hsa04110)
fun_interplay_plot(otu_data,otu_feature,gene_data,gene_feature_hsa04610)
fun_interplay_plot(otu_data,otu_feature,gene_data,gene_feature_hsa05150)
fun_interplay_plot(otu_data,otu_feature,gene_data,gene_feature_hsa00350)
fun_interplay_plot(otu_data,otu_feature,gene_data,gene_feature_hsa00830)
fun_interplay_plot(otu_data,otu_feature,gene_data,gene_feature_hsa00140)
# dev.off()


# pdf('./gene_differential_analysis/interplay_plot_for_abstract.pdf',width = 4,height = 4)
# fun_interplay_plot(otu_data,otu_feature,gene_data,gene_feature_hsa04110)
# dev.off()
