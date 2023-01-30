
rm(list = ls())
getwd()
load('step00_output.Rdata')
{
  library(GSVA)
  library(ggpubr)
  library(tidyr)
  library(dplyr)
  library(magrittr)
}
filter_genes <- function(genes,qt){
  genes.t <- t(genes)
  genes.sd <- transform(as.data.frame(genes.t), SD=apply(as.data.frame(genes.t),1,sd, na.rm = TRUE))
  ## select top genes with high SD (~ variability) across samples
  SD_quantile <- quantile(genes.sd$SD) ## identical to summary(genes.sd$SD)
  SD_cutoff <- SD_quantile[qt] ## 2nd quantile -- 25th quantile.
  genes.sd <- genes.sd[order(genes.sd$SD, decreasing = T),]
  top.variable.genes <- rownames(genes.sd[genes.sd$SD > SD_cutoff,])
  ## subset these genes from gene table
  select <- which(colnames(genes) %in% top.variable.genes)
  genes <- genes[,select]
}
FUN =  function(xxx) {
  (sum(xxx==0)/length(xxx))<0.5 }

exp_tpm <- read.table(file = './03transcriptomic_analysis/source_data/TCGA-ACC.htseq_tpm.tsv',header = T,row.names = 1,sep = '\t',check.names = F) %>% 
  .[,clinic$PATIENT_ID] %>% 
  t() %>% 
  .[,apply(., 2, FUN)] %>% 
  filter_genes(2) 
dim(exp_tpm)
#[1]  77 14267

write.table(log2(t(tpm)+1),file = './03transcriptomic_analysis/output_acc_filtered_mrna_exp_log_transformed.tsv',sep = '\t',row.names = T,col.names = NA)

# ESTIMATE ----------------------------------------------------------------
# rforge <- "http://r-forge.r-project.org"
# install.packages('estimate',repos = rforge,dependencies = T)
library('estimate')
library(ggunchained)
filterCommonGenes(input.f = './03transcriptomic_analysis/output_acc_filtered_mrna_exp_log_transformed.tsv',output.f = './03transcriptomic_analysis/output.gct',id='GeneSymbol') 
estimateScore('./03transcriptomic_analysis/output.gct','./03transcriptomic_analysis/score.gct')
estimate_score <- read.table('./03transcriptomic_analysis/score.gct',skip = 2,header = T,row.names = 1,check.names = F)[1:2,] %>% 
  mutate(Description=NULL) %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = 'PATIENT_ID') %>% 
  mutate(PATIENT_ID=gsub('[.]','-',.$PATIENT_ID)) %>% 
  merge(.,clinic[,c("PATIENT_ID","Without_combined","Likely_combined","PC_combined","Putative_combined")],by='PATIENT_ID') %>% 
  gather(.,category,value,-PATIENT_ID,-Without_combined,-Likely_combined,-PC_combined,-Putative_combined) %>%
  gather(.,data_type,MS,-PATIENT_ID,-value,-category) %>%
  mutate(category=as.factor(.$category)) %>% 
  mutate(MS=as.factor(.$MS)) %>% 
  mutate(data_type=as.factor(.$data_type)) %>% 
  mutate(value=as.numeric(.$value)) %>% 
  mutate(value=scale(.$value))

pdf('./03transcriptomic_analysis/output_estimate_score_ms1_ms2.pdf')
ggplot(estimate_score, aes(x=category, y=value,fill=MS)) +
  geom_split_violin(trim = FALSE,color = NA)+
  stat_summary(fun = mean,
               fun.min = function(x){quantile(x)[2]},
               fun.max = function(x){quantile(x)[4]},
               geom = "pointrange",
               #geom = 'errorbar',
               size=0.5,
               position = position_dodge(width = 0.2))+
  facet_wrap(~data_type,nrow = 2) +
  theme_bw()+
  thm+
  scale_fill_manual(values = c(MS1="#0072B5CC",MS2="#E18727CC"))+
  guides(fill=guide_legend(direction = "horizontal",nrow =1))+
  # ylim(0,10)+
  labs(x="", y = "Estimate Score (Normalized)", fill = "")+  
  stat_compare_means(aes(x = category , y = value,group=MS),label = "p.signif",size=8)
dev.off()


# (新增)TIL score and immune gene score ----------------------------------------------

# immune_infiltration_acc <- read.csv(file = './data/transcriptomic_analysis/infiltration_estimation_for_tcga.csv',header = T,row.names = 1,check.names = F) %>% 
#   .[substr(rownames(.),13,15)=='-01',] %>% 
#   set_rownames(substr(rownames(.),1,12)) %>% 
#   .[rownames(.)%in%clinic$PATIENT_ID,] 
#   # .[,str_detect(colnames(immune_infiltration_acc),'XCELL')]
# immune_infiltration_acc <- log2(immune_infiltration_acc+1)
# acc_TIL <- immune_infiltration_acc %>% 
#   tibble::rownames_to_column(var = 'PATIENT_ID') %>% 
#   merge(clinic[,c("PATIENT_ID","kmeans")],by='PATIENT_ID') %>% 
#   gather(TIL,value,-PATIENT_ID,-kmeans)
# TIL_cell <- unique(acc_TIL$TIL)
# p_value_1=c()
# for (i in 1:length(TIL_cell)){
#   data<-acc_TIL[acc_TIL$TIL==TIL_cell[i],]
#   P=wilcox.test(value~kmeans,data)
#   p_value_1[i]=P$p.value
#   p_value_1[i]<- round(p_value_1[i],3) #设置小数点位数
# }
# # p_value_1[p_value_1<0.05]
# acc_TIL_sig <- acc_TIL[acc_TIL$TIL%in%TIL_cell[p_value_1<0.05],]
# ggplot(acc_TIL_sig, aes(x=TIL, y=value,fill=kmeans)) +
#   geom_boxplot(aes(fill=kmeans))+ #绘制箱线图
#   theme_bw()+
#   scale_fill_manual(values = c(MS1="#0072B5CC",MS2="#E18727CC"))+
#   guides(fill=guide_legend(direction = "horizontal",nrow =1))+
#   labs(x="", y = "TIL", fill = "")+  # fill为修改图例标题
#   stat_compare_means(aes(x = TIL , y = value,group=kmeans),method = 'wilcox.test',label = "p.signif")


# exp <- read.table(file = './data/TCGA-ACC.htseq_tpm.tsv',header = T,row.names = 1,sep = '\t',check.names = F) %>% 
#   .[,clinic$PATIENT_ID]
# exp <- as.data.frame(log2(exp+1))


# 免疫相关基因表达 ----------------------------------------------------------------
thm <- theme_bw()+
  theme(panel.border = element_rect(),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black",size=0.5),
        text = element_text(size = 20),
        axis.text.x = element_text(size=25,color ='black'),
        axis.text.y = element_text(size = 20,color = 'black'),
        axis.title.y = element_text(size = 20),
        plot.subtitle  = element_text(face = 'italic',hjust = 0.5),
        legend.position='right',
        legend.text = element_text(size = 15,hjust = 0),
        legend.title = element_text(size = 20)
  )

##免疫抑制基因来源TIP:http://biocc.hrbmu.edu.cn/TIP/
immuno_inhibitory_gene_matrix <- read.table(file = './03transcriptomic_analysis/source_data/immune_signature_annotation.txt',sep = '\t',header = T)
immuno_inhibitory_gene <- intersect(unique(immuno_inhibitory_gene_matrix$GeneSymbol),colnames(exp_tpm))
exp_immuno <- log2(exp_tpm+1) %>% 
  .[,immuno_inhibitory_gene] %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = 'PATIENT_ID') %>% 
  merge(clinic[,c("PATIENT_ID","Without_combined","Likely_combined","PC_combined","Putative_combined")],by='PATIENT_ID') %>% 
  gather(.,gene,value,-PATIENT_ID,-Without_combined,-Likely_combined,-PC_combined,-Putative_combined)

pvalue_table=data.frame(NR=NA,LR=NA,CR=NA,PR=NA)
data_type=c("Without_combined","Likely_combined","PC_combined","Putative_combined")
for (i in 1:length(immuno_inhibitory_gene)){
  for (j in 1:4) {
    data<-exp_immuno[exp_immuno$gene==immuno_inhibitory_gene[i],] %>% 
      dplyr::mutate(type=.[,data_type[j]])
    # filter(gene==immunosuppression_gene[i])
    P=wilcox.test(value~type,data)
    pvalue=P$p.value
    pvalue_table[i,j]<- round(pvalue,3) #设置小数点位数
  }
}
rownames(pvalue_table) <- immuno_inhibitory_gene
write.table(pvalue_table,file = './03transcriptomic_analysis/output_pvalue_immune_gene_&_MS_four_data_types.tsv',sep = '\t',row.names = T,col.names = NA)

immuno_gene_sig <- immuno_inhibitory_gene[pvalue_table$NR<0.05&pvalue_table$LR<0.05&pvalue_table$CR<0.05&pvalue_table$PR<0.05]
exp_immuno <- log2(exp_tpm+1) %>% 
  .[,immuno_inhibitory_gene] %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = 'PATIENT_ID') %>% 
  merge(clinic[,c("PATIENT_ID","Without_combined","Likely_combined","PC_combined","Putative_combined")],by='PATIENT_ID') %>% 
  gather(.,gene,value,-PATIENT_ID,-Without_combined,-Likely_combined,-PC_combined,-Putative_combined) %>% 
  gather(.,data_type,MS,-PATIENT_ID,-value,-gene) %>%
  mutate(gene=as.factor(.$gene)) %>%
  mutate(MS=as.factor(.$MS)) %>% 
  dplyr::mutate(data_type=factor(.$data_type,levels = c("Without_combined","Likely_combined","PC_combined","Putative_combined")))

pdf('./03transcriptomic_analysis/output_immune_genes_&_MS_four_data_types.pdf',width = 20,height = 15)
ggplot(exp_immuno[exp_immuno$gene%in%immuno_gene_sig,], aes(x=gene, y=value,fill=MS)) +
  geom_split_violin(trim = FALSE,color = NA)+
  stat_summary(fun = mean,
               fun.min = function(x){quantile(x)[2]},
               fun.max = function(x){quantile(x)[4]},
               geom = "pointrange",
               #geom = 'errorbar',
               size=0.1,
               position = position_dodge(width = 0.2))+
  facet_wrap(~data_type,nrow = 4) +
  theme_bw()+
  thm+
  scale_fill_manual(values = c(MS1="#0072B5CC",MS2="#E18727CC"))+
  guides(fill=guide_legend(direction = "horizontal",nrow =1))+
  # ylim(0,10)+
  labs(x="", y = "Estimate Score (Normalized)", fill = "")+  
  stat_compare_means(aes(x = gene , y = value,group=MS),label = "p.signif",size=10)
dev.off()

# immune_gene_table <- read.table(file = './data/transcriptomic_analysis/GeneList.txt',sep = '\t',header = T) %>% 
#   filter(Category=='Antimicrobials')
# immune_gene <- unique(immune_gene_table$Symbol) %>% 
#   intersect(.,rownames(exp))
# acc_immune <- exp %>% 
#   .[immune_gene,] %>% 
#   t() %>% 
#   as.data.frame() %>% 
#   tibble::rownames_to_column(var = 'PATIENT_ID') %>% 
#   merge(clinic[,c("PATIENT_ID","kmeans")],by='PATIENT_ID') %>% 
#   gather(gene,value,-PATIENT_ID,-kmeans)
# p_value=c()
# for (i in 1:length(immune_gene)){
#   data<-acc_immune %>% 
#     filter(gene==immune_gene[i])
#   P=wilcox.test(value~kmeans,data)
#   p_value[i]=P$p.value
#   # p_value[i]<- round(p_value[i],3) #设置小数点位数
#  }
# p_value[p_value<0.05]
# immune_combine_1 <- t(immune_infiltration_acc) %>% as.data.frame() %>% .[unique(acc_TIL$TIL)[p_value_1<0.05],] %>% scale(scale = T,center = T)
# immune_combine_2 <- exp[immunosuppression_gene[p_value_2<0.05],] %>% na.omit() %>% scale(scale = T,center = T)
# immune_combine <- rbind(immune_combine_1,immune_combine_2)





# GSVA_22_immune_cells ----------------------------------------------------
# rm(list = ls())
{
  library(genefilter)
  library(GSVA)
  library(Biobase)
  library(stringr)
}
load('step00_output.Rdata')
exp_tpm <- read.table(file = './03transcriptomic_analysis/source_data/TCGA-ACC.htseq_tpm.tsv',header = T,row.names = 1,sep = '\t',check.names = F) %>% 
  .[,clinic$PATIENT_ID] %>% 
  t() %>% 
  .[,apply(., 2, FUN)] %>% 
  filter_genes(2) 
dim(exp_tpm)

uni_matrix <- as.matrix(log2(t(exp_tpm)+1))
gene_set<-read.csv("./03transcriptomic_analysis/source_data/gsva_immune_cells.csv")[, 1:2]#https://www.cell.com/cms/10.1016/j.celrep.2016.12.019/attachment/f353dac9-4bf5-4a52-bb9a-775e74d5e968/mmc3.xlsx
head(gene_set)
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
gsva_immu_cell<- gsva(as.matrix(uni_matrix), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE) #log值选择‘Gaussian’
gsva_immu_cell <- apply(gsva_immu_cell, 2, function (x) (max(x)-x)/(max(x)-min(x)))
# gsva_immu_cell <- log2(gsva_immu_cell+1)
immune_infiltrating <- gsva_immu_cell %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(.,'PATIENT_ID') %>% 
  merge(clinic[,c("PATIENT_ID","Without_combined","Likely_combined","PC_combined","Putative_combined")],by='PATIENT_ID') %>% 
  gather(.,immune_cell,value,-PATIENT_ID,-Without_combined,-Likely_combined,-PC_combined,-Putative_combined)


pvalue_table=data.frame(NR=NA,LR=NA,CR=NA,PR=NA)
logFC_table=data.frame(NR=NA,LR=NA,CR=NA,PR=NA)
data_type=c("Without_combined","Likely_combined","PC_combined","Putative_combined")
for (i in 1:length(unique(immune_infiltrating$immune_cell))){
  for (j in 1:4) {
    data<-immune_infiltrating[immune_infiltrating$immune_cell==unique(immune_infiltrating$immune_cell)[i],] %>% 
      dplyr::mutate(type=.[,data_type[j]])
    test=wilcox.test(value~type,data)
    FC=mean(filter(data,type=='MS1')[,'value'])/mean(filter(data,type=='MS2')[,'value'])
    logFC_table[i,j] <- log2(FC+1)
    pvalue=test$p.value
    pvalue_table[i,j]<- round(pvalue,3) #设置小数点位数
  }
}
rownames(logFC_table) <- unique(immune_infiltrating$immune_cell)
rownames(pvalue_table) <- unique(immune_infiltrating$immune_cell)
pvalue_table <- ifelse(pvalue_table<0.01,'**',ifelse(pvalue_table<0.05,'*',''))


# write.table(pvalue_table,file = './03transcriptomic_analysis/output_pvalue_immune_infiltrating_&_MS_four_data_types.tsv',sep = '\t',row.names = T,col.names = NA)



# boxplot
# pdf('gsva_22_immune_cells_ms1_ms2.pdf',width = 20,height = 10)
# ggplot(gsva_immu_box, aes(x=immu_cell, y=value,fill=kmeans)) +
#   geom_boxplot(width=0.4,position=position_dodge(0.8),aes(fill=kmeans))+ #绘制箱线图
#   theme_bw()+
#   theme(panel.border = element_rect(),
#         panel.grid.major = element_blank(),   
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black",size=0.5),
#         text = element_text(size = 20),
#         axis.text.x = element_text(size=15,color ='black',angle = 90,hjust = 1),
#         axis.text.y = element_text(size = 10,color = 'black'),
#         axis.title.y = element_text(size = 15),
#         legend.text = element_text(size = 15,hjust = 0),
#         legend.position = 'none',
#         plot.margin=margin(t = 0, r = 1, b = 0, l = 1, unit = "cm"))+
#   scale_fill_manual(values = c(MS1="#0072B5CC",MS2="#E18727CC"))+
#   guides(fill=guide_legend(direction = "horizontal",nrow =1))+
#   labs(x="", y = "GSVA Score Immune Cells Proportion %", fill = "")+  # fill为修改图例标题
#   stat_compare_means(aes(x = immu_cell , y = value,group=kmeans),label = "p.signif",size=5)
# dev.off()

# heatmap plot
{
  library(ComplexHeatmap)
  library(circlize)
}

ht_matrix <- as.matrix(t(logFC_table)) 
max(ht_matrix)
# [1] 4.925572
min(ht_matrix)
# [1] -1.703126

# clinic <- clinic[order(clinic$Without_combined),] #clinic按照Without_combined排序
# ht_matrix <- ht_matrix[,] #矩阵按照顺序排列



# colomn_group <- c(rep('MS1',34),rep('MS2',43))#列分组信息
# names(colomn_group) <- clinic$PATIENT_ID

# # top annotation
# top_anno = HeatmapAnnotation(
#   cluster = anno_block(gp = gpar(fill = c("#0072B5CC", "#E18727CC")),
#                        labels = c("MS1","MS2"),
#                        labels_gp = gpar(col = "white", fontsize = 15)))

library(circlize)
col_fun = colorRamp2(c(0.5, 1, 1.5), c("#E46726", "white", "#6D9EC1"))
cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(
    t(pvalue_table)[i, j], 
    x, y,
    gp = gpar(
      fontsize = 20
    ))
}

# plot heatmap
pdf('./03transcriptomic_analysis/output_heatmap_gsva_22_immune_cells_ms1_ms2.pdf',width = 20,height = 5)
Heatmap(ht_matrix,
        col = col_fun,
        cell_fun = cell_fun,
        cluster_rows = FALSE,
        show_column_dend = FALSE,
        show_row_names = T,
        column_names_rot = 90,
        show_column_names = T,
        column_names_gp = gpar(fontsize=10),
        # heatmap_width = unit(15,'cm'),
        heatmap_legend_param = list(title='LogFC(MS1/MS2)',title_position = "topcenter",legend_direction='horizontal'),
        rect_gp = gpar(col = "black", lty =1, lwd = 0.5),
        column_title_gp = gpar(fontsize=15))
dev.off()



# Heatmap(ht_matrix,
#         show_column_names = F,
#         row_names_gp = gpar(fontsize=15),
#         column_split = 2,
#         heatmap_width = unit(20,'cm'),
#         cluster_columns =cluster_within_group(ht_matrix,colomn_group),
#         show_heatmap_legend = F,
#         show_column_dend = F,
#         column_title = NULL,
#         top_annotation = top_anno)
# dev.off()

# # Immune checkpoint genes -------------------------------------------------
# 
# immu_check_genes <- c('PDCD1','CD274','CTLA4','TIGIT')
# exp <- read.table(file = './data/TCGA-ACC.htseq_tpm.tsv',header = T,row.names = 1,sep = '\t',check.names = F) %>% 
#   .[,clinic$PATIENT_ID] 
# exp <- as.data.frame(log2(exp+1))
# acc_immu_exp <- exp %>% 
#   .[immu_check_genes,] %>% 
#   t() %>% 
#   as.data.frame() %>% 
#   tibble::rownames_to_column(var = 'PATIENT_ID') %>% 
#   merge(clinic[,c("PATIENT_ID","kmeans")],by='PATIENT_ID') %>% 
#   gather(gene,value,-PATIENT_ID,-kmeans)
# 
# # plot
# pdf('./data/transcriptomic_analysis/immunecheck_molecules_ms1_ms2.pdf')
# ggplot(acc_immu_exp, aes(x=gene, y=value,fill=kmeans)) +
#   geom_boxplot(width=0.6,position=position_dodge(0.8),aes(fill=kmeans))+ #绘制箱线图
#   theme_bw()+
#   theme(panel.border = element_rect(),
#         panel.grid.major = element_blank(),   
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black",size=0.5),
#         text = element_text(size = 20),
#         axis.text.x = element_text(size=20,angle = 45,vjust = 0.5,color ='black'),
#         axis.text.y = element_text(size = 10,color = 'black'),
#         axis.title.y = element_text(size = 20),
#         legend.text = element_text(size = 20,hjust = 0),
#         legend.title = element_text(size = 20),
#         legend.position = 'top',
#         plot.margin=margin(t = 0, r = 2, b = 0, l = 2, unit = "cm"))+
#   scale_fill_manual(values = c(MS1="#0072B5CC",MS2="#E18727CC"))+
#   guides(fill=guide_legend(direction = "horizontal",nrow =1))+
#   labs(x="", y = "Gene Expression/log2(TPM+1)", fill = "")+  
#   stat_compare_means(aes(x=gene, y=value,group=kmeans),label = 'p.signif',size=8,label.y = 5)
# dev.off()
# 
# 


# DEG analysis ------------------------------------------------------------
rm(list = ls())
{
  library(tidyverse)
  library(ggplot2)
  library(cowplot)
  library(ggrepel)
  library(ggstatsplot)
  library(limma)
  library(ggsci)
  library(edgeR)
  library(org.Hs.eg.db)
  library(clusterProfiler)
}
load('step00_output.Rdata')

filter_genes <- function(genes,qt){
  genes.t <- t(genes)
  genes.sd <- transform(as.data.frame(genes.t), SD=apply(as.data.frame(genes.t),1,sd, na.rm = TRUE))
  ## select top genes with high SD (~ variability) across samples
  SD_quantile <- quantile(genes.sd$SD) ## identical to summary(genes.sd$SD)
  SD_cutoff <- SD_quantile[qt] ## 2nd quantile -- 25th quantile.
  genes.sd <- genes.sd[order(genes.sd$SD, decreasing = T),]
  top.variable.genes <- rownames(genes.sd[genes.sd$SD > SD_cutoff,])
  ## subset these genes from gene table
  select <- which(colnames(genes) %in% top.variable.genes)
  genes <- genes[,select]
}
FUN =  function(xxx) {
  (sum(xxx==0)/length(xxx))<0.5 }
exp_count <- read.table(file = './03transcriptomic_analysis/source_data/TCGA-ACC.htseq_counts.tsv',header = T,row.names = 1,sep = '\t',check.names = F) %>% 
  .[,clinic$PATIENT_ID] %>% 
  t() %>% 
  .[,apply(., 2, FUN)] %>% 
  filter_genes(2) 
dim(exp_count)
#[1]  77 14267
exp <- t(exp_count)

data_type <- c("Without_combined","Likely_combined","PC_combined","Putative_combined")
# kegg_result <- data.frame('ID'=NA,"Description"=NA,"GeneRatio"=NA,"BgRatio"=NA, "pvalue"=NA,"p.adjust"=NA,"qvalue"=NA,"geneID"=NA,"Count"=NA,"type"=NA)
kegg_result <- data.frame()
gsea_result <- data.frame()
# gsea_result <- data.frame( "ID"=NA,"Description"=NA,"setSize"=NA,"enrichmentScore"=NA,"NES"=NA,"pvalue"=NA,"p.adjust"=NA,"qvalue"=NA,"rank"=NA,"leading_edge"=NA,"core_enrichment"=NA,"type"=NA)
gsea_result_list <- list()
for (i in 1:4) {
  if(T){
    clinic$MS <- clinic[,data_type[i]]
    design <- model.matrix(~0+factor(clinic$MS)) %>%  # group matrix 
      set_colnames(c(levels(factor(clinic$MS)))) %>% 
      set_rownames(colnames(exp))
    y = voom(exp, design, plot = T) #important！
    contrast.matrix<-makeContrasts(paste0(c('MS1','MS2'),collapse = "-"),levels = design) 
    contrast.matrix
    fit <- lmFit(y,design) %>%  #step1 
      contrasts.fit(., contrast.matrix) %>% # step2 important!
      eBayes(.) # default no trend !!!
    nrDEG = topTable(fit, coef=1, n=Inf) %>% # step3
      na.omit(tempOutput) 
    head(nrDEG)
    
    logFC_cutoff <- 1  #logFC_cutoff <- with(DEG,mean(abs( logFC)) + 2*sd(abs( logFC)))
    DEG=nrDEG %>% 
      dplyr::mutate(change=as.factor(ifelse(.$P.Value < 0.05 & abs(.$logFC) > logFC_cutoff,
                                            ifelse(.$logFC > logFC_cutoff ,'UP','DOWN'),'NOT'))) %>% 
      dplyr::mutate(GENE_NAME=rownames(.))
    
    #KEGG
    
    DEG <- tibble::rownames_to_column(DEG,var = 'Gene_symbol')
    dff <- bitr(unique(DEG$Gene_symbol), fromType = "SYMBOL",
                toType = c( "ENTREZID"),
                OrgDb = org.Hs.eg.db)
    head(dff)
    DEG=merge(DEG,dff,by.x='Gene_symbol',by.y='SYMBOL')
    head(DEG)
    gene_up= DEG[DEG$change == 'UP','ENTREZID'] 
    kk.up <- enrichKEGG(gene         = gene_up,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05)
    # head(kk.up)[,1:6]
    x <- kk.up@result[1:10,] %>% 
      dplyr::mutate(type=rep(data_type[i]))
    kegg_result <- rbind(kegg_result,x)
    
    #GSEA
    geneList<-DEG$logFC
    names(geneList)=DEG$ENTREZID #使用转换好的ID
    geneList=sort(geneList,decreasing = T) #从高到低排序
    set.seed(123)
    KEGG_gseresult <- gseKEGG(geneList, keyType='kegg',organism = "hsa",nPerm = 1000, minGSSize = 25, maxGSSize = 300,pvalueCutoff = 0.1, pAdjustMethod = "BH",use_internal_data=FALSE,seed = T)
    KEGG_gseresult2 <- setReadable(KEGG_gseresult,     #将结果中的ENTREZID转化为symbol
                     OrgDb = "org.Hs.eg.db",
                     keyType = "ENTREZID")
    gsea_result_list[[i]] <- KEGG_gseresult2
    y <- KEGG_gseresult2 %>% .@result %>% 
      dplyr::mutate(type=rep(data_type[i]))
    gsea_result <- rbind(gsea_result,y)
    # y <- KEGG_gseresult@result
    write.table(y,file = paste0('./03transcriptomic_analysis/output_GSEA_result_',data_type[i],'.tsv'),sep = '\t',row.names = T,col.names = NA)
    
  }
}


# kegg plot
result1 <- kegg_result %>% .[-1,] %>% dplyr::mutate(type=factor(.$type,levels = c("Without_combined","Likely_combined","PC_combined","Putative_combined" )))
# com_path <- intersect(result1$Description[result1$type==data_type[1]],intersect(result1$Description[result1$type==data_type[2]],intersect(result1$Description[result1$type==data_type[3]],result1$Description[result1$type==data_type[4]])))
# result1 <- result1[result1$Description%in%com_path,]
pdf('./03transcriptomic_analysis/output_kegg_top10_four_data_types.pdf',width = 12,height = 15)
ggplot(result1,aes(x = type, y =reorder(Description,Count)))+ 
  geom_point(aes(size=Count,color=pvalue))+
  scale_size_continuous(range=c(4,10))+
  scale_colour_gradient(low="#f68084",high="#a6c0fe")+
  labs(
    color=expression(pvalue),
    size=" Count",
    x="",
    y='Kegg pathways'
  )+
  theme_bw()+
  thm
dev.off()
# gsea plot
result2 <- gsea_result %>% .[-1,]
{
  library(enrichplot)
  library(ggsci)
}

# for (j in 1:4) {
#   data <- filter(result2,type==data_type[1]) %>%
#     dplyr::filter(pvalue<0.05) %>%
#     dplyr::arrange(abs(NES))
#   pathway <- data$Description[1:5]
#   labels <- c()
#   for (i in 1:length(pathway)) {
#     z <- grep(pathway[i],result2[result2$type==data_type[1],'Description'])
#     labels <- c(labels,z)
#   }
#   pdf(paste0('./03transcriptomic_analysis/output_gsea_top5_',data_type[j],'.pdf'))
#   p <- gseaplot2(gsea_result_list[[j]],
#             labels,
#             pvalue_table = FALSE,
#             title = paste0('GSEA_',data_type[j]),
#             base_size = 10,
#             rel_heights = c(1.5,0.5,1),
#             color = pal_npg()(length(pathway)),
#             ES_geom = "line")
#   print(p)
#   dev.off()
# }

pathway=c('Cell cycle','p53 signaling pathway','TGF-beta signaling pathway','ECM-receptor interaction')
pdf('./03transcriptomic_analysis/output_gsea_p53_signaling_pathway.pdf')
gseaplot2(gsea_result_list[[3]],
          c(grep(pathway[2],gsea_result_list[[3]]@result[,'Description'])),
          pvalue_table = TRUE,
          title = paste0('GSEA_',data_type[3]),
          base_size = 10,
          rel_heights = c(1.5,0.5,1),
          color = pal_npg()(1),
          ES_geom = "line")
dev.off()
NES <- gsea_result_list[[3]]@result[grep(pathway[2],gsea_result_list[[3]]@result[,'Description']),'NES'];NES

