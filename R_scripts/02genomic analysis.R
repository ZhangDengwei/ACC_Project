rm(list = ls())
getwd()

# Landscape of mut and cnv(customed) --------------------------------------

load('step00_output.Rdata')
{
 require(maftools)
 library(tidyr)
  library(stringr)
 library(ggpubr)
  library(magrittr)
}

# require(TCGAbiolinks)
# library(tidyr)
# ACC_mutect2 <- GDCquery_Maf(tumor = "ACC", pipelines = "mutect2")
# write.csv(ACC_mutect2,file = 'ACC_mutect2.csv')

maf <- read.csv(file = './02genomic_analysis/source_data/ACC_mutect2.csv',header = T,row.names = 1,check.names = F) %>% 
  dplyr::mutate(Tumor_Sample_Barcode=str_sub(.$Tumor_Sample_Barcode,1,12)) %>% 
  .[.$Tumor_Sample_Barcode%in% clinic$PATIENT_ID,] #all patients=77

# clinical annotation
clinicalData <- clinic[,c('Without_combined',"Likely_combined","PC_combined","Putative_combined",'AGE_new','PATIENT_ID')]
colnames(clinicalData)[6] <- c("Tumor_Sample_Barcode")

# maf = read.maf(maf = mut,
#                clinicalData = clinicalData)
# plotmafSummary(maf = maf,   #使用plotmafSummary函数可视化maf对象汇总信息
#                rmOutlier = TRUE, 
#                addStat = 'median', 
#                dashboard = TRUE, 
#                titvRaw = FALSE)


# 在oncoprint中加入自定义CNV图景
interest_cytoband <- c('TERT','F12','CDK4','TERF2','CETN2','GPS1','BRD4','CTBP1','ANG','CCNE1','ZNRF3','DCTD','RERE','LSAMP','CDKN2A','RB1','PTPRD','PRKAR1A','CTSD','AGXT','ACLY','ADRA1D','RBFOX1','CCR6','DLG2','RAB11FIP4')
gistic <- read.table(file = './02genomic_analysis/source_data/TCGA-ACC.cnv.gistic.tsv',sep = '\t',header = T,check.names = F) %>% 
  .[.$Symbol%in%interest_cytoband,] %>% 
  .[match(interest_cytoband,.$Symbol),] %>% 
  mutate(Symbol=c('5p15.33','5q35.3','12q14.1','16q22.1','Xq28','17q25.3','19p13.12','4p16.3','14q11.2','19q12','22q12.1','4q35.1','1p36.23','3q13.31','9p21.3','13q14.2','9p23','17q24.2','11p15.5','2q37.3','17q21.31','20p12.1','16p13.3','6q26','11q14.1','17q11.2')) %>% 
  gather(.,key = 'PATIENT',value ='gistic_score',-1) %>% 
  mutate(gistic_score=ifelse(.$Symbol=='11p15.5' & .$gistic_score>0,0,.$gistic_score)) %>% 
  mutate(gistic_score=ifelse(.$Symbol=='14q11.2' & .$gistic_score<0,0,.$gistic_score)) %>% 
  mutate(CN=ifelse(.$gistic_score>0,'Amp',ifelse(.$gistic_score<0,'Del',NA))) %>% 
  mutate(gistic_score=NULL) %>% 
  na.omit() %>% 
  plyr::rename(c('PATIENT'='Tumor_Sample_Barcode')) %>% 
  .[.$Tumor_Sample_Barcode%in%clinic$PATIENT_ID,]

mut.plus.cnv = read.maf(maf = maf,
                        cnTable = gistic,
                        clinicalData = clinicalData,
                        verbose = F)
features <- c('TP53', 'ZNFR3', 'CTNNB1', 'PRKAR1A', 'CCNE1','TERF2','MEN1','RPL22','5p15.33','5q35.3','12q14.1','16q22.1','Xq28','17q25.3','19p13.12','4p16.3','14q11.2','19q12','22q12.1','4q35.1','1p36.23','3q13.31','9p21.3','13q14.2','9p23','17q24.2','11p15.5','2q37.3','17q21.31','20p12.1','16p13.3','6q26','11q14.1','17q11.2')#选择特定的基因或cytoband
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired') #自定义颜色
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)

oncoplot(maf = mut.plus.cnv, 
         genes = features, #interest genes
         removeNonMutated = F,
         colors = vc_cols, 
         sortByMutation = T,
         fontSize = 0.5,
         writeMatrix = T,  #export as ‘oncomatrix.txt’ file
         clinicalFeatures = c('Without_combined',"Likely_combined","PC_combined","Putative_combined",'AGE_new'),
         logColBar=T,
         sortByAnnotation = T)

# comparison of genes between ms1/2
data_type <- c('Without_combined',"Likely_combined","PC_combined","Putative_combined")
for (i in 1:4) {
  ms1_patient <- clinic$PATIENT_ID[clinic[,data_type[i]]=='MS1']
  ms2_patient <- clinic$PATIENT_ID[clinic[,data_type[i]]=='MS2']
  mut1 <- maf[maf$Tumor_Sample_Barcode%in%ms1_patient,] 
  mut2 <- maf[maf$Tumor_Sample_Barcode%in%ms2_patient,] 
  gistic1 <- gistic[gistic$Tumor_Sample_Barcode%in%ms1_patient,]
  gistic2 <- gistic[gistic$Tumor_Sample_Barcode%in%ms2_patient,]
  mut.plus.cnv1 = read.maf(maf = mut1,
                           cnTable = gistic1,
                           clinicalData = clinicalData,
                           verbose = F)
  mut.plus.cnv2 = read.maf(maf = mut2,
                           cnTable = gistic2,
                           clinicalData = clinicalData,
                           verbose = F)
  pt.vs.rt <- mafCompare(m1 = mut.plus.cnv1, m2 = mut.plus.cnv2, m1Name = 'MS1', m2Name = 'MS2',minMut = 0) #比较最少Mut个数为0的基因
  # forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.05, color = c('royalblue', 'maroon'), geneFontSize = 0.8)
  diff_mut_cnv_genes <- pt.vs.rt$results %>% 
    .[.$Hugo_Symbol%in%features,] #挑选热点基因
  write.table(diff_mut_cnv_genes,file = paste0('./02genomic_analysis/output_Diff_genes_of_mut_and_cnv(customed)_between_MS1_MS2_',data_type[i],'.tsv'),sep = '\t',row.names = F)
}

# functional enrichment of mut and cnv genes
pdf('output_OncogenicPathways_CR')
OncogenicPathways(maf = mut.plus.cnv1) 
PlotOncogenicPathways(maf = mut.plus.cnv1, pathways = "TP53")


# Visualization of mut and cnv(customed) by complexheatmap ----------------
rm(list = ls())
{
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(circlize)
}
load('step00_output.Rdata')

mat<-read.table(file="./02genomic_analysis/onco_matrix_heatmap.txt",check.names = F,
                header = T, row.names = 1,stringsAsFactors = FALSE,sep = "\t")
mat[mat==0]<-""
mat[1:3, 1:3]

vc_cols = RColorBrewer::brewer.pal(n = 10, name = 'Paired') #自定义颜色
vc_cols[9:10] <- c('#104E8B','#FF3030')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del',
  'Del',
  'Amp'
)
vc_cols

# col_fun <- colorRamp2(
#   c(0, 5, 10),
#   c("#ff7f00", "white", "#1f78b4")
# )

col <- vc_cols
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
  },
  # bug red
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["Missense_Mutation"], col = NA))
  },
  # small green
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["Splice_Site"], col = NA))
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
  },
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))
  },
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["In_Frame_Del"], col = NA))
  },
  Multi_Hit = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["Multi_Hit"], col = NA))
  },
  Amp = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Amp"], col = NA))
  },
  Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"),h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Del"], col = NA))
  }
)

clinic <- clinic[order(clinic$Without_combined),] #clinic按照kmeans排序
mat <- mat[,clinic$PATIENT_ID]#跟clinic相同顺序

row_group <- c(rep('Mutation',5),rep('CNV',26)) #行分组信息
names(row_group) <- rownames(mat)

# colomn_group <- c(rep('MS1',55),rep('MS2',22)) #列分组信息
# names(colomn_group) <- colnames(mat)

# column_order <- match(colnames(mat),clinic$PATIENT_ID) #列顺序
# mat <- mat[,column_order]


#顶部注释
top_anno = HeatmapAnnotation(NR = clinic$Without_combined,
                             LR=clinic$Likely_combined,
                             CR=clinic$PC_combined,
                             PR=clinic$Putative_combined,
                             Age=clinic$AGE,
                             Sex=clinic$SEX,
                             gap = unit(1, "mm"),
                             col = list(NR = c("MS1" = "#0072B5CC", "MS2" = "#E18727CC"),
                                        LR = c("MS1" = "#0072B5CC", "MS2" = "#E18727CC"),
                                        CR = c("MS1" = "#0072B5CC", "MS2" = "#E18727CC"),
                                        PR = c("MS1" = "#0072B5CC", "MS2" = "#E18727CC"),
                                        Sex=c('Male'='#8B4C39','Female'='#FF82AB'),
                                        Age=colorRamp2(c(0, 80), c("white", "red")))
)

#bottom annotation
bottom_anno <- HeatmapAnnotation(column_barplot = anno_oncoprint_barplot(height = unit(1,'cm'))) 

#left annotation
# data1 <- read.table(file = 'Diff_analysis_of_mut_and_cnv_between_MS1_MS2.tsv',header = T,sep = '\t',check.names = F)
# data1 <- data1[match(rownames(mat),data1$Hugo_Symbol),]
# bardata=data.frame(OR=data1[,'or'],row.names = rownames(mat))#取OR数据做条形图
# col_fun_row=colorRamp2(c(min(data1$pval),0.05,max(data1$pval)),c('#EC211C','#316CD7','blue'))#根据P值范围生成P值映射颜色范围
# # col_fun_row=colorRamp2(c(0,0.5,max(data1$pval)),c('#EC211C','#316CD7','blue'))#根据P值范围生成P值映射颜色范围
# barcol=col_fun_row(data1$pval)#生成所需基因的P值对应的颜色
# names(barcol)=rownames(bardata)
# pvalue_col <- ifelse(data1$pval<0.05,'significant','ns')
# names(pvalue_col) <- data1$Hugo_Symbol
# 
# left_anno = rowAnnotation('Odd Ratio'=row_anno_barplot(bardata, 
#                                width = unit(1.5, "cm"),#条形图高度
#                                bar_width = 0.8, #条形图条的宽度
#                                outline = FALSE, 
#                                gp = gpar(fill = barcol,lwd=0.1)),##条形图颜色,条形图边框宽度
#             'Pvalue'=pvalue_col,
#             # annotation_height=unit(c(1.5,0.1),'cm'),
#             annotation_width = unit(c(1.5,0.2), c("cm", "cm")),
#             width = unit(1.7, "cm"),
#             show_annotation_name = c(Pvalue = FALSE),
#             # show_legend = c("Pvalue" = FALSE),
#             col = list(Pvalue = c("significant" = "red", "ns" = 'gray95')))

# right annotation
data_type_1 <- c('Without_combined',"Likely_combined","PC_combined","Putative_combined")
data_type_2 <- c('NR','LR','CR','PR')
text=data.frame(row.names = rownames(mat))
for (i in 1:4) {
  data <- read.table(file = paste0('./02genomic_analysis/output_Diff_genes_of_mut_and_cnv(customed)_between_MS1_MS2_',data_type_1[i],'.tsv'),sep = '\t',header = T) %>%
    .[match(rownames(mat),.$Hugo_Symbol),]
  data <- data[match(rownames(mat),data$Hugo_Symbol),]
  text[,i]=round(data[,'pval'],2)
}
text <- ifelse(text<0.05,'*','')

right_anno <- rowAnnotation( 'NR'=anno_text(text[,'V1'],which = 'row',show_name = TRUE,
                                                  location = 0.5, just = "center",
                                                  gp = gpar(fill = NA, col = "red",fontsize=25, border = NA),
                                                  width = max_text_width(month.name)*0.5),
                             'LR'=anno_text(text[,'V2'],which = 'row',show_name = TRUE,
                                            location = 0.5, just = "center",
                                            gp = gpar(fill = NA, col = "red", fontsize=25,border = NA),
                                            width = max_text_width(month.name)*0.5),
                             'CR'=anno_text(text[,'V3'],which = 'row',show_name = TRUE,
                                            location = 0.5, just = "center",
                                            gp = gpar(fill = NA, col = "red",fontsize=25, border = NA),
                                            width = max_text_width(month.name)*0.5),
                             'PR'=anno_text(text[,'V4'],which = 'row',show_name = TRUE,
                                            location = 0.5, just = "center",
                                            gp = gpar(fill = NA, col = "red",fontsize=25, border = NA),
                                            width = max_text_width(month.name)*0.5),
                             annotation_name_rot = 0)

# plot heatmap
pdf('./02genomic_analysis/output_heatmap_genomics_ms1_ms2.pdf',width = 10)
oncoPrint(mat,
          alter_fun = alter_fun, 
          col = vc_cols, 
          column_title='77 ACC Samples',
          row_split = row_group,
          # column_split = colomn_group,
          row_gap = unit(1, "mm"),
          show_pct = F,
          row_names_side = "left",
          row_order = 1:nrow(mat),
          column_order = 1:ncol(mat),
          # column_order = column_order,
          top_annotation = top_anno,
          bottom_annotation =bottom_anno,
          # right_annotation = right_anno,
          show_heatmap_legend=T,
          remove_empty_columns = F, 
          remove_empty_rows = F)
dev.off()


# # Visualization of mut and cnv(customed) significantly by complexheatmap --------
mat1 <- mat[c('14q11.2','22q12.1','9p21.3','CTNNB1','TP53'),]

vc_cols = RColorBrewer::brewer.pal(n = 10, name = 'Paired') #自定义颜色
vc_cols[9:10] <- c('#104E8B','#FF3030')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del',
  'Del',
  'Amp'
)
vc_cols
col <- vc_cols
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
  },
  # bug red
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["Missense_Mutation"], col = NA))
  },
  # small green
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["Splice_Site"], col = NA))
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
  },
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))
  },
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["In_Frame_Del"], col = NA))
  },
  Multi_Hit = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["Multi_Hit"], col = NA))
  },
  Amp = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Amp"], col = NA))
  },
  Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"),h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Del"], col = NA))
  }
)


clinic <- clinic[order(clinic$Without_combined),] #clinic按照NR排序
mat1 <- mat1[,clinic$PATIENT_ID]#跟clinic相同顺序

row_group <- c(rep('CNV',3),rep('Mutation',2)) #行分组信息
names(row_group) <- rownames(mat1)

# colomn_group <- c(rep('MS1',55),rep('MS2',22)) #列分组信息
# names(colomn_group) <- colnames(mat)

# column_order <- match(colnames(mat),clinic$PATIENT_ID) #列顺序
# mat <- mat[,column_order]


#顶部注释
top_anno = HeatmapAnnotation(NR = clinic$Without_combined,
                             LR=clinic$Likely_combined,
                             CR=clinic$PC_combined,
                             PR=clinic$Putative_combined,
                             Age=clinic$AGE,
                             Sex=clinic$SEX,
                             gap = unit(1, "mm"),
                             col = list(NR = c("MS1" = "#0072B5CC", "MS2" = "#E18727CC"),
                                        LR = c("MS1" = "#0072B5CC", "MS2" = "#E18727CC"),
                                        CR = c("MS1" = "#0072B5CC", "MS2" = "#E18727CC"),
                                        PR = c("MS1" = "#0072B5CC", "MS2" = "#E18727CC"),
                                        Sex=c('Male'='#8B4C39','Female'='#FF82AB'),
                                        Age=colorRamp2(c(0, 80), c("white", "red")))
)

#bottom annotation
bottom_anno <- HeatmapAnnotation(column_barplot = anno_oncoprint_barplot(height = unit(1,'cm'))) 


## plot
pdf('./02genomic_analysis/output_heatmap_genomics_signifcant.pdf',width = 15,height = 3)
oncoPrint(mat1,
          alter_fun = alter_fun,
          row_order = 1:nrow(mat1),
          column_order = 1:ncol(mat1),
          show_pct = F,
          row_split = row_group,
          row_names_side = "left",
          top_annotation = top_anno,
          # right_annotation = right_anno,
          # column_split = colomn_group,
          column_gap = unit(3, "mm"),
          col = vc_cols)
dev.off()



# Landscape of mut and cnv ------------------------------------------------
# rm(list = ls())
# library(TCGAbiolinks)
# 
# # download masked copy number segment
# query <- GDCquery(project = "TCGA-ACC",
#                   data.category="Copy Number Variation",
#                   data.type="Masked Copy Number Segment",
#                   access = "open")
# GDCdownload(query,method = "api",files.per.chunk = 50)
# expdat <- GDCprepare(query = query)
# length(dir('Masked_Copy_Number_Segment') ) ## check 'Masked_Copy_Number_Segment' files件
# dir.create("rawdata_cnv_all") ## create a folder setting '.txt' for Masked_Copy_Number_Segment files
# 
# #用lapply将rawdata里面txt文件整理到一个文件夹里
# raw_txt <- dir("Masked_Copy_Number_Segment/")
# lapply(raw_txt, function(x){
#   mydir <- paste0("./Masked_Copy_Number_Segment/",x)
#   file <- list.files(mydir,pattern = "nocnv_grch38")
#   myfile <- paste0("./Masked_Copy_Number_Segment/",x,"/",file)
#   file.copy(myfile,"rawdata_cnv_all")  
# })
# 
# library(data.table)
# cnv_file <- dir("rawdata_cnv_all")
# cnv_file[1]
# file <- paste0("./rawdata_cnv_all/",cnv_file) ## 文件所在位置
# head(file )
# 
# # 批量读取txt文件
# file_list <- list()
# fred_cnv <- function(x){
#   file_list[[x]] <- data.table::fread(file = x,data.table = F)
#   file_list[[x]]   
# }
# cnv_df <- lapply(file, fred_cnv)
# cnv_df <- do.call(rbind,cnv_df)
# head(cnv_df)
# 
# # 读取metadata里面的注释信息
# metadata <- jsonlite::fromJSON("metadata.cart.2022-05-07.json")
# library(dplyr)
# metadata_id <- metadata %>% 
#   dplyr::select(c(file_name,associated_entities)) 
# meta_df <- do.call(rbind,metadata$associated_entities)
# head(metadata_id,2)
# 
# # 匹配样本
# cnv_df$Sample <- meta_df$entity_submitter_id[match(cnv_df$GDC_Aliquot,meta_df$entity_id)]
# length(unique(cnv_df$GDC_Aliquot))
# length(unique(cnv_df$Sample))
# head(cnv_df)
# cnv_df$Sample <- substring(cnv_df$Sample,1,16) #改名
# 
# # Gistic只需要maskedCNS的六列
# cnv_df <- cnv_df[,c('Sample','Chromosome','Start','End','Num_Probes','Segment_Mean')]
# head(cnv_df)
# dim(cnv_df)
# cnv_df <- cnv_df[grep("01A$",cnv_df$Sample),] # 只挑选肿瘤样本
# dim(cnv_df)
# head(cnv_df)
# cnv_df$Sample <- substring(cnv_df$Sample,1,12) 
# ms1_patient <- clinic$PATIENT_ID[clinic$kmeans=='MS1']
# ms2_patient <- clinic$PATIENT_ID[clinic$kmeans=='MS2']
# cnv_df_1 <- cnv_df[cnv_df$Sample%in%ms1_patient,] #MS1病人的cnv数据
# cnv_df_2 <- cnv_df[cnv_df$Sample%in%ms2_patient,] #MS2病人的cnv数据
# 
# # 整理好的数据
# write.table(cnv_df_1,"Masked Copy Number Segment(Tumor_MS1).txt",sep="\t",
#             quote = F,col.names = F,row.names = F)
# write.table(cnv_df_2,"Masked Copy Number Segment(Tumor_MS2).txt",sep="\t",
#             quote = F,col.names = F,row.names = F)
# write.table(cnv_df,"Masked Copy Number Segment(Tumor).txt",sep="\t",
#             quote = F,col.names = F,row.names = F)
# 
# # marker file数据下载和处理
# #GDC Reference Files---选择最新版本的“SNP6 GRCh38 Remapped Probeset File for Copy Number Variation Analysis”文件snp6.na35.remap.hg38.subset.txt.gz
# hg38_marker_file <- read.delim("./data/genomic_analysis/snp6.na35.remap.hg38.subset.txt.gz")
# head(hg38_marker_file )
# #注意“If you are using Masked Copy Number Segment for GISTIC analysis, please only keep probesets with freqcnv =FALSE”
# hg_marker_file <- hg38_marker_file[hg38_marker_file$freqcnv=="FALSE",]
# hg_marker_file <- hg_marker_file [,c(1,2,3)]
# head(hg_marker_file )
# write.table(hg_marker_file,file = "hg_marker_file.txt",sep = "\t",col.names = T,row.names = F)
# 
# ## 使用GISTIC2.0来识别体细胞拷贝数改变（SCNA），然后找到这些拷贝数显著变化的多基因区域。
# # https://cloud.genepattern.org/gp/pages/index.jsf 网址 username:LIYUQING password:
# # https://www.jianshu.com/p/755f8357c898 教程
# 
# ## 使用maftools读取和汇总gistic输出文件
# gistic.1 = readGistic(gisticAllLesionsFile = './data/genomic_analysis/Gistic 2.0 results/all_lesions_ms1.conf_99.txt',
#                       gisticAmpGenesFile = './data/genomic_analysis/Gistic 2.0 results/amp_genes_ms1.conf_99.txt', 
#                       gisticDelGenesFile = './data/genomic_analysis/Gistic 2.0 results/del_genes_ms1.conf_99.txt', 
#                       gisticScoresFile = './data/genomic_analysis/Gistic 2.0 results/scores_ms1.gistic', 
#                       isTCGA = T)
# gistic.2 = readGistic(gisticAllLesionsFile = './data/genomic_analysis/Gistic 2.0 results/all_lesions_ms2.conf_99.txt',
#                       gisticAmpGenesFile = './data/genomic_analysis/Gistic 2.0 results/amp_genes_ms2.conf_99.txt', 
#                       gisticDelGenesFile = './data/genomic_analysis/Gistic 2.0 results/del_genes_ms2.conf_99.txt', 
#                       gisticScoresFile = './data/genomic_analysis/Gistic 2.0 results/scores_ms2.gistic', 
#                       isTCGA = T)
# 
# #genome plot 部分突出了明显的扩增和删失区域
# gisticChromPlot(gistic = gistic.2, ref.build = 'hg38') ## 参考基因组选hg38
# 
# #Bubble plot 绘制显著改变的cytobands,以及改变的样本数和它所包含的基因数的函数。每一个气泡的大小是根据-log10转化q值。
# gisticBubblePlot(gistic = gistic.2)
# 
# #ms1,ms2分别合并mut和cnv
# load('step00_output.Rdata')
# maf <- read.csv(file = './data/genomic_analysis/ACC_mutect2.csv',header = T,check.names = F) %>% 
#   mutate(Tumor_Sample_Barcode=substr(.$Tumor_Sample_Barcode,1,12))
# 
# ms1_patient <- clinic$PATIENT_ID[clinic$kmeans=='MS1']
# ms2_patient <- clinic$PATIENT_ID[clinic$kmeans=='MS2']
# maf1 <- maf[maf$Tumor_Sample_Barcode%in%ms1_patient,] #MS1病人34
# maf2 <- maf[maf$Tumor_Sample_Barcode%in%ms2_patient,] #MS2病人43
# 
# mut.plus.cnv.1 <- read.maf(maf=maf1, 
#                               gisticAllLesionsFile = './data/genomic_analysis/Gistic 2.0 results/all_lesions_ms1.conf_99.txt',
#                               gisticAmpGenesFile = './data/genomic_analysis/Gistic 2.0 results/amp_genes_ms1.conf_99.txt', 
#                               gisticDelGenesFile = './data/genomic_analysis/Gistic 2.0 results/del_genes_ms1.conf_99.txt', 
#                               gisticScoresFile = './data/genomic_analysis/Gistic 2.0 results/scores_ms1.gistic')
# mut.plus.cnv.2 <- read.maf(maf=maf2, 
#                               gisticAllLesionsFile = './data/genomic_analysis/Gistic 2.0 results/all_lesions_ms2.conf_99.txt',
#                               gisticAmpGenesFile = './data/genomic_analysis/Gistic 2.0 results/amp_genes_ms2.conf_99.txt', 
#                               gisticDelGenesFile = './data/genomic_analysis/Gistic 2.0 results/del_genes_ms2.conf_99.txt', 
#                               gisticScoresFile = './data/genomic_analysis/Gistic 2.0 results/scores_ms2.gistic')
# # plot enriched pathways
# pdf('/data/genomic_analysis/oncogenic_pathways_ms1.pdf')
# OncogenicPathways(maf = mut.plus.cnv.1)
# dev.off()
# pdf('/data/genomic_analysis/oncogenic_pathways_ms2.pdf')
# OncogenicPathways(maf = mut.plus.cnv.2)
# dev.off()
# 
# #oncoplot
# gisticOncoPlot(gistic = gistic.2, top = 10,removeNonAltered = F)
# oncoplot(maf=mut.plus.cnv.2, borderCol=NULL, top=10)
# pt.vs.rt <- mafCompare(m1 = mut.plus.cnv.1, m2 = mut.plus.cnv.2, m1Name = 'MS1', m2Name = 'MS2',minMut = 5) #比较最少Mut个数为0的基因
# forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.05, color = c('royalblue', 'maroon'), geneFontSize = 0.2)
# diff_mut_cnv_genes <- pt.vs.rt$results
# write.table(diff_mut_cnv_genes,file = './data/genomic_analysis/Diff_genes_of_mut_and_cnv_between_MS1_MS2.tsv',sep = '\t',row.names = F)


# TMB and FGA Comparison between MS1 & MS2 -------------------------------------------
# install.packages("remotes")
# remotes::install_github("JanCoUnchained/ggunchained")
library(ggunchained)

# TMB
tmb = tmb(maf = read.maf(maf = maf), captureSize = 50, logScale = FALSE) %>% 
  as.data.frame() %>% 
  merge(.,clinic,by.x="Tumor_Sample_Barcode",by.y="PATIENT_ID")
outlier_values <- boxplot.stats(tmb$total_perMB)$out   #检测离群值
tmb <- tmb[!(tmb$total_perMB %in% outlier_values),]
tmb.x <- select(tmb,c(total_perMB,Without_combined,Likely_combined,PC_combined,Putative_combined)) %>%
  gather(.,data_type,MS,-total_perMB) %>% 
  dplyr::mutate(data_type=factor(.$data_type,levels = c('Without_combined','Likely_combined','PC_combined','Putative_combined')))

pdf('./02genomic_analysis/output_TMB_ms1_ms2.pdf')
ggplot(tmb.x, aes(x=data_type, y=total_perMB,fill=MS)) +
  #"trim"为TRUE(默认值),将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  # geom_jitter(shape=16, color="grey",size=2.0,position=position_jitter(0.2))+
  geom_split_violin(trim = FALSE,color = NA)+
  stat_summary(fun = mean,
               fun.min = function(x){quantile(x)[2]},
               fun.max = function(x){quantile(x)[4]},
               geom = "pointrange",
               #geom = 'errorbar',
               size=0.5,
               position = position_dodge(width = 0.2))+
  theme_bw()+
  thm+
  scale_fill_manual(values = c(MS1="#0072B5CC",MS2="#E18727CC"))+
  guides(fill=guide_legend(direction = "horizontal",nrow =1))+
  labs(x="", y = "Tumor Muation Burden(Total_perMB)", fill = "")+  # fill为修改图例标题
  stat_compare_means(aes(x = data_type, y = total_perMB,group=MS),method = "wilcox.test",label = "p.signif",size=8)
dev.off()

# Fraction Genome Altered(FGA)
FGA <- select(clinic,c(Fraction.Genome.Altered,Without_combined,Likely_combined,PC_combined,Putative_combined)) %>%
  gather(.,data_type,MS,-Fraction.Genome.Altered) %>% 
  dplyr::mutate(data_type=factor(.$data_type,levels = c('Without_combined','Likely_combined','PC_combined','Putative_combined')))
pdf('./02genomic_analysis/output_fraction_geome_altered_ms1_ms2.pdf')
ggplot(FGA, aes(x=data_type, y=Fraction.Genome.Altered,fill=MS)) +
  geom_split_violin(trim = FALSE,color = NA)+
  stat_summary(fun = mean,
               fun.min = function(x){quantile(x)[2]},
               fun.max = function(x){quantile(x)[4]},
               geom = "pointrange",
               # geom = 'errorbar',
               size=0.5,
               position = position_dodge(width = 0.2))+
  theme_bw()+
  thm+
  scale_fill_manual(values = c(MS1="#0072B5CC",MS2="#E18727CC"))+
  guides(fill=guide_legend(direction = "horizontal",nrow =1))+
  labs(x="", y = "Fraction Genome Altered", fill = "")+  # fill为修改图例标题
  stat_compare_means(aes(x = data_type, y = Fraction.Genome.Altered,group=MS),method = "wilcox.test",label = "p.signif",size=8)
dev.off()

