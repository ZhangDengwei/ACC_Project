getwd()
rm(list = ls())

##ACC转录组表达数据
{library(stringr)
  library(org.Hs.eg.db)
  library(clusterProfiler)}
dat1 <- read.table(file = 'TCGA-ACC.htseq_counts.tsv',header = T,row.names = 1,check.names = F,sep = '\t') #from xena
dat1 <- 2^(dat1)-1
rownames(dat1) <- str_sub(rownames(dat1),1,str_locate(rownames(dat1),'\\.')[1]-1)
dat2 <- read.table(file = 'TCGA-ACC.htseq_fpkm.tsv',header = T,row.names = 1,check.names = F,sep = '\t')  #from xena
dat2 <- 2^(dat2)-1
rownames(dat2) <- str_sub(rownames(dat2),1,str_locate(rownames(dat2),'\\.')[1]-1)

list <- bitr(unique(rownames(dat1)), 
             fromType = "ENSEMBL",
             toType = c( "SYMBOL"),
             OrgDb = org.Hs.eg.db)
dim(list)

#count矩阵
dat1$Symbol <- list[match(rownames(dat1),list[,'ENSEMBL']),'SYMBOL']
dat1 <- na.omit(dat1)
dat1=aggregate(.~Symbol,mean,data=dat1) #相同的行合并并取均值
rownames(dat1) <- dat1$Symbol
dat1 <- dat1[,-1] #去除symbol列
colnames(dat1) <- substr(colnames(dat1),1,12)
write.table(dat1,file = 'TCGA-ACC.htseq_counts.tsv',sep = '\t',row.names = T,col.names = NA)

#fpkm矩阵
dat2$Symbol <- list[match(rownames(dat2),list[,'ENSEMBL']),'SYMBOL']
dat2 <- na.omit(dat2)
dat2=aggregate(.~Symbol,mean,data=dat2) #相同的行合并并取均值
rownames(dat2) <- dat2$Symbol
dat2 <- dat2[,-1] #去除symbol列
colnames(dat2) <- substr(colnames(dat2),1,12)
write.table(dat2,file = 'TCGA-ACC.htseq_fpkm.tsv',sep = '\t',row.names = T,col.names = NA)

#tpm矩阵
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
dat3 <- apply(dat2,2,fpkmToTpm)
colSums(dat3)
write.table(dat3,file = 'TCGA-ACC.htseq_tpm.tsv',sep = '\t',row.names = T,col.names = NA)


##ACC CNV数据
acc_cnv <- read.table(file = 'TCGA-ACC.gistic.tsv',header = T,row.names = 1,check.names = F,sep = '\t') #from xena
rownames(acc_cnv) <- str_sub(rownames(acc_cnv),1,str_locate(rownames(acc_cnv),'\\.')[1]-1)
list <- bitr(unique(rownames(acc_cnv)), 
             fromType = "ENSEMBL",
             toType = c( "SYMBOL"),
             OrgDb = org.Hs.eg.db)
acc_cnv$Symbol <- list[match(rownames(acc_cnv),list[,'ENSEMBL']),'SYMBOL']
acc_cnv <- na.omit(acc_cnv)
acc_cnv <- acc_cnv[!duplicated(acc_cnv$Symbol),]
rownames(acc_cnv) <- acc_cnv$Symbol
acc_cnv <- acc_cnv[,-91] #去除symbol列
colnames(acc_cnv) <- substr(colnames(acc_cnv),1,12)
acc_cnv <- tibble::rownames_to_column(acc_cnv,var = 'Symbol')
write.table(acc_cnv,file = 'TCGA-ACC.cnv.gistic.tsv',sep = '\t',row.names = T,col.names = NA)


##ACC 临床数据准备
{library(tidyr)
library(dplyr)
library(magrittr)}
#dir()：获取文件夹下所有文件;pattern：匹配正则表达式;full.names：是否包含文件绝对路径
getwd()
file_names <- dir('./00data_preparation/source_data/', pattern = "*.txt", full.names = T)
df <- read.table(file_names[1],header = T,sep = '\t',check.names = F)
for (i in 2:length(file_names)) {
  df <- merge(df,read.table(file_names[i],header = T,sep = '\t',check.names = F),by='PATIENT_ID',all=T)}
clinic <- df %>%
  drop_na(id) %>%
  transform(.,Without_combined=as.factor(paste0('MS',.$Without_combined))) %>%
  transform(.,Likely_combined=as.factor(paste0('MS',.$Likely_combined))) %>%
  transform(.,PC_combined=as.factor(paste0('MS',.$PC_combined))) %>%
  transform(.,Putative_combined=as.factor(paste0('MS',.$Putative_combined)))  %>% 
  transform(.,Strigent_combined=as.factor(paste0('MS',.$Strigent_combined)))
clinic[clinic=='[Not Available]'] <- NA


# save(clinic,file = 'step00_output.Rdata')
# write.table(clinic,file = 'ACC_clinic.tsv',sep = '\t',row.names = T,col.names = NA)

###卡方检验
# load('step00_output.Rdata')
# library(stringr)
# clinic$HISTORY_ADRENAL_HORMONE_EXCESS_new <- ifelse(str_detect(clinic$HISTORY_ADRENAL_HORMONE_EXCESS,'Cortisol'),'Cortisol',
#                                                    ifelse(str_detect(clinic$HISTORY_ADRENAL_HORMONE_EXCESS,'None'),'None','Others'))
# head(clinic$HISTORY_ADRENAL_HORMONE_EXCESS_new)
# name <- 'HISTORY_ADRENAL_HORMONE_EXCESS_new'
# table(clinic[clinic$kmeans=='MS1',name]);table(clinic[clinic$kmeans=='MS2',name])
# x <- as.matrix(cbind(as.array(table(clinic[clinic$kmeans=='MS1',name])),as.array(table(clinic[clinic$kmeans=='MS2',name]))))
# # x <- matrix(c(15,14,5,18,12,9),ncol = 2)
# x
# if (sum(x<5)>0){
#   a <- fisher.test(x)
#   print(a)
# }else{
#   a <- chisq.test(x)
#   print(a)
# }
# 


