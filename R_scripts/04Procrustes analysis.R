rm(list = ls())
getwd()
{library(vegan)
library(coda.base)
library(ggplot2)}
load('step00_output.Rdata')
filter_genes <- function(genes, qt){
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

## microbiome matrix
acc_micro_matrix <- read.table(file = './data/ACC_77_filter_144_genus_TMM_normalization.tsv',header = T,row.names = 1,check.names = F,sep = '\t')
colnames(acc_micro_matrix) <- clinic$PATIENT_ID[match(colnames(acc_micro_matrix),clinic$names)]

microbes <- acc_micro_matrix
microbes <- t(microbes)
dim(microbes)
#[1]  77 144
ms2_patient <- clinic$PATIENT_ID[clinic$kmeans=='MS2']
microbes.ms2 <- microbes[ms2_patient,]
dim(microbes.ms2)
#[1]  43 144

## mrna exp matrix
acc_exp <- read.table(file = './data/TCGA-ACC.htseq_tpm.tsv',header = T,row.names = 1,check.names = F,sep = '\t')
genes <- t(acc_exp[,rownames(microbes)])
FUN =  function(xxx) {
  (sum(xxx==0)/length(xxx))<0.5 }
genes <- genes[,apply(genes, 2, FUN)]
genes <- filter_genes(genes,2) #过滤低表达基因
dim(genes)
#[1]    77 14267
genes.ms2 <- genes[ms2_patient,]
dim(genes.ms2)
#[1]    43 14267


# MS1 & MS2 Patients ------------------------------------------------------

## Procrustes analysis

load('step06_output.Rdata')
micro.dist <- vegan::vegdist(microbes,method = 'robust.aitchison')
gene.dist <- vegan::vegdist(genes,method = 'bray')
mds.m <- monoMDS(micro.dist)
mds.g <- monoMDS(gene.dist)

pro.m.g <- procrustes(mds.g,mds.m,symmetric = TRUE)
summary(pro.m.g)

# plot(pro.m.g,kind=2)
# residuals(pro.m.g)

# M2 calculating
set.seed(123)
pro.m.g_t <- protest(mds.g,mds.m,permutations = 9999)
pro.m.g_t

# result
pro.m.g_t$ss #偏差平方和（M2统计量）
pro.m.g_t$signif #对应p值结果

# 获得x和y轴的坐标及旋转过的坐标
Pro_Y <- cbind(data.frame(pro.m.g$Yrot),data.frame(pro.m.g$X))
Pro_X <- data.frame(pro.m.g$rotation)

# plot
pdf('procrusters_analysis_all_patients.pdf')
ggplot(data = Pro_Y)+
  geom_segment(aes(x=X1,y=X2,xend=(X1+MDS1)/2,yend=(X2+MDS2)/2),
               arrow = arrow(length = unit(0,'cm')),
               color='grey',size=0.5)+
  geom_segment(aes(x=(X1+MDS1)/2,y=(X2+MDS2)/2,
                   xend=MDS1,yend=MDS2),
               arrow = arrow(length = unit(0.2,'cm')),
               color='grey',size=0.5)+
  geom_point(aes(X1,X2,shape='17'),color='#9BBB59',size=4)+  #shape放在aes里面要加引号‘’
  geom_point(aes(MDS1,MDS2,shape='1'),color='#9BBB59',size=4)+
  scale_shape_manual(name = "",
                     values = c('17' = 17,'1' = 1),
                     labels = c('Microbiota abundance(Bray-Curtis)', 'Host gene expression(Bray-Curtis)'))+
  annotate('text',label=paste0('        Procrusters analysis\n
           M2=',round(pro.m.g_t$ss,3), ', p-value=',round(pro.m.g_t$signif,3)),
           x=-0.2,y=0.1,size=4)+
  labs(x='Dimension 1',y='Dimension 2')+
  theme_bw()+
  theme(legend.position = 'bottom')
dev.off()

## Mantel
micro.dist <- vegan::vegdist(microbes,method = 'bray')
gene.dist <- vegan::vegdist(genes,method = 'robust.aitchison')

set.seed(123)
man.m.g <- mantel(micro.dist,gene.dist,method = 'spearman', permutations = 9999, na.rm = TRUE)
man.m.g


# MS2 Patients ------------------------------------------------------------


## Procrustes analysis
micro.dist.2 <- vegan::vegdist(microbes.ms2,method = 'robust.aitchison')
gene.dist.2 <- vegan::vegdist(genes.ms2,method = 'bray')
mds.m.2 <- monoMDS(micro.dist.2)
mds.g.2 <- monoMDS(gene.dist.2)

pro.m.g.2 <- procrustes(mds.g.2,mds.m.2,symmetric = TRUE)
summary(pro.m.g.2)

#M2统计量的显著性检验
set.seed(123)
pro.m.g_t.2 <- protest(mds.g.2,mds.m.2,permutations = 9999)
pro.m.g_t.2

#提取结果
pro.m.g_t.2$ss #偏差平方和（M2统计量）
pro.m.g_t.2$signif #对应p值结果

#获得x和y轴的坐标及旋转过的坐标
Pro_Y.2 <- cbind(data.frame(pro.m.g.2$Yrot),data.frame(pro.m.g.2$X))
Pro_X.2 <- data.frame(pro.m.g.2$rotation)

# plot
pdf('procrusters_analysis_ms2_patients.pdf')
ggplot(data = Pro_Y.2)+
  geom_segment(aes(x=X1,y=X2,xend=(X1+MDS1)/2,yend=(X2+MDS2)/2),
               arrow = arrow(length = unit(0,'cm')),
               color='grey',size=0.5)+
  geom_segment(aes(x=(X1+MDS1)/2,y=(X2+MDS2)/2,
                   xend=MDS1,yend=MDS2),
               arrow = arrow(length = unit(0.2,'cm')),
               color='grey',size=0.5)+
  geom_point(aes(X1,X2,shape='17'),color='#9BBB59',size=4)+  #shape放在aes里面要加引号‘’
  geom_point(aes(MDS1,MDS2,shape='1'),color='#9BBB59',size=4)+
  scale_shape_manual(name = "",
                     values = c('17' = 17,'1' = 1),
                     labels = c('Microbiota abundance(Bray-Curtis)', 'Host gene expression(Bray-Curtis)'))+
  annotate('text',label=paste0('        Procrusters analysis\n
           M2=',round(pro.m.g_t.2$ss,3), ', p-value=',signif(pro.m.g_t.2$signif,3)),
           x=-0.2,y=0.1,size=4)+
  labs(x='Dimension 1',y='Dimension 2')+
  theme_bw()+
  theme(legend.position = 'bottom')
dev.off()

# Mantel
man.m.g.2 <- mantel(micro.dist.2,gene.dist.2,method = 'spearman', permutations = 9999, na.rm = TRUE)
man.m.g.2


# save(pro.m.g,pro.m.g_t,man.m.g,pro.m.g.2,pro.m.g_t.2,man.m.g.2,file = 'step06_output.Rdata')
