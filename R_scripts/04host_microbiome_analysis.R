rm(list = ls())
getwd()

# procrustes analysis -----------------------------------------------------

{
  library(vegan)
  library(coda.base)
  library(magrittr)
  library(ggsci)
  library(ggplot2)
  }
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
FUN =  function(xxx) {
  (sum(xxx==0)/length(xxx))<0.5 }

# microbiome matrix
microbiome_NR <- read.csv(file = './01clinical_analysis/source_data/Voom-SNM-ACC.csv',header = T,row.names = 1) %>% set_rownames(clinic$PATIENT_ID[match(rownames(.),clinic$names)])  
microbiome_LR <- read.csv(file = './01clinical_analysis/source_data/Voom-SNM-Filter-Likely-ACC.csv',header = T,row.names = 1)  %>% set_rownames(clinic$PATIENT_ID[match(rownames(.),clinic$names)])  
microbiome_CR <- read.csv(file = './01clinical_analysis/source_data/Voom-SNM-Filter-Plate_Center-ACC.csv',header = T,row.names = 1)  %>% set_rownames(clinic$PATIENT_ID[match(rownames(.),clinic$names)])  
microbiome_PR <- read.csv(file = './01clinical_analysis/source_data/Voom-SNM-Filter-Putative-ACC.csv',header = T,row.names = 1)  %>% set_rownames(clinic$PATIENT_ID[match(rownames(.),clinic$names)])  

# microbes <- acc_micro_matrix

# dim(microbes)
#[1]  77 1795

# ms2_patient <- clinic$PATIENT_ID[clinic$kmeans=='MS2']
# microbes.ms2 <- microbes[ms2_patient,]
# dim(microbes.ms2)
# #[1]  43 144


# mRNA exp matrix
exp <- read.table(file = './03transcriptomic_analysis/source_data/TCGA-ACC.htseq_tpm.tsv',header = T,row.names = 1,sep = '\t',check.names = F) %>%
  .[,clinic$PATIENT_ID] %>%
  t() %>%
  .[,apply(., 2, FUN)] %>%
  filter_genes(2)
# exp <- log2(exp+1)
# exp <- read.table(file = './03transcriptomic_analysis/source_data/TCGA-ACC.htseq_tpm.tsv',header = T,row.names = 1,sep = '\t',check.names = F) %>% 
#   .[,clinic$PATIENT_ID]

dim(exp)
#[1]  77 14267




# acc_exp <- read.table(file = './data/TCGA-ACC.htseq_tpm.tsv',header = T,row.names = 1,check.names = F,sep = '\t')
# genes <- t(acc_exp[,rownames(microbes)])
# FUN =  function(xxx) {
#   (sum(xxx==0)/length(xxx))<0.5 }
# genes <- genes[,apply(genes, 2, FUN)]
# genes <- filter_genes(genes,2) #过滤低表达基因
# dim(genes)
# #[1]    77 14267
# genes.ms2 <- genes[ms2_patient,]
# dim(genes.ms2)
# #[1]    43 14267

m2 <- c()
pvalue <- c()
pvalue.man <- c()
microbiome_list <- list(microbiome_NR,microbiome_LR,microbiome_CR,microbiome_PR)
data_type <- c( "Without_combined","Likely_combined","PC_combined","Putative_combined")
procrustes_result <- data.frame(X1=NA,X2=NA,MDS1=NA,MDS2=NA,type=NA)
for (i in 1:4) {
  if(T){
    # ms1_patient <- clinic$PATIENT_ID[clinic[,data_type[1]]=='MS1']
    genes <- exp
    # genes <- exp[,ms1_patient] %>% t() %>% .[,apply(., 2, FUN)] %>% filter_genes(2) 
    
    microbes <- microbiome_list[[i]] %>% .[rownames(genes),]
    microbes[microbes<0] <- 0
    # genes <- as.data.frame(exp)[ms1_patient,]
    # genes <- genes[match(rownames(microbes),rownames(genes)),]
    
    
    microbe.dist <- vegan::vegdist(microbes,method = 'bray')
    gene.dist <- vegan::vegdist(genes,method = 'bray')
    mds.m <- monoMDS(microbe.dist)
    mds.g <- monoMDS(gene.dist)
    
    pro.m.g <- procrustes(mds.g,mds.m,symmetric = TRUE)
    summary(pro.m.g)
    
    # M2 calculating
    set.seed(123)
    pro.m.g_t <- protest(mds.g,mds.m,permutations = 9999)
    # pro.m.g_t
    m2 <- c(m2,pro.m.g_t$ss)
    pvalue <- c(pvalue,round(pro.m.g_t$signif,5))
    # result
    # pro.m.g_t$ss #偏差平方和（M2统计量）
    # pro.m.g_t$signif #对应p值结果
    
    # 获得x和y轴的坐标及旋转过的坐标
    Pro_Y <- cbind(data.frame(pro.m.g$Yrot),data.frame(pro.m.g$X))
    Y <- Pro_Y %>% dplyr::mutate(type=rep(data_type[i]))
    Pro_X <- data.frame(pro.m.g$rotation)
    procrustes_result <- rbind(procrustes_result,Y)
    
    # mantel analysis
    micro.dist <- vegan::vegdist(microbes,method = 'bray')
    gene.dist <- vegan::vegdist(genes,method = 'bray')
    set.seed(123)
    man.m.g <- mantel(micro.dist,gene.dist,method = 'spearman', permutations = 9999, na.rm = TRUE)
    pvalue.man <- c(pvalue.man,round(man.m.g$signif,5))
  }
}

result <- procrustes_result %>% .[-1,] %>% dplyr::mutate(type=as.factor(.$type)) %>% dplyr::mutate(type=factor(.$type,levels = c("Without_combined","Likely_combined","PC_combined","Putative_combined")))

# plot
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
pdf('./04host_microbiome_analysis/output_procrusters_analysis_all_patients.pdf',width = 15,height = 5)
ggplot(data = result,aes(color=type))+
  geom_segment(aes(x=X1,y=X2,xend=(X1+MDS1)/2,yend=(X2+MDS2)/2),
               arrow = arrow(length = unit(0,'cm')),
               size=0.5)+
  geom_segment(aes(x=(X1+MDS1)/2,y=(X2+MDS2)/2,
                   xend=MDS1,yend=MDS2),
               arrow = arrow(length = unit(0.2,'cm')),
               size=0.5)+
  geom_point(aes(X1,X2,shape='17'),size=4)+  #shape放在aes里面要加引号‘’
  geom_point(aes(MDS1,MDS2,shape='1'),size=4)+
  scale_color_npg()+
  thm+
  # scale_color_manual(values=c("Without_combined" = "red", "Likely_combined" = "blue", "PC_combined" = "green" , 'Putative_combined' = "orange"))+
  scale_shape_manual(name = "",
                     values = c('17' = 17,'1' = 1),
                     labels = c('Microbiota abundance(Bray-Curtis)', "Host gene expression(Aitchison's)"))+
  # annotate('text',label=paste0('        Procrusters analysis\n
  #          M2=',round(pro.m.g_t$ss,3), ', p-value=',round(pro.m.g_t$signif,3)),
  #          x=-0.2,y=0.1,size=4)+
  labs(x='Dimension 1',y='Dimension 2')+
  guides(color=guide_legend(direction = "horizontal",nrow=2),
         shape=guide_legend(direction = "horizontal",nrow=2))+
  theme_bw()+
  facet_wrap(.~ type,nrow = 2)+
  theme(legend.position = 'bottom')
dev.off()

m2
#[1] 0.9410252 0.9407154 0.9242710 0.9064394
pvalue
#[1] 0.0194 0.0190 0.0070 0.0022
pvalue.man
#[1] 0.0483 0.0036 0.0001 0.0004




# 
# 
# ## Procrustes analysis
# micro.dist.2 <- vegan::vegdist(microbes.ms2,method = 'robust.aitchison')
# gene.dist.2 <- vegan::vegdist(genes.ms2,method = 'bray')
# mds.m.2 <- monoMDS(micro.dist.2)
# mds.g.2 <- monoMDS(gene.dist.2)
# 
# pro.m.g.2 <- procrustes(mds.g.2,mds.m.2,symmetric = TRUE)
# summary(pro.m.g.2)
# 
# #M2统计量的显著性检验
# set.seed(123)
# pro.m.g_t.2 <- protest(mds.g.2,mds.m.2,permutations = 9999)
# pro.m.g_t.2
# 
# #提取结果
# pro.m.g_t.2$ss #偏差平方和（M2统计量）
# pro.m.g_t.2$signif #对应p值结果
# 
# #获得x和y轴的坐标及旋转过的坐标
# Pro_Y.2 <- cbind(data.frame(pro.m.g.2$Yrot),data.frame(pro.m.g.2$X))
# Pro_X.2 <- data.frame(pro.m.g.2$rotation)
# 
# # plot
# pdf('procrusters_analysis_ms2_patients.pdf')
# ggplot(data = Pro_Y.2)+
#   geom_segment(aes(x=X1,y=X2,xend=(X1+MDS1)/2,yend=(X2+MDS2)/2),
#                arrow = arrow(length = unit(0,'cm')),
#                color='grey',size=0.5)+
#   geom_segment(aes(x=(X1+MDS1)/2,y=(X2+MDS2)/2,
#                    xend=MDS1,yend=MDS2),
#                arrow = arrow(length = unit(0.2,'cm')),
#                color='grey',size=0.5)+
#   geom_point(aes(X1,X2,shape='17'),color='#9BBB59',size=4)+  #shape放在aes里面要加引号‘’
#   geom_point(aes(MDS1,MDS2,shape='1'),color='#9BBB59',size=4)+
#   scale_shape_manual(name = "",
#                      values = c('17' = 17,'1' = 1),
#                      labels = c('Microbiota abundance(Bray-Curtis)', 'Host gene expression(Bray-Curtis)'))+
#   annotate('text',label=paste0('        Procrusters analysis\n
#            M2=',round(pro.m.g_t.2$ss,3), ', p-value=',signif(pro.m.g_t.2$signif,3)),
#            x=-0.2,y=0.1,size=4)+
#   labs(x='Dimension 1',y='Dimension 2')+
#   theme_bw()+
#   theme(legend.position = 'bottom')
# dev.off()
# 
# # Mantel
# man.m.g.2 <- mantel(micro.dist.2,gene.dist.2,method = 'spearman', permutations = 9999, na.rm = TRUE)
# man.m.g.2
# 
# 
# # save(pro.m.g,pro.m.g_t,man.m.g,pro.m.g.2,pro.m.g_t.2,man.m.g.2,file = 'step06_output.Rdata')



# sparseCCA ---------------------------------------------------------------

rm(list = ls()) ## check to make sure no precomputed dataframe before clearing
getwd()
load('step00_output.Rdata')
{
  library(PMA)
  library(data.table)
  library(tidyr)
  library(msigdbr)
  library(data.table)
}

# Functions ---------------------------------------------------------------
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
get_avg_features <- function(cca_cov, CCA.K){
  num_features <- 0
  for(k in 1:CCA.K){
    num_features <- num_features + length(which(cca_cov[,k]!=0))
  }
  avg_features <- num_features/CCA.K
}
FUN =  function(xxx) {
  (sum(xxx==0)/length(xxx))<0.5 }
msigdb_enrichment <- function(pathway_DB, pathways, background_genes, genes_of_interest){
  
  enrichment_list <- list()
  for(i in 1:length(pathways)){
    
    pathway <- pathways[i]
    
    ## genes in this pathway
    pathway_gene_set <- pathway_DB[pathway_DB$gs_name == pathway,]$human_gene_symbol
    length(pathway_gene_set) 
    ## If the criteria for min and max #genes in a given pathway is not satified, 
    ## skip testing the current pathway
    if(length(pathway_gene_set) < min_genes || length(pathway_gene_set) > max_genes) next
    
    ## The contingency table
    ## Inspired by: https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
    ##                    genes_of_interest genes_NOT_of_interest  Total
    ## In_pathway               x                 m - x              m
    ## Not_in_pathway         k - x             n - (k - x)          n
    ## Total                    k                 m + n - k          m + n
    
    ## m, #overlapping genes in this pathway and background genes for CRC
    m <- length(intersect(background_genes,pathway_gene_set))
    m ## 11
    
    ## n, #genes in background but not in pathway 
    n <- length(setdiff(background_genes,pathway_gene_set))
    n #12502
    
    ## x, #genes of interest in pathway 
    x <- length(intersect(pathway_gene_set,genes_of_interest))
    x ## 1
    ## If the overlap between genes of interest and the pathway is less than the overlap cut-off, 
    ## then skip testing the current pathway
    if(x < overlap_genes) next 
    
    ## Extract list of genes in the genes of interest that are included in the pathway. 
    gene_names_in_pathway = noquote(paste(intersect(pathway_gene_set,genes_of_interest), collapse = ","))
    
    ## k, total #genes in genes list of interest
    k <- length(genes_of_interest)
    k
    
    ## Build the contigency table
    contingency_table <- matrix(c(x, k-x, m-x, n-(k-x)), #matrix is filled along the column by default
                                nrow = 2, ncol = 2, 
                                dimnames = list(c("In_pathway","Not_in_pathway"),
                                                c("Genes_of_interest","Genes_NOT_of_interest"))
    )
    
    contingency_table
    
    fisher_result <- fisher.test( contingency_table, alternative = "greater")
    
    ## save details in a dataframe
    enrichment_result_df <- data.frame( pathway = pathway,
                                        BG_genes = length(background_genes),
                                        genes_in_pathway = length(pathway_gene_set),
                                        genes_in_path_and_BG = m,
                                        genes_of_interest = k,
                                        genes_of_interest_in_pathway = x,
                                        gene_names_in_pathway = gene_names_in_pathway,
                                        ## fill in contingency table entries: x, k-x, m-x, n-(k-x) for z-score computation.
                                        cont_n1 = x,
                                        cont_n2 = k-x,
                                        cont_n3 = m-x,
                                        cont_n4 = n-(k-x),
                                        CI_95 = paste0(signif(fisher_result$conf.int[1],5),"-",signif(fisher_result$conf.int[2],5)),
                                        odds_ratio = unname(fisher_result$estimate),
                                        p_val = fisher_result$p.value
    )
    
    enrichment_list[[i]] <- enrichment_result_df
    
  }
  
  return(enrichment_list)
  
}
perform_enrichment <- function(input_dir, filenames, msigdb_collection_name, pathway_DB, pathways, background_genes, output_dir){
  enriched_list <- list()
  count <- 0
  
  ## perform enrichment for each component 
  for(i in filenames){
    
    ## debug
    # i <- "gene_taxa_component_1.txt"
    
    count <- count + 1
    
    print(paste0("Enrichment for component: ", i));flush.console()
    ## load genes list for a given component
    sparseCCA_genes <- read.table(paste0(input_dir,"/",i),sep="\t",head=T, row.names = 1, check.names = F, stringsAsFactors = F)
    sparseCCA_genes <- sparseCCA_genes$gene
    genes_of_interest <- sparseCCA_genes
    
    ## perform enrichment using pathways in msigdb collection
    msigdb_pathways <- msigdb_enrichment(pathway_DB, pathways, background_genes, genes_of_interest) 
    ## convert the list to dataframe
    msigdb_pathways_df <- do.call(rbind.data.frame, msigdb_pathways)
    ## Add current component as a column to enrichment result
    msigdb_pathways_df$Component <- rep(paste0(i),nrow(msigdb_pathways_df))
    ## Add collection name as a column to the enrichemnt result
    msigdb_pathways_df$Collection <- msigdb_collection_name
    
    ## Make last column(i.e., component) as first column
    msigdb_pathways_df <- msigdb_pathways_df[,c(ncol(msigdb_pathways_df)-1, ncol(msigdb_pathways_df),1:(ncol(msigdb_pathways_df)-2))]
    
    ## update column name to reflect component, except pathway colname
    # colnames(msigdb_pathways_df[,-3]) <- paste(colnames(msigdb_pathways_df[,-3]), count, sep = "_")
    
    ## sort pathways by pval
    msigdb_pathways_df <- msigdb_pathways_df[order(msigdb_pathways_df$p_val, decreasing = F),]
    length(msigdb_pathways_df$pathway)
    
    ## save dataframe of pathways for a component to file
    if(!is.null(output_dir)){
      filename <- strsplit(i,"\\.")[[1]][1]
      write.table(msigdb_pathways_df, file = paste0(output_dir,"/msigdb_",filename,".txt") , sep="\t", row.names = F)
    }
    ## add df to list
    enriched_list[[count]] <- msigdb_pathways_df
  }
  
  ## This is the list of dataframes, where each dataframe holds enrichment result for a component.  
  return(enriched_list)
}

# Prepare data ------------------------------------------------------------

## load gene expression and microbiome tables


# microbe abundance matrix
microbiome_NR <- read.csv(file = './01clinical_analysis/source_data/Voom-SNM-ACC.csv',header = T,row.names = 1) %>% set_rownames(clinic$PATIENT_ID[match(rownames(.),clinic$names)])  
microbiome_LR <- read.csv(file = './01clinical_analysis/source_data/Voom-SNM-Filter-Likely-ACC.csv',header = T,row.names = 1)  %>% set_rownames(clinic$PATIENT_ID[match(rownames(.),clinic$names)])  
microbiome_CR <- read.csv(file = './01clinical_analysis/source_data/Voom-SNM-Filter-Plate_Center-ACC.csv',header = T,row.names = 1)  %>% set_rownames(clinic$PATIENT_ID[match(rownames(.),clinic$names)])  
microbiome_PR <- read.csv(file = './01clinical_analysis/source_data/Voom-SNM-Filter-Putative-ACC.csv',header = T,row.names = 1)  %>% set_rownames(clinic$PATIENT_ID[match(rownames(.),clinic$names)])  
microbiome_list <- list(microbiome_NR,microbiome_LR,microbiome_CR,microbiome_PR)

data_type <- c( "Without_combined","Likely_combined","PC_combined","Putative_combined")
# ms1_patient <- clinic$PATIENT_ID[clinic[,data_type[1]]=='MS1']
# microbes <- microbiome_list[[1]] %>% .[ms1_patient,] %>% as.matrix()
# dim(microbes)

# mRNA exp matrix
exp <- read.table(file = './03transcriptomic_analysis/source_data/TCGA-ACC.htseq_tpm.tsv',header = T,row.names = 1,sep = '\t',check.names = F) %>% 
  .[,clinic$PATIENT_ID] %>% 
  t() %>% 
  .[,apply(., 2, FUN)] %>% 
  filter_genes(2)
dim(exp)

## Ensure same sampleIDs in both genes and microbes data
# stopifnot(all(rownames(genes) == rownames(microbes)))

for (k in 4:4) {
  # data preparation
  ms1_patient <- clinic$PATIENT_ID[clinic[,data_type[k]]=='MS1']
  microbes <- microbiome_list[[k]] %>% .[ms1_patient,] %>% as.matrix()
  genes <- exp %>% .[ms1_patient,]
  stopifnot(all(rownames(genes) == rownames(microbes)))
  
  #automatically tuning parameters
  X <- genes
  Z <- microbes
  perm.out <- CCA.permute(X,Z,typex="standard",typez="standard",nperms=10)
  bestpenaltyX <- perm.out$bestpenaltyx
  bestpenaltyZ <- perm.out$bestpenaltyz
  
  # Sparse CCA using selected tuning param
  CCA.out <- CCA(X,Z,typex="standard",typez="standard",K=10,
                 penaltyx=bestpenaltyX,penaltyz=bestpenaltyZ,
                 v=perm.out$v.init)
  rownames(CCA.out$u) <- colnames(X)
  rownames(CCA.out$v) <- colnames(Z)
  
  CCA_var_genes <- X %*% CCA.out$u ## canonical variance for genes 
  CCA_var_microbes <- Z %*% CCA.out$v ## canonical variance for microbes
  
  CCA.out$cors ## canonical correlation for each component:
  
  cca.k <- 10
  avg_genes <- get_avg_features(CCA.out$u, cca.k) ## average number of genes and microbes in resulting components
  avg_genes
  
  avg.microbes <- get_avg_features(CCA.out$v, cca.k)
  avg.microbes
  
  # Test significance using LOOCV
  X <- genes;Y <- microbes;cca.k = 10
  scoresXcv <- matrix(nrow = nrow(X), ncol = cca.k)
  scoresYcv <-  matrix(nrow = nrow(Y), ncol = cca.k)
  corr_pval <- c();corr_r <- c()
  for(i in 1:nrow(genes)){ #n = no. of samples
    #compute weights with sample i held out:
    res <- CCA(X[-i,],Y[-i,], penaltyx=bestpenaltyX, penaltyz=bestpenaltyZ, K=cca.k, trace = F) ## default niter = 15 which is spit out when trace = T (default)
    ###compute scores for i'th sample for each component (pair of canonical variables)
    for(j in 1:cca.k){
      print(paste0("i = ", i," K = ", j)); flush.console()
      scoresXcv[i,j] <- X[i,]%*%res$u[,j]
      scoresYcv[i,j] <- Y[i,]%*%res$v[,j]
    }
  }
  for(j in 1:cca.k){
    # plot(scoresXcv,scoresYcv)
    corr <- cor.test(scoresXcv[,j],scoresYcv[,j])
    corr_pval[j] <- corr$p.value
    corr_r[j] <- corr$estimate
  }
  
  corr_padj <- p.adjust(corr_pval, method = "BH")
  
  ## Spit out significant components
  sig <- which(corr_padj < 0.1)
  dirname <- paste0("./04host_microbiome_analysis/sparse_cca_analysis/sig_gene_taxa_components_",data_type[k],"_",bestpenaltyX,"_", bestpenaltyZ,"_padj/")
  ifelse(!dir.exists(dirname), dir.create(dirname), FALSE)
  
  for(i in sig){
    selected_X <- which(CCA.out$u[,i]!=0) 
    selected_X <- rownames(CCA.out$u)[selected_X]
    coeff_X <- unname(CCA.out$u[selected_X,i])
    selected_Z <- which(CCA.out$v[,i]!=0)
    selected_Z <- rownames(CCA.out$v)[selected_Z]
    coeff_Z <- unname(CCA.out$v[selected_Z,i])
    ## Make all vectors of same length to avoid repetition of elements from shorter vectors.
    n <- max(length(selected_X), length(selected_Z))
    length(selected_X) <- n                      
    length(selected_Z) <- n
    length(coeff_X) <- n
    length(coeff_Z) <- n
    selected_XZ <- as.data.frame(cbind(gene = selected_X, gene_coeff = coeff_X,
                                       taxa = selected_Z, taxa_coeff = coeff_Z))
    write.table(selected_XZ, file=paste0(dirname,"gene_taxa_component_",i,".txt"), sep = "\t", col.names = NA)
  }
  
  # Enrichment_sparseCCA
  
  ##background genes
  background.genes  <- colnames(genes)
  background.genes <- background.genes[order(background.genes)]
  
  msigdb_C2 <- msigdbr(species = "Homo sapiens", category = "C2")
  msigdb_C2 <- as.data.frame(msigdb_C2)
  table(msigdb_C2$gs_subcat)
  
  ## Only keep canonical pathways, or in other words, remove CGP (Chemical and genetic perturbations)
  msigdb_C2_CP <- msigdb_C2[msigdb_C2$gs_subcat!= "CGP",]
  table(msigdb_C2_CP$gs_subcat)
  dim(msigdb_C2_CP)
  #[1] 168417     15
  length(unique(msigdb_C2_CP$gs_name))
  
  ## Pathway DB of interest: KEGG PID
  path_DB <- c("CP:KEGG","CP:PID")
  msigdb_C2_CP <- msigdb_C2_CP[(msigdb_C2_CP$gs_subcat %in% path_DB), ]
  table(msigdb_C2_CP$gs_subcat)
  
  ## unique pathways in DBs of interest
  msigdb_C2_CP_unique_pathways <- unique(msigdb_C2_CP$gs_name)
  length(msigdb_C2_CP_unique_pathways) 
  #[1] 186
  
  ## Set filtering criteria for all pathway DBs
  min_genes <- 25
  max_genes <- 300
  overlap_genes <- 5
  
  ## Perform enrichment per component
  msigdb_collection_name <- "C2_CP" 
  input_enrichment <- paste0("sig_gene_taxa_components_",data_type[k],"_",bestpenaltyX,"_", bestpenaltyZ,"_padj")
  input_dir <- paste0("./04host_microbiome_analysis/sparse_cca_analysis/",input_enrichment)
  filenames <- list.files(input_dir)
  stopifnot(length(filenames) > 0) # make sure the input dir path is set correctly
  pathway_DB <- msigdb_C2_CP
  pathways <- msigdb_C2_CP_unique_pathways
  background_genes <- background.genes
  output_dir <- paste0("./04host_microbiome_analysis/enrichment_sparseCCA_msigDB/per_component_kegg_",data_type[k],"/")
  ## create the output dir if it does not exist
  ifelse(!dir.exists(output_dir), dir.create(output_dir, recursive = T), FALSE)
  
  ## Perform enrichment
  enrichment_list <- perform_enrichment(input_dir, filenames, 
                                        msigdb_collection_name, pathway_DB, pathways, background_genes, 
                                        output_dir)
  
  ## Pool enrichment output across components.
  enrichment_case <- do.call(rbind,enrichment_list)
  dim(enrichment_case)
  # [1] 314   16
  
  any(duplicated(enrichment_case$pathway))
  
  ## order by p-value
  enrichment_case <- enrichment_case[order(enrichment_case$p_val),]
  
  ##filter duplicates (keeps pathways with smallest p-value across components)
  enrichment_case_no_dups <- enrichment_case[!duplicated(enrichment_case$pathway),]
  dim(enrichment_case_no_dups) 
  # [1] 116  16
  
  ## MHT correction -- FDR
  enrichment_case_no_dups$p_adj <- p.adjust(enrichment_case_no_dups$p_val, method = "BH")
  enrichment_case_no_dups <- enrichment_case_no_dups[order(enrichment_case_no_dups$p_adj),]
  length(which(enrichment_case_no_dups$p_adj < 0.1))
  # [1] 101
  
  # dir.create(paste0("./04host_microbiome_analysis/enrichment_sparseCCA_msigDB/enrichment_kegg_",data_type[k],"/"), recursive = T)
  write.table(enrichment_case_no_dups,file = paste0("./04host_microbiome_analysis/enrichment_sparseCCA_msigDB/enrichment_kegg_",data_type[k],'.txt'),sep="\t", row.names = T,col.names = NA)
  ## write to file
  # filename <- paste0("ACC_ms2_msigdb.txt")
  # filepath <- paste0("./data/sparse_cca_analysis/ACC/output_sparseCCA_V2/final_gene_taxa_components/enrichment_sparseCCA_msigDB/",msigdb_collection_name,"/",filename)
  # write.table(enrichment_case, file = paste0(filepath) , sep="\t", row.names = F)
  
  
}

# Visualization of specific pathway ----------------------------------------------------
pathway_name <- 'KEGG_CELL_CYCLE'
pathway_name <- 'PID_P53_REGULATION_PATHWAY'
max_min<-function(x){
  return((x-min(x))/(max(x)-min(x)))}
enrichment <- read.table("./04host_microbiome_analysis/enrichment_sparseCCA_msigDB/enrichment_kegg_Likely_combined.txt",sep="\t",head=T, row.names = 1, check.names = F, stringsAsFactors = F)
component_name <- enrichment$Component[enrichment$pathway==pathway_name]
component <- read.table("./04host_microbiome_analysis/sparse_cca_analysis/sig_gene_taxa_components_Likely_combined_0.7_0.7_padj/gene_taxa_component_1.txt",sep="\t",head=T, row.names = 1, check.names = F, stringsAsFactors = F)


microbiota_name <- as.character(na.omit(component$taxa))
microbiota_coeff_abs <- as.numeric(abs(na.omit(component$taxa_coeff)))
names(microbiota_coeff_abs) <- microbiota_name
microbiota_coeff_abs <- max_min(microbiota_coeff_abs) 
sort(microbiota_coeff_abs,decreasing = T)[1:10]


gene_name <- enrichment$gene_names_in_pathway[enrichment$pathway==pathway_name]
gene_name <- as.character(unlist(strsplit(gene_name,',')))
component.1 <- component[component$gene %in% gene_name,]
gene_coeff_abs <- abs(na.omit(component.1$gene_coeff))
names(gene_coeff_abs) <- gene_name
gene_coeff_abs <- max_min(gene_coeff_abs)
sort(gene_coeff_abs,decreasing = T)[1:10]



# plot
#spearman correlation
ms1_patient <- clinic$PATIENT_ID[clinic[,data_type[2]]=='MS1']
microbes <- microbiome_list[[2]] %>% .[ms1_patient,] %>% as.matrix()
genes <- exp %>% .[ms1_patient,]
cor <- cbind(as.data.frame(microbes[ms1_patient,names(sort(microbiota_coeff_abs,decreasing = T)[1:10])]) %>% set_colnames(str_sub(colnames(.),str_locate(colnames(.),'g__')[,1],str_length(colnames(.)))),
            genes[ms1_patient,names(sort(gene_coeff_abs,decreasing = T)[1:10])])
microbiota_name <- str_sub(microbiota_name,str_locate(microbiota_name,'g__')[,1],str_length(microbiota_name))
library(ggpubr)
library(gridExtra)
library(Hmisc)
cor_res <- list()
for (i in 1:10) {
  for (j in 11:20) {
    model <- rcorr(cor[,i],cor[,j], type  = "pearson")
    # model <- cor.test(cor[,ncol(cor)],cor[,i],method = 'spearman')
    text <- paste0('Pearson,','rho=',signif(model$r,3),',p-value=',signif(model$P,3))
    p <- ggplot(cor,aes_string(x =colnames(cor)[i],y =colnames(cor)[j]))+
      geom_point(size=5,color = '#E18727CC')+
      geom_smooth(method = 'lm', color = 'black',size=2)+
      stat_cor(data=cor, method = "pearson")+
      thm+
      theme_bw()
    cor_res[[i+j]] <- p
  }
  # cor_res
  # model <- rcorr(as.matrix(cor), type  = "pearson")
  # # model <- cor.test(cor[,ncol(cor)],cor[,i],method = 'spearman')
  # text <- paste0('Pearson,','rho=',signif(model$r,3),',p-value=',signif(model$p,3))
  # p <- ggplot(cor,aes_string(x =colnames(cor)[i],y ='P53'))+
  #   geom_point(size=5,color = '#E18727CC')+
  #   geom_smooth(method = 'lm', color = 'black',size=2)+
  #   stat_cor(data=cor, method = "pearson")+
  #   theme_bw()
  # cor_res[[i]] <- p
}

p <- grid.arrange(grobs=cor_res,ncol=5)
pdf('./04host_microbiome_analysis/output_correlation_cell_cycle_pathway.pdf',width = 20,height = 20)
grid.arrange(grobs=cor_res,ncol=5)
dev.off()







# ms1_patient <- clinic$PATIENT_ID[clinic$kmeans=='MS2']
# pathway_score<- t(gsva_kegg_score[,ms2_patient])
# cor_res <- list()
# cor <- cbind(as.data.frame(microbes[ms2_patient,microbiota_name]),as.data.frame(pathway_score[,pathway_name]))
# colnames(cor)[ncol(cor)] <- c('P53')
# library(ggpubr)
# library(gridExtra)
# for (i in 1:(ncol(cor)-1)) {
#   model <- cor.test(cor[,ncol(cor)],cor[,i],method = 'spearman')
#   text <- paste0('Spearman,','rho=',signif(model$estimate,3),',p-value=',signif(model$p.value,3))
#   p <- ggplot(cor,aes_string(x =colnames(cor)[i],y ='P53'))+
#     geom_point(size=5,color = '#E18727CC')+
#     geom_smooth(method = 'lm', color = 'black',size=2)+
#     stat_cor(data=cor, method = "spearman")+
#     theme_bw()
#   cor_res[[i]] <- p
# }
# p <- grid.arrange(grobs=cor_res,ncol=5)
# p




# 
# X <- genes
# Z <- microbes
# 
# perm.out <- CCA.permute(X,Z,typex="standard",typez="standard",nperms=10)
# print(perm.out)
# plot(perm.out)
# bestpenaltyX <- perm.out$bestpenaltyx
# bestpenaltyZ <- perm.out$bestpenaltyz
# 
# 
# CCA.out <- CCA(X,Z,typex="standard",typez="standard",K=10,
#            penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz,
#            v=perm.out$v.init)


# Sparse CCA tuning parameters -------------------------------------------------

## select tuning parameters using grid-search
if(T){
  X <- genes
  Y <- microbes
  scoreXcv <- c()
  scoreYcv <- c()
  penaltyX <- seq(0.1,0.4,length=10)
  penaltyY <- seq(0.15,0.4,length=10)
  corr_ACC <- matrix(nrow = length(penaltyX), ncol =  length(penaltyY))
  num_samples <- nrow(genes)
  start_time <- Sys.time()
  for( i in 1:length(penaltyX)){
    for(j in 1:length(penaltyY)){
      # print(paste0("Index: i = ",i,", j =", j)); flush.console()
      for(k in 1:num_samples){
        print(paste0("Index: i = ",i,", j =", j," k = ",k)); flush.console()
        #compute weights with sample k held out:
        res <- CCA(X[-k,],Y[-k,], penaltyx = penaltyX[i], penaltyz = penaltyY[j], K=1, niter = 5, trace = F)
        ## Compute scores for k'th sample for first pair of canonical variables
        ## Take weight of features (res$u and res$v) computed using all except 
        ## the kth sample and multiple it by values for the kth sample in the 
        ## feature matrix X and Y. 
        ## %*% implies matrix multiplication. 
        scoreXcv[k] <- X[k,]%*%res$u ## single value
        scoreYcv[k] <- Y[k,]%*%res$v ## single value
      }
      ## correlation between scores for X and Y for all held out samples.
      corr_ACC[i,j] = cor(scoreXcv,scoreYcv) 
    }
  }
  row.names(corr_ACC) <- as.character(penaltyX)
  colnames(corr_ACC) <- as.character(penaltyY)
  
  corr_ACC_df <- as.data.frame(corr_ACC)
  rownames(corr_ACC_df)
  colnames(corr_ACC_df)
}

## save to file
# save(corr_ACC,corr_ACC_df, file = paste0("./data/analysis/",'ACC',"/output_sparseCCA_V2/grid_search/corr_ACC_0.1_0.4_V4.RData"))

## load precomputed grid-search output
# load( file = paste0("./data/sparse_cca_analysis/ACC/output_sparseCCA_V2/grid_search/corr_ACC_0.1_0.4_V4.RData"))
# penaltyX <- seq(0.1,0.4,length=10)
# penaltyY <- seq(0.15,0.4,length=10)
# 
# # find index with max absolute corr
# bestpenalty <- which(abs(corr_ACC) == max(abs(corr_ACC)), arr.ind = TRUE)
# bestpenalty
# 
# bestpenaltyX <- penaltyX[bestpenalty[1]]
# bestpenaltyX 
# bestpenaltyY <- penaltyY[bestpenalty[2]]
# bestpenaltyY
# 
# ## order abs corr to get top 5 corr
# index <- order(abs(corr_ACC), decreasing = T)
# abs(corr_ACC)[index][1:5] ## top 5 absolute corr
# [1] 0.9196399 0.9195552 0.9176733 0.9176085 0.9164174

# cca.k = 10 ## number of desired components


# Sparse CCA using selected tuning param ----------------------------------
# X <- genes
# Z <- microbes
# CCA.out <-  CCA(X,Z,typex="standard",typez="standard",K=10,
#                 penaltyx=bestpenaltyX,penaltyz=bestpenaltyY) ## standardize=T by default

## add rownames to output factors
rownames(CCA.out$u) <- colnames(X)
rownames(CCA.out$v) <- colnames(Z)

## compute contribution of selected features to each of the samples.
CCA_var_genes <- X %*% CCA.out$u ## canonical variance for genes 
CCA_var_microbes <- Z %*% CCA.out$v ## canonical variance for microbes

## canonical correlation for each component:
CCA.out$cors

## average number of genes and microbes in resulting components
avg_genes <- get_avg_features(CCA.out$u, cca.k)
avg_genes

avg.microbes <- get_avg_features(CCA.out$v, cca.k)
avg.microbes


# Test significance using LOOCV -------------------------------------------

X <- genes
Y <- microbes
cca.k = 10
scoresXcv <- matrix(nrow = nrow(X), ncol = cca.k)
scoresYcv <-  matrix(nrow = nrow(Y), ncol = cca.k)
corr_pval <- c()
corr_r <- c()
for(i in 1:nrow(genes)){ #n = no. of samples
  #compute weights with sample i held out:
  res <- CCA(X[-i,],Y[-i,], penaltyx=bestpenaltyX, penaltyz=bestpenaltyZ, K=cca.k, trace = F) ## default niter = 15 which is spit out when trace = T (default)
  ###compute scores for i'th sample for each component (pair of canonical variables)
  for(j in 1:cca.k){
    print(paste0("i = ", i," K = ", j)); flush.console()
    scoresXcv[i,j] <- X[i,]%*%res$u[,j]
    scoresYcv[i,j] <- Y[i,]%*%res$v[,j]
  }
}

## Test for each components
for(j in 1:cca.k){
  # plot(scoresXcv,scoresYcv)
  corr <- cor.test(scoresXcv[,j],scoresYcv[,j])
  corr_pval[j] <- corr$p.value
  corr_r[j] <- corr$estimate
}

corr_pval
which(corr_pval < 0.1)
which(corr_pval < 0.05)

corr_padj <- p.adjust(corr_pval, method = "BH")
corr_padj
which(corr_padj < 0.05)
length(which(corr_padj < 0.05))

## LOOCV corr
corr_r

## Spit out significant components
sig <- which(corr_padj < 0.1)
dirname <- paste0("./04host_microbiome_analysis/sparse_cca_analysis/sig_gene_taxa_components_",bestpenaltyX,"_", bestpenaltyZ,"_padj/",data_type[k])
ifelse(!dir.exists(dirname), dir.create(dirname), FALSE)

for(i in sig){
  selected_X <- which(CCA.out$u[,i]!=0) 
  selected_X <- rownames(CCA.out$u)[selected_X]
  coeff_X <- unname(CCA.out$u[selected_X,i])
  selected_Z <- which(CCA.out$v[,i]!=0)
  selected_Z <- rownames(CCA.out$v)[selected_Z]
  coeff_Z <- unname(CCA.out$v[selected_Z,i])
  ## Make all vectors of same length to avoid repetition of elements from shorter vectors.
  n <- max(length(selected_X), length(selected_Z))
  length(selected_X) <- n                      
  length(selected_Z) <- n
  length(coeff_X) <- n
  length(coeff_Z) <- n
  selected_XZ <- as.data.frame(cbind(gene = selected_X, gene_coeff = coeff_X,
                                     taxa = selected_Z, taxa_coeff = coeff_Z))
  write.table(selected_XZ, file=paste0(dirname,"gene_taxa_component_",i,".txt"), sep = "\t", col.names = NA)
}


# Enrichment_sparseCCA ----------------------------------------------------

{
  library(msigdbr)
  library(data.table)
  }
## Do enrichment per component
msigdb_enrichment <- function(pathway_DB, pathways, background_genes, genes_of_interest){
  
  enrichment_list <- list()
  for(i in 1:length(pathways)){
    
    pathway <- pathways[i]
    
    ## genes in this pathway
    pathway_gene_set <- pathway_DB[pathway_DB$gs_name == pathway,]$human_gene_symbol
    length(pathway_gene_set) 
    ## If the criteria for min and max #genes in a given pathway is not satified, 
    ## skip testing the current pathway
    if(length(pathway_gene_set) < min_genes || length(pathway_gene_set) > max_genes) next
    
    ## The contingency table
    ## Inspired by: https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
    ##                    genes_of_interest genes_NOT_of_interest  Total
    ## In_pathway               x                 m - x              m
    ## Not_in_pathway         k - x             n - (k - x)          n
    ## Total                    k                 m + n - k          m + n
    
    ## m, #overlapping genes in this pathway and background genes for CRC
    m <- length(intersect(background_genes,pathway_gene_set))
    m ## 11
    
    ## n, #genes in background but not in pathway 
    n <- length(setdiff(background_genes,pathway_gene_set))
    n #12502
    
    ## x, #genes of interest in pathway 
    x <- length(intersect(pathway_gene_set,genes_of_interest))
    x ## 1
    ## If the overlap between genes of interest and the pathway is less than the overlap cut-off, 
    ## then skip testing the current pathway
    if(x < overlap_genes) next 
    
    ## Extract list of genes in the genes of interest that are included in the pathway. 
    gene_names_in_pathway = noquote(paste(intersect(pathway_gene_set,genes_of_interest), collapse = ","))
    
    ## k, total #genes in genes list of interest
    k <- length(genes_of_interest)
    k
    
    ## Build the contigency table
    contingency_table <- matrix(c(x, k-x, m-x, n-(k-x)), #matrix is filled along the column by default
                                nrow = 2, ncol = 2, 
                                dimnames = list(c("In_pathway","Not_in_pathway"),
                                                c("Genes_of_interest","Genes_NOT_of_interest"))
    )
    
    contingency_table
    
    fisher_result <- fisher.test( contingency_table, alternative = "greater")
    
    ## save details in a dataframe
    enrichment_result_df <- data.frame( pathway = pathway,
                                        BG_genes = length(background_genes),
                                        genes_in_pathway = length(pathway_gene_set),
                                        genes_in_path_and_BG = m,
                                        genes_of_interest = k,
                                        genes_of_interest_in_pathway = x,
                                        gene_names_in_pathway = gene_names_in_pathway,
                                        ## fill in contingency table entries: x, k-x, m-x, n-(k-x) for z-score computation.
                                        cont_n1 = x,
                                        cont_n2 = k-x,
                                        cont_n3 = m-x,
                                        cont_n4 = n-(k-x),
                                        CI_95 = paste0(signif(fisher_result$conf.int[1],5),"-",signif(fisher_result$conf.int[2],5)),
                                        odds_ratio = unname(fisher_result$estimate),
                                        p_val = fisher_result$p.value
    )
    
    enrichment_list[[i]] <- enrichment_result_df
    
  }
  
  return(enrichment_list)
  
}
perform_enrichment <- function(input_dir, filenames, msigdb_collection_name, pathway_DB, pathways, background_genes, output_dir){
  enriched_list <- list()
  count <- 0
  
  ## perform enrichment for each component 
  for(i in filenames){
    
    ## debug
    # i <- "gene_taxa_component_1.txt"
    
    count <- count + 1
    
    print(paste0("Enrichment for component: ", i));flush.console()
    ## load genes list for a given component
    sparseCCA_genes <- read.table(paste0(input_dir,"/",i),sep="\t",head=T, row.names = 1, check.names = F, stringsAsFactors = F)
    sparseCCA_genes <- sparseCCA_genes$gene
    genes_of_interest <- sparseCCA_genes
    
    ## perform enrichment using pathways in msigdb collection
    msigdb_pathways <- msigdb_enrichment(pathway_DB, pathways, background_genes, genes_of_interest) 
    ## convert the list to dataframe
    msigdb_pathways_df <- do.call(rbind.data.frame, msigdb_pathways)
    ## Add current component as a column to enrichment result
    msigdb_pathways_df$Component <- rep(paste0(i),nrow(msigdb_pathways_df))
    ## Add collection name as a column to the enrichemnt result
    msigdb_pathways_df$Collection <- msigdb_collection_name
    
    ## Make last column(i.e., component) as first column
    msigdb_pathways_df <- msigdb_pathways_df[,c(ncol(msigdb_pathways_df)-1, ncol(msigdb_pathways_df),1:(ncol(msigdb_pathways_df)-2))]
    
    ## update column name to reflect component, except pathway colname
    # colnames(msigdb_pathways_df[,-3]) <- paste(colnames(msigdb_pathways_df[,-3]), count, sep = "_")
    
    ## sort pathways by pval
    msigdb_pathways_df <- msigdb_pathways_df[order(msigdb_pathways_df$p_val, decreasing = F),]
    length(msigdb_pathways_df$pathway)
    
    ## save dataframe of pathways for a component to file
    if(!is.null(output_dir)){
      filename <- strsplit(i,"\\.")[[1]][1]
      write.table(msigdb_pathways_df, file = paste0(output_dir,"/msigdb_",filename,".txt") , sep="\t", row.names = F)
    }
    ## add df to list
    enriched_list[[count]] <- msigdb_pathways_df
  }
  
  ## This is the list of dataframes, where each dataframe holds enrichment result for a component.  
  return(enriched_list)
}

##background genes
background.genes  <- colnames(genes)
background.genes <- background.genes[order(background.genes)]
length(background.genes)

## Get MSigDB pathways
## http://software.broadinstitute.org/gsea/msigdb/collection_details.jsp#C2

msigdb_C2 <- msigdbr(species = "Homo sapiens", category = "C2")
msigdb_C2 <- as.data.frame(msigdb_C2)
table(msigdb_C2$gs_subcat)

## Only keep canonical pathways, or in other words, remove CGP (Chemical and genetic perturbations)
msigdb_C2_CP <- msigdb_C2[msigdb_C2$gs_subcat!= "CGP",]
table(msigdb_C2_CP$gs_subcat)
dim(msigdb_C2_CP)
#[1] 168417     15
length(unique(msigdb_C2_CP$gs_name))

## Pathway DB of interest: KEGG
path_DB <- c("CP:KEGG")
msigdb_C2_CP <- msigdb_C2_CP[(msigdb_C2_CP$gs_subcat %in% path_DB), ]
table(msigdb_C2_CP$gs_subcat)

## unique pathways in DBs of interest
msigdb_C2_CP_unique_pathways <- unique(msigdb_C2_CP$gs_name)
length(msigdb_C2_CP_unique_pathways) 
#[1] 186

## Set filtering criteria for all pathway DBs
min_genes <- 25
max_genes <- 300
overlap_genes <- 5

## Perform enrichment per component
msigdb_collection_name <- "C2_CP" 
input_enrichment <- "sig_gene_taxa_components_0.1_0.1_padj"
input_dir <- paste0("./04host_microbiome_analysis/sparse_cca_analysis/",input_enrichment)
filenames <- list.files(input_dir)
stopifnot(length(filenames) > 0) # make sure the input dir path is set correctly
pathway_DB <- msigdb_C2_CP
pathways <- msigdb_C2_CP_unique_pathways
background_genes <- background.genes
output_dir <- paste0("./04host_microbiome_analysis/enrichment_sparseCCA_msigDB/")
## create the output dir if it doesnot exist
ifelse(!dir.exists(output_dir), dir.create(output_dir, recursive = T), FALSE)

## Perform enrichment
enrichment_list <- perform_enrichment(input_dir, filenames, 
                                      msigdb_collection_name, pathway_DB, pathways, background_genes, 
                                      output_dir)

## Pool enrichment output across components.
enrichment_case <- do.call(rbind,enrichment_list)
dim(enrichment_case)
# [1] 314   16

any(duplicated(enrichment_case$pathway))

## order by p-value
enrichment_case <- enrichment_case[order(enrichment_case$p_val),]

##filter duplicates (keeps pathways with smallest p-value across components)
enrichment_case_no_dups <- enrichment_case[!duplicated(enrichment_case$pathway),]
dim(enrichment_case_no_dups) 
# [1] 116  16

## MHT correction -- FDR
enrichment_case_no_dups$p_adj <- p.adjust(enrichment_case_no_dups$p_val, method = "BH")
enrichment_case_no_dups <- enrichment_case_no_dups[order(enrichment_case_no_dups$p_adj),]
length(which(enrichment_case_no_dups$p_adj < 0.1))
# [1] 101

## write to file
# filename <- paste0("ACC_ms2_msigdb.txt")
# filepath <- paste0("./data/sparse_cca_analysis/ACC/output_sparseCCA_V2/final_gene_taxa_components/enrichment_sparseCCA_msigDB/",msigdb_collection_name,"/",filename)
# write.table(enrichment_case, file = paste0(filepath) , sep="\t", row.names = F)

#ggplot dotplot for enrichment result
kegg <- enrichment_case_no_dups
kegg$GeneRatio<-kegg$genes_of_interest_in_pathway/kegg$genes_in_pathway
kegg <- kegg[order(kegg$GeneRatio,decreasing = T),]

#intersection between results from GSEA and KEGG
load('step03_output.Rdata')
gsea_result <- KEGG_gseresult2@result %>% 
  filter(NES>0) %>% 
  mutate(Description=toupper(.$Description)) %>% 
  mutate(Description=gsub('[ ]','_',.$Description)) %>% 
  mutate(Description=paste0('KEGG_',.$Description))

inter_pathway <- intersect(gsea_result$Description,kegg$pathway)
length(inter_pathway)

kegg_intersect <- kegg[kegg$pathway%in%inter_pathway,]
kegg_intersect <- kegg_intersect[order(kegg_intersect$p_adj),]
kegg_intersect <- kegg_intersect[!str_detect(kegg_intersect$pathway,c('CANCER|CARCINOMA|LEUKEMIA')),]

# plot
pdf('enrichment_gsea_sparse.pdf')
ggplot(kegg_intersect,aes(x=GeneRatio,y=pathway))+
  geom_point(aes(size=genes_of_interest_in_pathway,color=p_adj))+
  scale_color_gradient(low="#FF3030",high ="#104E8B")+
  theme_bw()
dev.off()


# Visualization of P53 ----------------------------------------------------

pathway_name <- 'KEGG_P53_SIGNALING_PATHWAY'
FUNCTION<-function(x){
  return((x-min(x))/(max(x)-min(x)))}

component_name <- kegg_intersect$Component[kegg_intersect$pathway==pathway_name]
component <- read.table(paste0("./data/sparse_cca_analysis/ACC/output_sparseCCA_V2/grid_search/gene_taxa_components/sig_gene_taxa_components_0.2_0.15_padj/",component_name),sep="\t",head=T, row.names = 1, check.names = F, stringsAsFactors = F)
microbiota_name <- as.character(na.omit(component$taxa))
microbiota_coeff_abs <- as.numeric(abs(na.omit(component$taxa_coeff)))
names(microbiota_coeff_abs) <- microbiota_name
microbiota_coeff_abs
microbiota_coeff_abs <- FUNCTION(microbiota_coeff_abs) 
sort(microbiota_coeff_abs,decreasing = T)[1:10]

gene_name <- kegg_intersect$gene_names_in_pathway[kegg_intersect$pathway==pathway_name]
gene_name <- as.character(unlist(strsplit(gene_name,',')))
component.1 <- component[component$gene %in% gene_name,]
gene_coeff_abs <- abs(na.omit(component.1$gene_coeff))
names(gene_coeff_abs) <- gene_name
gene_coeff_abs
gene_coeff_abs <- FUNCTION(gene_coeff_abs)
sort(gene_coeff_abs,decreasing = T)[1:10]

# #spearman correlation
# ms2_patient <- clinic$PATIENT_ID[clinic$kmeans=='MS2']
# pathway_score<- t(gsva_kegg_score[,ms2_patient])
# cor_res <- list()
# cor <- cbind(as.data.frame(microbes[ms2_patient,microbiota_name]),as.data.frame(pathway_score[,pathway_name]))
# colnames(cor)[ncol(cor)] <- c('P53')
# library(ggpubr)
# library(gridExtra)
# for (i in 1:(ncol(cor)-1)) {
#   model <- cor.test(cor[,ncol(cor)],cor[,i],method = 'spearman')
#   text <- paste0('Spearman,','rho=',signif(model$estimate,3),',p-value=',signif(model$p.value,3))
#   p <- ggplot(cor,aes_string(x =colnames(cor)[i],y ='P53'))+
#     geom_point(size=5,color = '#E18727CC')+
#     geom_smooth(method = 'lm', color = 'black',size=2)+
#     stat_cor(data=cor, method = "spearman")+
#     theme_bw()
#   cor_res[[i]] <- p
# }
# p <- grid.arrange(grobs=cor_res,ncol=5)
# p










