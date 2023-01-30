rm(list = ls())
getwd()
load('step00_output.Rdata')

# Clinical parameters & MS ------------------------------------------------
{
  library(stringr)
  library(dplyr)
  library(plyr)
  library(magrittr)
  }

variable <- c('SEX',
              'AGE_new',
              # 'RACE',
              'Ethnicity.Race',
              'Clinical.Stage',
              'T.stage',
              'Lymph.node',
              'ATYPICAL_MITOTIC_FIGURES',
              'CAPSULAR_INVASION',
              'Clinical.Status.3.Mo.Post.Op',
              'Distant.Metastasis',
              'CYTOPLASM_PRESENCE_LESS_THAN_EQUAL_25_PERCENT',
              'DIFFUSE_ARCHITECTURE',
              'Excessive.Hormone',
              'LATERALITY',
              'METASTATIC_SITE_1',
              'METASTATIC_SITE_2',
              'MITOTIC_RATE',
              'NECROSIS',
              'NUCLEAR_GRADE',
              'PHARMACEUTICAL_TX_ADJUVANT',
              'PHARM_TX_MITOTANE_ADJUVANT',
              'PHARM_TX_MITOTANE_THERAPUTIC_AT_REC',
              'RADIATION_TREATMENT_ADJUVANT',
              'Surgical.margin',
              'SINUSOID_INVASION',
              'PRIMARY_THERAPY_OUTCOME',
              'Neoplasm.Status',
              'WEISS_VENOUS_INVASION',
              'TUMOR_WEIGHT_new')
data_type <- c('Without_combined',"Likely_combined","PC_combined","Putative_combined",'Strigent_combined')

# Tables
for (j in 1:5) {
  # result <- data.frame(MS1=NA,MS2=NA)
  result <- data.frame()
  for (i in 1:length(variable)) {
    name <- variable[i]
    clinic <- dplyr::mutate(clinic,MS=clinic[,data_type[j]])
    dat <- data.frame(table(clinic[,c('MS',name)])) %>% 
      ddply(.(MS),transform,percent=Freq/sum(Freq)*100) 
    colnames(dat)[1:2] <- c('Var1','Var2')
    x <- matrix(dat$Freq,ncol=2) %>% as.data.frame() %>% set_colnames(c('MS1','MS2'))
    x <- tibble::rownames_to_column(x,var = 'parameter') %>% dplyr::mutate(parameter=unique(dat$Var2))
    result <- rbind(result,x)
  }
  write.table(result,file =paste0( './01clinical_analysis/output_ACC_params_&_MS_',data_type[j],'(table).tsv'),sep = '\t',row.names = T,col.names = NA)
}

# P-value table
result <- data.frame(Parameter=variable,
                     Without_pvalue=NA,
                     Likely_pvalue=NA,
                     PC_pvalue=NA,
                     Putative_pvalue=NA,
                     Strigent_pvalue=NA)
for (i in 1:length(variable)) {
  for (j in 1:5) {
    name <- variable[i]
    clinic <- dplyr::mutate(clinic,MS=clinic[,data_type[j]])
    dat <- data.frame(table(clinic[,c('MS',name)])) %>% 
      ddply(.(MS),transform,percent=Freq/sum(Freq)*100) 
    colnames(dat)[1:2] <- c('Var1','Var2')
    x <- matrix(dat$Freq,ncol=2)
    if (sum(x<5)>0){
      a <- fisher.test(x)
    }else{
      a <- chisq.test(x)
    }
    result[i,j+1] <- round(a$p.value,2)
  }
}
result[result<0.05] <- '*'
write.table(result,file = './01clinical_analysis/output_ACC_params_MS(chi-squared test).tsv',sep = '\t',row.names = T,col.names = NA)


# Multivariate COX --------------------------------------------------------

{
  library(survival)
  library(dplyr)
}
result <- data.frame()
res.cox.all <- coxph(formula = Surv(OS.time,OS)~ 
                       # Without_combined +
                       # Likely_combined+
                       # PC_combined+
                       Putative_combined+
                       Clinical.Stage+
                       # Excessive.Hormone+
                       ATYPICAL_MITOTIC_FIGURES+
                       Surgical.margin+
                       Clinical.Status.3.Mo.Post.Op+
                       Distant.Metastasis+
                       # Neoplasm.Status+
                       PRIMARY_THERAPY_OUTCOME,
                       data = clinic)
x <- summary(res.cox.all)
result_cox <- data.frame(HR=signif(x$conf.int[,1],digits=2),
                         HR.confit.lower=signif(x$conf.int[,"lower .95"],2),
                         HR.confit.upper=signif(x$conf.int[,"upper .95"],2),
                         p.value =signif(x$coefficients[,"Pr(>|z|)"],digits = 3))
result_cox <- result_cox %>%
  mutate('HR(95%CI)'=paste(.$HR,"(",result_cox$HR.confit.lower,"-",.$HR.confit.upper,")",sep = "")) %>%
  dplyr::mutate(p.value=as.numeric(.$p.value)) %>%
  dplyr::mutate(p.value= ifelse(.$p.value<0.001,'***',ifelse(.$p.value<0.01,'**',ifelse(.$p.value<0.05,'*','NS'))))%>%
  tibble::rownames_to_column(var = 'Parameter')
result_cox[,2:4] <- lapply(result_cox[,2:4],as.numeric)
result <- rbind(result,result_cox)
write.table(result,file = './01clinical_analysis/ouput_multicox_result.tsv',sep = '\t',row.names = T,col.names = NA)


dat <- read.table(file = './01clinical_analysis/ouput_multicox_result.tsv',sep = '\t',header = T)
dat$` ` <- paste(rep(" ", nrow(dat)), collapse = " ")
dat$Parameter <- ifelse(is.na(dat$HR), 
                       dat$Parameter,
                      paste0("   ", dat$Parameter))
dat <- dat[,c(1,7,2:6)]

tm <- forest_theme(base_size = 15)
pdf('./01clinical_analysis/output_multicox_result.pdf',width = 20,height = 20)
forest(dat[,c(1:2,6:7)],
       est = dat$HR,
       lower = dat$HR.confit.lower, 
       upper = dat$HR.confit.upper,
       sizes = 1,
       ci_column = 2,
       ref_line = 1,
       # arrow_lab = c("Placebo Better", "Treatment Better"),
       xlim = c(0, 100),
       ticks_at = c(1,50),
       footnote = "",
       theme = tm)
dev.off()

# Signature microbiome -----------------------------------------------------
rm(list = ls())
load('step00_output.Rdata')
library(openxlsx)

micro_sig <- read.xlsx('./01clinical_analysis/source_data/Core_features_coefficients.xlsx',sheet = 2,colNames = FALSE,check.names = F)

microbiome_NR <- read.csv(file = './01clinical_analysis/source_data/Voom-SNM-ACC.csv',header = T,row.names = 1)  %>% .[,micro_sig$X2] %>% set_colnames(str_sub(colnames(.),str_locate(colnames(.),'g__')[,1],str_length(colnames(.))))%>% tibble::rownames_to_column(var = 'names')
microbiome_LR <- read.csv(file = './01clinical_analysis/source_data/Voom-SNM-Filter-Likely-ACC.csv',header = T,row.names = 1)  %>% .[,micro_sig$X2] %>% set_colnames(str_sub(colnames(.),str_locate(colnames(.),'g__')[,1],str_length(colnames(.))))%>% tibble::rownames_to_column(var = 'names')
microbiome_CR <- read.csv(file = './01clinical_analysis/source_data/Voom-SNM-Filter-Plate_Center-ACC.csv',header = T,row.names = 1)  %>% .[,micro_sig$X2] %>% set_colnames(str_sub(colnames(.),str_locate(colnames(.),'g__')[,1],str_length(colnames(.))))%>% tibble::rownames_to_column(var = 'names')
microbiome_PR <- read.csv(file = './01clinical_analysis/source_data/Voom-SNM-Filter-Putative-ACC.csv',header = T,row.names = 1)  %>% .[,micro_sig$X2] %>% set_colnames(str_sub(colnames(.),str_locate(colnames(.),'g__')[,1],str_length(colnames(.))))%>% tibble::rownames_to_column(var = 'names')

clinic.NR <- merge(clinic,microbiome_NR,by='names') 
clinic.LR <- merge(clinic,microbiome_LR,by='names') 
clinic.CR <- merge(clinic,microbiome_CR,by='names')
clinic.PR <- merge(clinic,microbiome_PR,by='names') 


# AUC curve ---------------------------------------------------------------
{
  library(pROC)
  library(rms)
  library(survival)
  library(ggsci)
  library(RColorBrewer)
  }

roc_result <- data.frame(genus=micro_sig$X1,
                         NR=NA,
                         LR=NA,
                         CR=NA,
                         PR=NA)
clinic_type <- list(clinic.NR,clinic.LR,clinic.CR,clinic.PR)
data_type <- c('NR','LR','CR','PR')
for (j in 1:4) {
  roc.list <- roc(OS ~ g__Ilyobacter+g__Isosphaera+g__Singulisphaera+g__Thermopetrobacter+g__Chelativorans+g__Pseudorhodobacter+g__Wenxinia+g__Acidocella+g__Diaphorobacter+g__Vitreoscilla+g__Bacteriovorax+g__Desulfuromonas+g__Proteus+g__Terrimicrobium+g__Haloferula, 
                  smooth=T,
                  data = clinic_type[[j]])
  for(i in 1:15){
    # x=c('Micro subtype, AUC=','Stage, AUC=','Micro-subtype & Stage, AUC=','Nomogram(without MS), AUC=','Nomogram, AUC=')
    # x=paste0(micro_sig[,'X1'],'=')
    # x.1 <- paste0(x[i],round((roc.list[[i]])$auc,2))
    roc_result[i,data_type[j]] <- round((roc.list[[i]])$auc,2)
    # labels <- c(labels,x.1)
  }
}
write.table(roc_result,file = './01clinical_analysis/output_AUC_value_of_each_signature_microbiome.tsv',sep = '\t',row.names = T,col.names = NA)

for (i in 1:4) {
  signature_combine <- psm(Surv(OS.time,OS) ~ g__Ilyobacter+g__Isosphaera+g__Singulisphaera+g__Thermopetrobacter+g__Chelativorans+g__Pseudorhodobacter+g__Wenxinia+g__Acidocella+g__Diaphorobacter+g__Vitreoscilla+g__Bacteriovorax+g__Desulfuromonas+g__Proteus+g__Terrimicrobium+g__Haloferula, 
                           data =  clinic_type[[i]], dist='lognormal') 
  clinic$signature_combine <- predict(signature_combine)[match(clinic$names,clinic_type[[i]]$names)]
  colnames(clinic)[which(colnames(clinic)=='signature_combine')] <- paste0('signature_combine_',data_type[i])
}

roc.list <- roc(OS ~ signature_combine_NR+signature_combine_LR+signature_combine_CR+signature_combine_PR, 
                smooth=T,
                data = clinic)
labels=c()#blank list
for(i in 1:4){
  x=c('15 signature combine(NR), AUC=','15 signature combine(LR), AUC=','15 signature combine(CR), AUC=','15 signature combine(PR), AUC=')
  # x=paste0(micro_sig[,'X1'],'=')
  x.1 <- paste0(x[i],round((roc.list[[i]])$auc,2))
  labels <- c(labels,x.1)
};labels

pdf('./01clinical_analysis/output_auc_four_data_type_signature_combine.pdf',width = 11)
ggroc(roc.list,linetype = 1, size = 1.2) +
  scale_colour_manual(name = "",
                      labels =labels,
                      # values = colorRampPalette(brewer.pal(5,'Paired'))(20),
                      values=pal_npg()(4))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        # panel.background = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  thm+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1),
               color="grey", linetype="dashed")
dev.off()


# Time dependent ROC ------------------------------------------------------
# rm(list = ls())
{ 
  library(survival)
  library(timeROC)
  }
load('step00_output.Rdata')
clinic.roc <- clinic %>% 
  dplyr::mutate(OS.time= round(.$OS.time/30,2)) %>% 
  dplyr::mutate(NR=ifelse(.$Without_combined=='MS2',1,0)) %>% 
  dplyr::mutate(LR=ifelse(.$Likely_combined=='MS2',1,0)) %>% 
  dplyr::mutate(CR=ifelse(.$PC_combined=='MS2',1,0)) %>% 
  dplyr::mutate(PR=ifelse(.$Putative_combined=='MS2',1,0)) %>% 
  dplyr::mutate(stage.new=ifelse(.$Clinical.Stage=='Stage IV',3,
                                 ifelse(.$Clinical.Stage=='Stage III',2,
                                        ifelse(.$Clinical.Stage=='Stage II',1,0)))) %>% 
  dplyr::mutate(T.stage.new=gsub('T','',.$T.stage)) %>% 
  dplyr::mutate(T.stage.new=as.numeric(.$T.stage.new)-1)

# plot
dev.off()
res.cox.stage <- coxph(Surv(OS.time,OS)~T.stage.new,data=clinic.roc)
clinic.roc$cox.stage <- predict(res.cox.stage)
ROC.stage <- timeROC(T=clinic.roc$OS.time,delta=clinic.roc$OS,marker=clinic.roc$cox.stage,iid = T,cause=1,weighting="marginal",times=quantile(clinic.roc$OS.time,probs=seq(0.1,0.9,0.05)))
plotAUCcurve(ROC.stage,col = 'grey')

palette <- pal_npg()(6)
pval <- data.frame()
for (i in 1:4) {
  data_type=c("Without_combined","Likely_combined","PC_combined","Putative_combined")
  clinic.roc <- dplyr::mutate(clinic.roc,MS=clinic.roc[,data_type[i]])
  res.cox <- coxph(Surv(OS.time,OS)~MS+T.stage.new,
                            data=clinic.roc)
  clinic.roc$TYPE <- predict(res.cox)
  ROC <- timeROC(T=clinic.roc$OS.time,delta=clinic.roc$OS,marker=clinic.roc$TYPE,iid = T,cause=1,weighting="marginal",times=quantile(clinic.roc$OS.time,probs=seq(0.1,0.9,0.05)))
  print(plotAUCcurve(ROC,add = T,col = palette[i]))
  outcome <- compare(ROC,ROC.stage,adjusted = T,abseps = 1e-06)
  table <- as.data.frame(outcome$p_values_AUC)
  pval <- rbind(pval,table)
}
pval
abline(v=22.07,lwd=1,lty=2,col="#BC3C29FF")#add a vertical line
abline(v=26.96,lwd=1,lty=2,col="#BC3C29FF")
abline(v=31.894,lwd=1,lty=2,col="#BC3C29FF")
abline(v=98.408,lwd=1,lty=2,col="#BC3C29FF")
legend("bottomleft",c('MS(NR)+Stage','MS(LR)+Stage','MS(CR)+Stage','MS(PR)+Stage','Stage Only'),col=c(pal_npg()(4),'grey'),lty=1,lwd=3)


# Random Forest -----------------------------------------------------------
rm(list = ls())
library(randomForest)
library(tidyverse)
load('step00_output.Rdata')

microbiome_NR <- read.csv(file = './01clinical_analysis/source_data/Voom-SNM-ACC.csv',header = T,row.names = 1)   
microbiome_LR <- read.csv(file = './01clinical_analysis/source_data/Voom-SNM-Filter-Likely-ACC.csv',header = T,row.names = 1)  
microbiome_CR <- read.csv(file = './01clinical_analysis/source_data/Voom-SNM-Filter-Plate_Center-ACC.csv',header = T,row.names = 1)  
microbiome_PR <- read.csv(file = './01clinical_analysis/source_data/Voom-SNM-Filter-Putative-ACC.csv',header = T,row.names = 1)  

clinic <- clinic[match(rownames(microbiome_NR),clinic$names),]
variable <- c('Clinical.Stage',
              'ATYPICAL_MITOTIC_FIGURES',
              'Clinical.Status.3.Mo.Post.Op',
              'Distant.Metastasis',
              'Excessive.Hormone',
              'Surgical.margin',
              'PRIMARY_THERAPY_OUTCOME',
              'Neoplasm.Status')
result <- data.frame(parameter=variable,NR=NA,LR=NA,CR=NA,PR=NA)
for (i in 1:4) {
  dat <- list(microbiome_NR,microbiome_LR,microbiome_CR,microbiome_PR)[[i]]
  for (j in 1:length(variable)) {
    mat <-  dat  %>% dplyr::mutate(group=as.factor(clinic[,variable[j]])) %>% drop_na(group)
    set.seed(430)
    rf = randomForest(group~ .,data=mat, ntree=500, importance=TRUE, proximity=TRUE)
    res1 <- as.data.frame(rf$importance) %>% .[order(.$MeanDecreaseAccuracy,decreasing = T),]
    res2 <- as.data.frame(rf$importance) %>% .[order(.$MeanDecreaseGini,decreasing = T),]
    x <- rownames(res1)[1:10] %>% str_sub(.,str_locate(.,'g__')[,1],str_length(.))
    y <- rownames(res2)[1:10]%>% str_sub(.,str_locate(.,'g__')[,1],str_length(.))
    result[j,i+1] <- paste(intersect(x,y),collapse=";")
  }
}
write.table(result,file = './01clinical_analysis/output_random_forest.tsv',sep = '\t',row.names = T,col.names = NA)


# MS & PHENOTYPE ----------------------------------------------------------
rm(list = ls())
getwd()

{
  library(plyr)
  library(dplyr)
  library(ggplot2)
  library(ggthemes)
  library(ggpubr)
  library(RColorBrewer)
  library(tidyverse)
}
detach("package:randomForest", unload=TRUE) #exclude the influence of the package to next plotting
load('step00_output.Rdata')

clinic$mRNA_K4_new <- ifelse(str_detect(clinic$mRNA_K4,'proliferation'),NA,clinic$mRNA_K4)
pheno <- c('Histology','C1A_C1B','mRNA_K4_new','miRNA_cluster','MethyLevel','SCNA_cluster','protein_cluster','COC')
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
        legend.position='bottom',
        legend.text = element_text(size = 15,hjust = 0),
        legend.title = element_text(size = 20)
  )

# plot for phenos
data_type <- c("Without_combined","Likely_combined","PC_combined","Putative_combined")
for (j in 1:4) {
  clinic$MS <- clinic[,data_type[j]]
  for (i in 1:length(pheno)) {
    dat <- data.frame(table(clinic[,c('MS',pheno[i])])) %>% 
      ddply(.(MS),transform,percent=Freq/sum(Freq)*100) 
    colnames(dat)[1:2] <- c('Var1','Var2')
    x <- matrix(dat$Freq,ncol=2)
    if (sum(x<5)>0){
      pvalue <- fisher.test(x)$p.value
      print(pvalue)
    }else{
      pvalue <- chisq.test(x)$p.value
      print(pvalue)
    }
    legend_n <- length(unique(dat$Var2))
    pdf(paste0('./01clinical_analysis/output_MS&',pheno[i],'_',data_type[j],'.pdf'))
    p <- dat%>%
      drop_na() %>% 
      ggplot(aes(fill=Var2,x=Var1,y=Freq))+
      geom_bar(position = 'fill',stat = 'identity',width = 0.8)+
      guides(fill=guide_legend(title=pheno[i],direction = "horizontal",
                               nrow=ifelse(legend_n>4,3,ifelse(legend_n>1,2,1)))) +
      # scale_fill_brewer(palette = 'Set1')+
      # scale_fill_nejm(alpha = 0.7)+
      scale_fill_manual(values = c(pal_nejm(alpha = 0.9)(3)[2:3],pal_nejm(alpha = 0.9)(1),pal_nejm(alpha = 0.9)(8)[4:8]))+
      scale_y_continuous(labels = scales::percent)+
      coord_fixed(ratio = 5/2)+
      labs(x='',y='Proportion %',
           subtitle = paste0("P ", ifelse(pvalue<0.001, "< 0.001", paste0("= ",round(pvalue,3)))))+
      thm
    print(p)
    dev.off()
  }
}


# MS & ADS
for (i in 1:4) {
  clinic$MS <- clinic[,data_type[i]]
  pdf(paste0('./01clinical_analysis/output_MS&ADS','_',data_type[i],'.pdf'))
  p <- ggplot(data = clinic,aes(x=MS, y=ADS)) + 
    geom_violin(trim=T) +
    geom_jitter(shape=16, color="grey",size=2.0,position=position_jitter(0.2))+
    geom_boxplot(width=0.1,position=position_dodge(0.8),aes(fill=MS))+
    thm+
    theme(legend.position = 'right')+
    scale_fill_manual(values = c(MS1="#0072B5CC",MS2="#E18727CC"))+
    guides(fill=guide_legend(direction = "horizontal",ncol =1))+
    labs(x="", y = "ADS ", fill = "")+  
    stat_compare_means(aes(x = MS , y = ADS),label = "p.format",size=6)
  print(p)
  dev.off()
}
