rm(list = ls())
getwd()

### load requirements
{
  library(dplyr)
  library(tidyr)
  library(grid)
  library(forestploter)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(magrittr)
}

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



### 1. Clustering and clinical parameters
variable <- c('SEX',
              'AGE_new',
              'RACE',
              'STAGE',
              'T_STAGE',
              'N_STAGE',
              'M_STAGE',
              'ATYPICAL_MITOTIC_FIGURES',
              'CAPSULAR_INVASION',
              'CLINICAL_STATUS_WITHIN_3_MTHS_SURGERY',
              'CYTOPLASM_PRESENCE_LESS_THAN_EQUAL_25_PERCENT',
              'DIFFUSE_ARCHITECTURE',
              'HISTORY_ADRENAL_HORMONE_EXCESS',
              'LATERALITY',
              'MITOTIC_RATE',
              'NECROSIS',
              'NUCLEAR_GRADE_III_IV',
              'PHARMACEUTICAL_TX_ADJUVANT',
              'PHARM_TX_MITOTANE_ADJUVANT',
              'PHARM_TX_MITOTANE_THERAPUTIC_AT_REC',
              'RADIATION_TREATMENT_ADJUVANT',
              'RESIDUAL_TUMOR',
              'SINUSOID_INVASION',
              'TREATMENT_OUTCOME_FIRST_COURSE',
              'TUMOR_STATUS',
              'WEISS_VENOUS_INVASION',
              'TUMOR_WEIGHT_new')
clinic <- df_cluster_clinical %>% 
  dplyr::mutate(HISTORY_ADRENAL_HORMONE_EXCESS=ifelse(HISTORY_ADRENAL_HORMONE_EXCESS=='None','Absent','Present')) %>% 
  dplyr::mutate(AGE_new = ifelse(AGE>median(AGE),'High','Low')) %>% 
  dplyr::mutate(TUMOR_WEIGHT_new = ifelse(TUMOR_WEIGHT>median(TUMOR_WEIGHT),'High','Low')) 

# 计算卡方检验，并整理结果
result_list <- list()
for (var in variable) {  # 排除 MS 列
  
  # 创建交叉表
  tbl <- table(clinic$MS, clinic[[var]])  
  
  # 检查是否使用Fisher检验
  use_fisher <- any(tbl < 5)  # 如果有频数小于 5，则使用 Fisher
  test <- if (use_fisher) fisher.test(tbl) else chisq.test(tbl)
  
  # 获取临床变量和类别
  categories <- rownames(tbl)
  groups <- colnames(tbl)
  
  # 创建结果数据框
  result_df <- data.frame(
    Variable = rep(var, length(categories) * length(groups)),
    Category = rep(categories, each = length(groups)),
    Group = rep(groups, length(categories)),
    Count = as.vector(tbl),
    Test_Type = ifelse(use_fisher, "Fisher", "Chi-Square"),  # 统计方法
    Statistic = ifelse(use_fisher, NA, test$statistic),  # 统计量
    P_Value = round(test$p.value,2)  # P 值
  )
  
  # 将当前变量的结果添加到结果列表中
  result_list[[var]] <- result_df
}
result <- bind_rows(result_list) # 合并所有结果为一个数据框
result_wide <- result %>%
  pivot_wider(names_from = Category, values_from = Count) %>%
  dplyr::select(Variable, Group, MS1, MS2, Test_Type, Statistic, P_Value) %>% 
  dplyr::mutate(P_Value=ifelse(P_Value<0.05,'*',P_Value))

# write.table(result_wide,file = './cluster_and_clinical_parameters/ACC_params_MS(chi-squared test).tsv',sep = '\t',row.names = F,col.names = TRUE)


### 2. Univarite and Multivariate COX 
clinic <- df_cluster_clinical %>% 
  dplyr::mutate(HISTORY_ADRENAL_HORMONE_EXCESS=ifelse(HISTORY_ADRENAL_HORMONE_EXCESS=='None','Absent','Present')) %>% 
  dplyr::mutate(AGE_new = ifelse(AGE>median(AGE),'High','Low')) %>% 
  dplyr::mutate(TUMOR_WEIGHT_new = ifelse(TUMOR_WEIGHT>median(TUMOR_WEIGHT),'High','Low')) 

# 单因素Cox回归
univ_cox_results <- lapply(c("MS", "AGE", "T_STAGE",'N_STAGE','M_STAGE','SINUSOID_INVASION'), function(var) {
  formula <- as.formula(paste("Surv(OS_Months, OS_Status) ~", var))
  cox_model <- coxph(formula, data = clinic)
  summary_model <- summary(cox_model)
  
  # 提取HR和95%置信区间
  result <- data.frame(
    Variable = var,
    HR = exp(coef(cox_model)),
    Lower_95 = exp(coef(cox_model) - 1.96 * summary_model$coefficients[, "se(coef)"]),
    Upper_95 = exp(coef(cox_model) + 1.96 * summary_model$coefficients[, "se(coef)"]),
    P_Value = summary_model$sctest[3],
    Model = "Univariate"
  )
  return(result)
})
univ_cox_results <- do.call(rbind, univ_cox_results) %>%  # 将单因素结果合并为一个数据框
  dplyr::mutate(Variable=rownames(.))

# 多因素Cox回归
multi_cox <- coxph(formula = Surv(OS_Months,OS_Status)~ MS+AGE+T_STAGE+N_STAGE+M_STAGE+SINUSOID_INVASION, data = clinic)
multi_cox_summary <- summary(multi_cox)
multi_cox_results <- data.frame(
  HR = exp(coef(multi_cox)),
  Lower_95 = exp(coef(multi_cox) - 1.96 * multi_cox_summary$coefficients[, "se(coef)"]),
  Upper_95 = exp(coef(multi_cox) + 1.96 * multi_cox_summary$coefficients[, "se(coef)"]),
  P_Value = signif(multi_cox_summary$coefficients[,"Pr(>|z|)"],digits = 3),
  Model = "Multivariate"
) %>% 
  dplyr::mutate(Variable=rownames(.)) %>% 
  .[,colnames(univ_cox_results)]

# 合并单因素和多因素Cox回归结果
cox_results <- rbind(univ_cox_results, multi_cox_results)
# write.table(cox_results,file = './cluster_and_clinical_parameters/unicox_multicox_result.tsv',sep = '\t',row.names = F,col.names = T)

# 绘制森林图
cox_results$Variable <- factor(cox_results$Variable, levels = rev(c("MSMS2",'AGE','T_STAGET2','T_STAGET3','T_STAGET4','N_STAGEN1','M_STAGEM1','SINUSOID_INVASIONSinusoid Invasion Present')))
cox_results$HR_CI <- paste0(round(cox_results$HR, 2), " (", 
                            round(cox_results$Lower_95, 2), "-", 
                            round(cox_results$Upper_95, 2), ")")
# pdf('./cluster_and_clinical_parameters/unicox_multicox_result.pdf',width = 14,height = 6)
p_value_to_stars <- function(p) {
  if (is.na(p)) {
    return("")
  } else if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}
ggplot(cox_results, aes(x = HR, y = Variable, xmin = Lower_95, xmax = Upper_95, color = Model)) +
  geom_point(size = 8, shape = 15) +  # 点的大小和形状
  geom_errorbarh(height = 0.4) +  # 水平误差条
  scale_x_continuous(trans = 'log10') +  # 对数坐标轴
  scale_color_manual(values = c("#0072B5CC", "#E18727CC")) +  # 设置颜色
  theme_minimal(base_size = 18) +  # 使用简洁主题并调整基本字体大小
  theme(legend.position = "top", 
        legend.title = element_blank(), 
        strip.background = element_blank(),  # 移除分面标题背景
        panel.spacing = unit(1, "lines")) +  # 控制面板间距
  labs(
    x = "Hazard Ratio (HR)",
    y = "Variables",
    title = "Forest Plot for Univariate and Multivariate Cox Regression"
  ) + 
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +  # HR=1的虚线
  facet_grid(cols = vars(Model), scales = "free_y") +  # 按照单因素和多因素分面显示
  geom_text(aes(x = max(cox_results$Upper_95) * 1.1, y = Variable, label = sapply(P_Value, p_value_to_stars)), 
            size = 6, hjust = 0)   # p值显示在图像右侧竖直线
# dev.off()



