# 清理环境
rm(list=ls())

library(readxl)
library(tidyverse)
library(mgcv)
library(psych)
library(gamlss)
library(scales)
library(tableone)
library(openxlsx)
library(ggplot2)
library(patchwork) 

# 设置路径
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- '/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psy/interfileFolder'
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psy/interfileFolder_back12before"
  functionFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psycode/Rcode/functions"
  resultFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psy/results_corr"
  FigureFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psy/results_corr"
} else {
  datapath <- 'D:/datasets/yunfu/raw_data'
  FigureFolder <- 'D:/datasets/yunfu/figures/fig1'
  interfileFolder <- "D:/datasets/yunfu/interfile_folder/temp"
  functionFolder <- "D:/code/EF_Normative_Model/functions"
  resultFolder <- "D:/datasets/yunfu/results/psy_corr"
}
source(paste0(functionFolder, "/ordinalcorr_new.R"))

# 读取数据
PHQ_data <- read_xlsx(paste0(interfileFolder, '/PHQ.xlsx'))
GNGd_data <- read_rds(paste0(interfileFolder, '/GNGd_prime.deviations.rds'))
back1_data <- read_rds(paste0(interfileFolder, '/back1Acc.deviations.rds'))
back2_data <- read_rds(paste0(interfileFolder, '/back2Acc.deviations.rds'))

# 数据预处理
PHQ_data <- PHQ_data %>%
  select(用户ID, PHQ_sum) %>%  
  mutate(用户ID = as.character(用户ID))

GNGd_data <- GNGd_data %>% mutate(x__ID = as.character(x__ID))
back1_data <- back1_data %>% mutate(x__ID = as.character(x__ID))
back2_data <- back2_data %>% mutate(x__ID = as.character(x__ID))

# 找到每个数据集中以deviationZ结尾的列名
gngd_dev_col <- names(GNGd_data)[str_ends(names(GNGd_data), "deviationZ")]
back1_dev_col <- names(back1_data)[str_ends(names(back1_data), "deviationZ")]
back2_dev_col <- names(back2_data)[str_ends(names(back2_data), "deviationZ")]

# 创建用于分析的数据框
# GNGd数据
gngd_analysis_data <- GNGd_data %>%
  select(x__ID, all_of(gngd_dev_col)) %>%
  rename(用户ID = x__ID) %>%
  inner_join(PHQ_data, by = "用户ID")

# back1数据
back1_analysis_data <- back1_data %>%
  select(x__ID, all_of(back1_dev_col)) %>%
  rename(用户ID = x__ID) %>%
  inner_join(PHQ_data, by = "用户ID")

# back2数据
back2_analysis_data <- back2_data %>%
  select(x__ID, all_of(back2_dev_col)) %>%
  rename(用户ID = x__ID) %>%
  inner_join(PHQ_data, by = "用户ID")

# 检查数据
print(paste("GNGd分析数据行数:", nrow(gngd_analysis_data)))
print(paste("back1分析数据行数:", nrow(back1_analysis_data)))
print(paste("back2分析数据行数:", nrow(back2_analysis_data)))

# 进行相关分析
# 注意：需要根据你的实际情况调整协变量、平滑变量等参数
# 这里假设使用年龄作为协变量和平滑变量，knots=5

# GNGd相关分析
gngd_results <- ordinalcorr(
  dependentvar = "PHQ_sum",
  dataname = "gngd_analysis_data",
  interest.indep.var = gngd_dev_col,
  covariates = "",  # 如果有协变量，填写如 "age + gender"
  smoothvar = "PHQ_sum",  # 如果需要控制非线性趋势的变量
  knots = 5,
  stats_only = TRUE
)

# back1相关分析
back1_results <- ordinalcorr(
  dependentvar = "PHQ_sum",
  dataname = "back1_analysis_data",
  interest.indep.var = back1_dev_col,
  covariates = "",  # 如果有协变量，填写如 "age + gender"
  smoothvar = "PHQ_sum",  # 如果需要控制非线性趋势的变量
  knots = 5,
  stats_only = TRUE
)

# back2相关分析
back2_results <- ordinalcorr(
  dependentvar = "PHQ_sum",
  dataname = "back2_analysis_data",
  interest.indep.var = back2_dev_col,
  covariates = "",  # 如果有协变量，填写如 "age + gender"
  smoothvar = "PHQ_sum",  # 如果需要控制非线性趋势的变量
  knots = 5,
  stats_only = TRUE
)

# 整理结果
all_results <- rbind(
  cbind(gngd_results, dataset = "GNGd"),
  cbind(back1_results, dataset = "back1"),
  cbind(back2_results, dataset = "back2")
)

# 转换为数据框并保存
results_df <- as.data.frame(all_results)
results_df <- results_df %>%
  mutate(
    across(everything(), as.character),
    gam.smooth.t = as.numeric(gam.smooth.t),
    gam.smooth.pvalue = as.numeric(gam.smooth.pvalue),
    anova.cov.pvalue = as.numeric(anova.cov.pvalue),
    partialRsq = as.numeric(partialRsq),
    correstimate = as.numeric(correstimate),
    samplesize = as.numeric(samplesize),
    beta = as.numeric(beta)
  )

# 查看结果
print(results_df)

# 保存结果
write.xlsx(results_df, paste0(resultFolder, "/PHQ_deviation_correlations.xlsx"))

# 打印摘要信息
cat("\n=== 相关分析结果摘要 ===\n")
for(i in 1:nrow(results_df)) {
  cat(sprintf("数据集: %s\n", results_df$dataset[i]))
  cat(sprintf("变量: %s\n", results_df$interest.indep.var[i]))
  cat(sprintf("相关系数: %.3f\n", results_df$correstimate[i]))
  cat(sprintf("p值: %.3f\n", results_df$anova.cov.pvalue[i]))
  cat(sprintf("偏R²: %.3f\n", results_df$partialRsq[i]))
  cat("---\n")
}