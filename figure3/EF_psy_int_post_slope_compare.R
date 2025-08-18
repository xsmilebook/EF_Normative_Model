rm(list = ls())
library(ggplot2)
library(tidyverse)
library(mgcv)
library(readxl)
library(parallel)
library(gamlss)
library(scales)
library(tableone)
library(openxlsx)
library(showtext)
# 文件路径设置
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/EFdevelopment/Crsen_survey_GD_YunFu/data/EF_results'
  demopath <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/EFdevelopment/Crsen_survey_GD_YunFu/data/rawdata_results0616'
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/EFdevelopment/Crsen_survey_GD_YunFu/interfileFolder"
  functionFolder <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/EFdevelopment/Crsen_survey_GD_YunFu/Rcode_EFnorms/functions"
  
} else {
  datapath <- '/Users/tanlirou/Documents/YF_EF_psy/EF_psy/correlation/data'
  FigureFolder <- '/Users/tanlirou/Documents/YF_EF_psy/EF_psy/code_pure_202507/results/Figure3_int_delet_Z'
  interfileFolder <- "/Users/tanlirou/Documents/YF_EF_psy/EF_psy/correlation/data"
  functionFolder <- "/Users/tanlirou/Documents/YF_EF_psy/EF_psy/code_pure_202507/functions"
  resultFolder <- "/Users/tanlirou/Documents/YF_EF_psy/EF_psy/code_pure_202507/results/Figure3_int_delet_Z"
}

beta_table <- read.csv("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/code_pure_202507/results/Figure2_corr_delete_Z/corr_results_with_bonf.csv")
# 读取数据
GNGd_data <- readRDS(paste0(datapath, '/GNGd_prime.deviations.rds'))
back1_data <- readRDS(paste0(datapath, '/back1Acc.deviations.rds'))
back2_data <- readRDS(paste0(datapath, '/back2Acc.deviations.rds'))
GNGd_data$Sex <- as.factor(GNGd_data$Sex)
back1_data$Sex <- as.factor(back1_data$Sex)
back2_data$Sex <- as.factor(back2_data$Sex)
#source("~/Documents/EF_yunfu_check/code/functions/gamsmooth.R")
source(paste0(functionFolder, '/gam_varyingcoefficients.R'))
# 定义变量
psyc_variables_continous <- c("SDQ_PB_sum", "SDQ_H_sum", "SDQ_CP_sum", "SDQ_PP_sum", "SDQ_ES_sum")
EFvars.set <- matrix(c("d_prime_deviationZ", "GNGd",
                       "Oneback_acc_deviationZ", "back1",
                       "Twoback_acc_deviationZ", "back2"), byrow=TRUE, ncol=2, dimnames=list(NULL, c("varname", "dataname")))
EFvars.set <- as.data.frame(EFvars.set)
#####
# GNGd_data[, paste0(psyc_variables_continous, "_z")] <- scale(GNGd_data[, psyc_variables_continous])
# back1_data[, paste0(psyc_variables_continous, "_z")] <- scale(back1_data[, psyc_variables_continous])
# back2_data[, paste0(psyc_variables_continous, "_z")] <- scale(back2_data[, psyc_variables_continous])
standardize_clean <- function(df, vars) {
  for (var in vars) {
    x <- df[[var]]
    x <- x[!is.na(x)]
    mu <- mean(x)
    sd_val <- sd(x)
    valid_index <- which(abs(df[[var]] - mu) <= 3 * sd_val)
    z_varname <- paste0(var, "_z")
    df[[z_varname]] <- NA
    df[[z_varname]][valid_index] <- scale(df[[var]][valid_index])
  }
  return(df)
}

original_vars <- c("SDQ_PB_sum", "SDQ_H_sum", "SDQ_CP_sum", "SDQ_PP_sum", "SDQ_ES_sum")

GNGd_data  <- standardize_clean(GNGd_data,  original_vars)
back1_data <- standardize_clean(back1_data, original_vars)
back2_data <- standardize_clean(back2_data, original_vars)

psyc_variables_continous <- paste0(original_vars, "_z")

# 设置参数
knots <- 3
set_fx <- FALSE 
increments <- 1000
draws <- 1000
return_posterior_coefficients <- T


# 初始化结果存储列表
interaction_results <- list()

# 主循环 - 循环计算所有模型并保存结果
for (i in 1:nrow(EFvars.set)) {
  int_var <- EFvars.set$varname[i]
  dataname <- paste0(EFvars.set$dataname[i], "_data")
  
  for (j in 1:length(psyc_variables_continous)) {
    dependentvar <- psyc_variables_continous[j]
    smooth_var <- "Age_year"
    
    cat(paste("Processing:", dataname, "~", int_var, "by", smooth_var, "for dependent var:", dependentvar, "\n"))
    
    # 调用函数，返回后验斜率数据
    result <- gam.varyingcoefficients(
      dependentvar = dependentvar,
      dataname = dataname,
      smooth_var = smooth_var,
      int_var = int_var,
      covariates = "Sex",
      knots = knots,
      set_fx = set_fx,
      increments = increments,
      draws = draws,
      return_posterior_coefficients = return_posterior_coefficients
    )
    
    # 保存结果
    model_name <- paste(EFvars.set$dataname[i], psyc_variables_continous[j], sep = "_")
    interaction_results[[model_name]] <- result
  }
}

# 保存所有结果为RDS文件
saveRDS(interaction_results, file = paste0(resultFolder, "/all_interaction_results.rds"))


interaction_results <- readRDS(paste0(resultFolder, "/all_interaction_results.rds"))
task_mapping <- c("GNGd" = "Go/No-Go", "back1" = "1-back", "back2" = "2-back")
variable_mapping <- c(
  "SDQ_ES_sum_z" = "Emotional Symptoms", "SDQ_PP_sum_z" = "Peer Problems",
  "SDQ_CP_sum_z" = "Conduct Problems", "SDQ_H_sum_z" = "Hyperactivity",
  "SDQ_PB_sum_z" = "Prosocial Behavior"
)

# 定义统一的主题
custom_theme <- theme_minimal() +
  theme(
    axis.line = element_line(size = 0.25, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 9),
    plot.title = element_text(size = 9, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.25),
    axis.ticks.length = unit(0.05, "cm")
  )

# 遍历所有模型并生成图表
for (model_name in names(interaction_results)) {
  result <- interaction_results[[model_name]]
  
  # 从模型名称中提取任务和变量信息
  parts <- strsplit(model_name, "_")[[1]]
  current_task <- parts[1]
  current_variable_raw <- paste(parts[2:length(parts)], collapse = "_")
  
  # 获取绘图所需的数据
  current_slope_data <- result[[2]]
  
  # 获取总体beta值
  slope_value <- beta_table %>%
    filter(dataname == current_task, 
           parcel == current_variable_raw) %>% # <--- 修正：使用带"_z"的变量名进行匹配
    pull(beta) # <--- 修正：列名是beta，不是slope
  
  # 如果找不到beta值，跳过
  if (length(slope_value) == 0) {
    message(paste("找不到", current_task, "和", current_variable_raw, "的beta值，跳过绘图."))
    next
  }
  
  # 计算每个年龄点的斜率中位数、95% CI，并进行显著性检验
  slope_summary <- current_slope_data %>%
    group_by(Age_year) %>%
    summarise(
      median_slope = median(slope, na.rm = TRUE),
      lower_95CI = quantile(slope, probs = 0.025, na.rm = TRUE),
      upper_95CI = quantile(slope, probs = 0.975, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      is_significant = (slope_value < lower_95CI | slope_value > upper_95CI),
      abs_diff = abs(median_slope - slope_value)
    )
  
  # 绘制主图：斜率随年龄的变化趋势
  slope_plot_main <- ggplot(slope_summary, aes(x = Age_year)) +
    geom_line(aes(y = median_slope), color = "black", size = 0.5) +
    geom_ribbon(aes(ymin = lower_95CI, ymax = upper_95CI), alpha = 0.3, fill = "grey80") +
    geom_hline(yintercept = slope_value, linetype = "dashed", color = "red", size = 0.4) +
    labs(
      title = paste0(variable_mapping[current_variable_raw], " ~ ", task_mapping[current_task]),
      x = "",
      y = "Slope"
    ) +
    scale_x_continuous(name = "", limits = c(11, 18), breaks = seq(11, 18, by = 1)) +
    custom_theme
  
  # 绘制显著性条形图（子图）
  slope_plot_bar <- ggplot(slope_summary, aes(x = Age_year, y = 1)) +
    geom_col(aes(fill = abs_diff), data = filter(slope_summary, is_significant == TRUE),
             width = 1, color = "transparent") +
    geom_col(data = filter(slope_summary, is_significant == FALSE), fill = "white",
             width = 1, color = "transparent") +
    scale_fill_gradient(low = "#d9e0f0", high = "#A4C5DF", na.value = "white", guide = "none") +
    scale_x_continuous(name = "Age_year", limits = c(11, 18), breaks = seq(11, 18, by = 1)) +
    theme_minimal() +
    theme(
      axis.title.y = element_blank(), axis.text.y = element_blank(),
      axis.ticks.y = element_blank(), axis.line.y = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(t = -10, r = 0, b = 0, l = 0, unit = "pt")
    )
  
  # 合并主图和子图并保存
  final_plot <- plot_grid(slope_plot_main, slope_plot_bar, ncol = 1, rel_heights = c(0.8, 0.2), align = "v")
  
  ggsave(
    filename = paste0(FigureFolder, "/slope_comparison_", current_task, "_", current_variable_raw, ".pdf"),
    plot = final_plot,
    width = 6, height = 6.5, units = "cm"
  )
}
