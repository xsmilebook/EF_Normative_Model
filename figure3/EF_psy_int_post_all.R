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
  FigureFolder <- '/Users/tanlirou/Documents/YF_EF_psy/EF_psy/correlation_change/resutls250515'
  interfileFolder <- "/Users/tanlirou/Documents/YF_EF_psy/EF_psy/correlation/data"
  functionFolder <- "/Users/tanlirou/Documents/YF_EF_psy/EF_psy/correlation_change/functions"
  resultFolder <- "/Users/tanlirou/Documents/YF_EF_psy/EF_psy/correlation_change/resutls250515"
}

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
###3# 对 SDQ 连续变量进行标准化：mean = 0, sd = 1
GNGd_data[, paste0(psyc_variables_continous, "_z")] <- scale(GNGd_data[, psyc_variables_continous])
back1_data[, paste0(psyc_variables_continous, "_z")] <- scale(back1_data[, psyc_variables_continous])
back2_data[, paste0(psyc_variables_continous, "_z")] <- scale(back2_data[, psyc_variables_continous])

psyc_variables_continous <- paste0(psyc_variables_continous, "_z")
# 设置参数
knots <- 3
set_fx <- FALSE 
increments <- 1000
draws <- 1000
return_posterior_coefficients <- T

# 初始化存储交互效应和 p 值的列表
int.results.df <- data.frame(
  dependentvar = character(),
  Age = numeric(),
  parcel = character(),
  dataname = character(),
  Interaction = numeric(),
  anova.cov.pvalue = numeric(),
  anova.pvalue = numeric(),
  anovap.fdr = numeric(),
  sig = logical(),
  stringsAsFactors = FALSE
)

# 初始化存储显著交互效应的列表
significant_interactions <- data.frame(dependentvar = character(),
                                       Age = numeric(),
                                       parcel = character(),
                                       dataname = character(),
                                       Interaction = numeric(),
                                       slope_data = I(list()),
                                       stringsAsFactors = FALSE)
nosignificant_interactions <- data.frame(dependentvar = character(),
                                       Age = numeric(),
                                       parcel = character(),
                                       dataname = character(),
                                       Interaction = numeric(),
                                       slope_data = I(list()),
                                       stringsAsFactors = FALSE)

# 主循环
for (i in 1:nrow(EFvars.set)) {
  int_var <- EFvars.set$varname[i]
  dataname <- paste0(EFvars.set$dataname[i], "_data")
  p_values_task <- c()
  anova.pvalues_task <- c()
  Interactions_task <- c()
  slope_data_list <- list()
  
  for (j in 1:length(psyc_variables_continous)) {
    dependentvar <- psyc_variables_continous[j]
    smooth_var <- "Age_year"
    result <- gam.varyingcoefficients(dependentvar, dataname, smooth_var, int_var, covariates = "Sex", knots, set_fx, increments, draws, return_posterior_coefficients)
    
    # 提取交互效应的 p 值和估计值
    p_value <- as.numeric(result[[1]][1, "anova.int.pvalue"])
    Interaction <- as.numeric(result[[1]][1, "IntpartialRsq"])
    anova.pvalues <- as.numeric(result[[1]][1, "anova.pvalues"])
    slope_data <- result[[2]] 
    slope_data_list[[j]] <- data.frame(Age_year = slope_data[[smooth_var]], slope = slope_data$slope)
    
    # 将值存入列表
    p_values_task <- c(p_values_task, p_value)
    anova.pvalues_task <- c(anova.pvalues_task, anova.pvalues)
    Interactions_task <- c(Interactions_task, Interaction)
    
    
    temp_df <- data.frame(
      dependentvar = dependentvar,
      Age = smooth_var,
      parcel = int_var,
      dataname = EFvars.set$dataname[i],
      Interaction = Interaction,
      anova.cov.pvalue = p_value,
      anova.pvalues = anova.pvalues,
      anovap.fdr = NA,
      sig = FALSE,
      stringsAsFactors = FALSE
    )
    
    int.results.df <- rbind(int.results.df, temp_df)
  }
  
  # FDR 校正并更新显著性标记
  fdr_corrected_p <- p.adjust(anova.pvalues_task, method = "fdr")
  int.results.df[int.results.df$parcel == int_var, "anovap.fdr"] <- fdr_corrected_p
  int.results.df[int.results.df$parcel == int_var, "sig"] <- (fdr_corrected_p < 0.05)
  
  sig_indices <- which(fdr_corrected_p < 0.05)
  if (length(sig_indices) > 0) {
    for (idx in sig_indices) {
      significant_interactions <- rbind(significant_interactions, 
                                        data.frame(dependentvar = psyc_variables_continous[idx], 
                                                   Age = smooth_var,
                                                   parcel = int_var, 
                                                   dataname = EFvars.set$dataname[i], 
                                                   Interaction = Interactions_task[idx],
                                                   slope_data = I(list(slope_data_list[[idx]]))))
    }
  }
  nosig_indices <- which(fdr_corrected_p > 0.05)
  if (length(nosig_indices) > 0) {
    for (idx in nosig_indices) {
      nosignificant_interactions <- rbind(nosignificant_interactions, 
                                        data.frame(dependentvar = psyc_variables_continous[idx], 
                                                   Age = smooth_var,
                                                   parcel = int_var, 
                                                   dataname = EFvars.set$dataname[i], 
                                                   Interaction = Interactions_task[idx],
                                                   slope_data = I(list(slope_data_list[[idx]]))))
    }
  }
  write.xlsx(int.results.df, file = paste0(resultFolder, "/interaction_results_full.xlsx"), row.names = FALSE)
}
write_rds(significant_interactions, file = paste0(resultFolder, "/int_significant_YF.rds"))
write_rds(nosignificant_interactions, file = paste0(resultFolder, "/int_Nosignificant_YF.rds"))

#int.results.df <- read_xlsx(paste0(resultFolder, "/interaction_results_full.xlsx"))

# 绘制热图
lwth <- min(int.results.df$Interaction, na.rm = TRUE)
upth <- max(int.results.df$Interaction, na.rm = TRUE)
y_levels <- c("SDQ_PB_sum_z", "SDQ_H_sum_z", "SDQ_CP_sum_z", "SDQ_PP_sum_z", "SDQ_ES_sum_z")

for (i in 1:nrow(EFvars.set)) {
  EFvar.tmp <- EFvars.set$varname[i]
  dataname <- EFvars.set$dataname[i]
  corr.result.tmp <- int.results.df[which(int.results.df$parcel == EFvar.tmp & int.results.df$dataname == dataname), ]
  corr.result.tmp$sig <- as.logical(corr.result.tmp$sig)  # 转换 sig 列为逻辑值
  corr.result.tmp.sig <- corr.result.tmp[corr.result.tmp$sig == TRUE, ]

  # 绘制热图
  Fig <- ggplot() +
    geom_tile(data = corr.result.tmp, aes(x = parcel, y = dependentvar, fill = Interaction), color = "white") +
    geom_text(data = corr.result.tmp.sig, aes(x = parcel, y = dependentvar, label = "*"), vjust = 0.75, hjust = 0.5, size = 9) +
    scale_fill_distiller(type = "seq", palette = "Reds", limits = c(lwth, upth), direction = 1) +
    scale_y_discrete(limits = y_levels,
                     labels = c("SDQ_PB_sum_z" = "Prosocial Behavior","SDQ_H_sum_z" = "Hyperactivity",
                                "SDQ_CP_sum_z" = "Conduct Problems","SDQ_PP_sum_z" = "Peer Problems",
                                "SDQ_ES_sum_z" = "Emotional Symptoms","SDQ_sum_z" = "Total SDQ Score")) +
    scale_x_discrete(labels = c("d_prime_deviationZ" = "Go/No-go",
                                "Oneback_acc_deviationZ" = "1-back",
                                "Twoback_acc_deviationZ" = "2-back")) +
    labs(title = paste0("Interaction Effects in ", dataname),
         x = "Tasks", y = "Mental Health") +
    theme(axis.line = element_blank(),
          aspect.ratio = 1.2,
          axis.text.x = element_text(size = 12, hjust = 0.5),
          axis.text.y = element_text(size = 12, hjust = 0.5, vjust = 0.5),
          axis.title = element_text(size = 18),
          plot.title = element_text(size = 18, hjust = 0.5),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          panel.background = element_rect(fill = NA),
          panel.grid.major = element_line(linewidth = 0),
          panel.grid.minor = element_line(linewidth = 1))

  print(Fig)
  ggsave(paste0(FigureFolder, "/corr_int_", EFvar.tmp, "_tasks.pdf"), plot = Fig, width = 16, height = 12, units = "cm")
}


# 确定颜色的上下限
lwth <- min(int.results.df$Interaction, na.rm = TRUE)
upth <- max(int.results.df$Interaction, na.rm = TRUE)

# 自定义量表顺序
y_levels <- c("SDQ_sum_z", "SDQ_ES_sum_z", "SDQ_PP_sum_z", "SDQ_CP_sum_z", "SDQ_H_sum_z", "SDQ_PB_sum_z")  # 量表顺序

# 创建整合的任务数据框
combined_data <- int.results.df
combined_data$sig <- as.logical(combined_data$sig)  # 转换 sig 列为逻辑值
combined_data$parcel <- factor(combined_data$parcel, levels = c( "d_prime_deviationZ", "Oneback_acc_deviationZ", "Twoback_acc_deviationZ"), 
                               labels = c("Go/No-go ", "1-back", "2-back"))  # 设置任务的顺序和标签

# 定义自定义颜色梯度
custom_colors <- c( "#f5dfdb", "#edb8b0", "#e69191", "#c25759")  # 你可以替换成你喜欢的颜色

# 绘制整合的热图
Fig <- ggplot() +
  geom_tile(data = combined_data, aes(x = parcel, y = dependentvar, fill = Interaction), color = "white") +  
  geom_text(data = combined_data[combined_data$sig == TRUE, ], 
            aes(x = parcel, y = dependentvar, label = "*"), 
            vjust = 0.8, hjust = 0.5, size = 9) +
  scale_fill_gradientn(colors = custom_colors, limits = c(lwth, upth)) +  # 使用自定义颜色
  scale_y_discrete(limits = y_levels,  
                   labels = c("SDQ_sum_z" = "SDQ Total Score","SDQ_ES_sum_z" = "Emotional Symptoms", 
                              "SDQ_PP_sum_z" = "Peer Problems","SDQ_CP_sum_z" = "Conduct Problems",
                              "SDQ_H_sum_z" = "Hyperactivity","SDQ_PB_sum_z" = "Prosocial Behavior")) +
  labs(title = "Correlation of Executive Functions and Mental Health across Age", 
       x = "Tasks", y = "Mental Health") +
  theme(axis.line = element_blank(),
        aspect.ratio = 1.6,
        axis.text.x = element_text(size = 14, hjust = 0.5),
        axis.text.y = element_text(size = 14, hjust = 0.5, vjust = 0.5),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 22, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linewidth = 0),
        panel.grid.minor = element_line(linewidth = 1))

# 打印并保存图像
print(Fig)
ggsave(paste0(FigureFolder, "/corr_int.pdf"), plot = Fig, width = 18, height =24 , units = "cm")



# 确定颜色的上下限
lwth <- min(int.results.df$Interaction, na.rm = TRUE)
upth <- max(int.results.df$Interaction, na.rm = TRUE)

# 自定义量表顺序
y_levels <- c("APSS_sum_z", "SDQ_PB_sum_z", "SDQ_H_sum_z", "SDQ_CP_sum_z", "SDQ_PP_sum_z", "SDQ_ES_sum_z", "SDQ_sum")  # 量表顺序

# 创建整合的任务数据框
combined_data <- int.results.df
combined_data$sig <- as.logical(combined_data$sig)  # 转换 sig 列为逻辑值
combined_data$parcel <- factor(combined_data$parcel, levels = c("Twoback_acc_deviationZ", "Oneback_acc_deviationZ",  "d_prime_deviationZ"), 
                               labels = c("2-back", "1-back", "Go/No-go"))  # 设置任务的顺序和标签

# # 绘制3/4环形热图
# figcircle <- ggplot(combined_data, aes(x = parcel, y = dependentvar, fill = Interaction)) +  
#   geom_tile() +  
#   scale_fill_gradientn(colors = c("#f5dfdb", "#edb8b0", "#e69191", "#c25759"),
#                        name = "IntpartialRsq")  +  
#   coord_polar() +  
#   theme_bw()+  
#   theme(panel.border = element_blank(),  # 隐藏图形的边框        
#         panel.grid = element_blank(),   # 隐藏背景网格线        
#         axis.title = element_blank(),
#         legend.position = c(-0.1, 0.1),  # 图例位置，左下角
#         legend.title = element_text(size = 10),  # 图例标题字体大小
#         legend.text = element_text(size = 10)) + # 隐藏轴标题  
#   geom_text(data = combined_data[combined_data$sig == TRUE, ], 
#             aes(x = parcel, y = dependentvar, label = "*"), 
#             vjust = 0.8, hjust = 0.5, size = 7) +
#   scale_y_discrete(expand=expansion(mult=c(0.3,0)),
#                    limits = y_levels,  
#                    labels = c("SDQ_sum" = "SDQ Total Score","SDQ_ES_sum" = "Emotional Symptoms", 
#                               "SDQ_PP_sum" = "Peer Problems","SDQ_CP_sum" = "Conduct Problems",
#                               "SDQ_H_sum" = "Hyperactivity","SDQ_PB_sum" = "Prosocial Behavior"))+  
#   scale_x_discrete(expand=expansion(mult=c(1,0)))
# print(figcircle)
# ggsave(paste0(FigureFolder, "/corr_int_circle.tiff"), plot = figcircle, width = 16, height =24 , units = "cm")

####post plot
# 遍历所有显著交互效应
for (k in 1:nrow(significant_interactions)) {
  # 获取当前显著结果的数据
  current_slope_data <- significant_interactions$slope_data[[k]]
  current_task <- significant_interactions$dataname[k]
  current_variable <- significant_interactions$dependentvar[k]
  
  # 计算中位数和置信区间
  slope_summary <- current_slope_data %>%
    group_by(Age_year) %>%
    summarise(
      median_slope = median(slope, na.rm = TRUE),
      lower_95CI = quantile(slope, probs = 0.025, na.rm = TRUE),
      upper_95CI = quantile(slope, probs = 0.975, na.rm = TRUE)
    )
  
  # 绘制斜率随年龄变化的趋势
  slope_plot <- ggplot(slope_summary, aes(x = Age_year)) +
    geom_line(aes(y = median_slope), color = "black", size = 1.5) + 
    geom_ribbon(aes(ymin = lower_95CI, ymax = upper_95CI), alpha = 0.3, fill = "grey80") + 
    labs(title = paste0(current_task, " (", current_variable, ")"), 
         x = "", 
         y = "Slope") +
    scale_x_continuous(name="", limits = c(11, 18),breaks = seq(11, 18, by = 1))+
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.5, color = "black"),  
      axis.text = element_text(size = 20,color = "black",family = "Arial"), 
      axis.title = element_text(size = 20,family = "Arial"), 
      plot.title = element_text(size = 22, hjust = 0.5,family = "Arial"),
      panel.grid.major = element_blank(),    # 删除主网格线
      panel.grid.minor = element_blank(),    # 删除次网格线
      panel.border = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.5),  # 设置小线段颜色和大小
      axis.ticks.length = unit(0.2, "cm") 
    )
  
  print(slope_plot)
  # 保存图像
  ggsave(paste0(FigureFolder, "/significant_slope_change_", current_task, "_", current_variable, ".pdf"), plot = slope_plot, width = 16, height = 12, units = "cm")
}


significant_interactions <- read_rds(paste0(resultFolder, "/int_significant_YF.rds"))
nosignificant_interactions <- read_rds(paste0(resultFolder, "/int_Nosignificant_YF.rds"))
for (k in 1:nrow(significant_interactions)) {
  # 获取当前显著结果的数据
  current_slope_data <- significant_interactions$slope_data[[k]]
  current_task <- significant_interactions$dataname[k]
  current_variable <- significant_interactions$dependentvar[k]
  
  modified_task_name <- ifelse(current_task %in% names(task_mapping), task_mapping[current_task], current_task)
  modified_variable_name <- ifelse(current_variable %in% names(variable_mapping), variable_mapping[current_variable], current_variable)
  
  # 计算中位数和置信区间
  slope_summary <- current_slope_data %>%
    group_by(Age_year) %>%
    summarise(
      median_slope = median(slope, na.rm = TRUE),
      lower_95CI = quantile(slope, probs = 0.025, na.rm = TRUE),
      upper_95CI = quantile(slope, probs = 0.975, na.rm = TRUE)
    )
  
  # 绘制斜率随年龄变化的趋势
  slope_plot <- ggplot(slope_summary, aes(x = Age_year)) +
    geom_line(aes(y = median_slope), color = "black", size = 0.5) + 
    geom_ribbon(aes(ymin = lower_95CI, ymax = upper_95CI), alpha = 0.3, fill = "grey80") + 
    labs(title = paste0(modified_variable_name, " ~ ", modified_task_name), 
         x = "", 
         y = "Slope") +
    scale_x_continuous(name="", limits = c(11, 18),breaks = seq(11, 18, by = 1))+
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.25, color = "black"),  
      axis.text = element_text(size = 9,color = "black"), 
      axis.title = element_text(size = 9), 
      plot.title = element_text(size = 9, hjust = 0.5),
      panel.grid.major = element_blank(),    # 删除主网格线
      panel.grid.minor = element_blank(),    # 删除次网格线
      panel.border = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.25),  # 设置小线段颜色和大小
      axis.ticks.length = unit(0.05, "cm") 
    )
  
  print(slope_plot)
  # 保存图像
  ggsave(paste0(FigureFolder, "/significant_slope_change_", current_task, "_", current_variable, ".pdf"), plot = slope_plot, width = 6, height = 5.2, units = "cm")
}


# 遍历所有显著交互效应
task_mapping <- c(
  "GNGd" = "Go/No-Go",
  "back1" = "1-back",
  "back2" = "2-back"
)

variable_mapping <- c(
  "SDQ_ES_sum_z" = "Emotional Symptoms",
  "SDQ_PP_sum_z" = "Peer Problems",
  "SDQ_CP_sum_z" = "Conduct Problems",
  "SDQ_H_sum_z" = "Hyperactivity",
  "SDQ_PB_sum_z" = "Prosocial Behavior"
)

# 合并显著和不显著交互效应的斜率数据
all_slope_data <- c(significant_interactions$slope_data, nosignificant_interactions$slope_data)

# 获取所有 Y 轴刻度标签的最大宽度
all_y_labels <- unlist(lapply(all_slope_data, function(df) {
  y_values <- df$slope
  y_breaks <- pretty(y_values)  # 生成 Y 轴刻度
  max(nchar(sprintf("%.2f", y_breaks)))  # 格式化标签并找到最长的标签宽度
}))
max_y_label_width <- max(all_y_labels, na.rm = TRUE)  # 全局最大标签宽度

# 定义统一的 Y 轴标签格式化函数
format_y_labels <- function(y) {
  sprintf(paste0("%", max_y_label_width, ".2f"), y)  # 格式化标签为固定宽度
}

# # 定义统一的主题
# custom_theme <- theme_minimal() +
#   theme(
#     axis.line = element_line(size = 0.25, color = "black"),  
#     axis.text.x = element_text(size = 9, color = "black"), 
#     axis.text.y = element_text(size = 9, color = "black"),  # 使用等宽字体
#     axis.title = element_text(size = 9), 
#     plot.title = element_text(size = 9, hjust = 0.5),
#     panel.grid.major = element_blank(),    
#     panel.grid.minor = element_blank(),    
#     panel.border = element_blank(),
#     axis.ticks = element_line(color = "black", size = 0.25),  
#     axis.ticks.length = unit(0.05, "cm") 
#   )
for (k in 1:nrow(significant_interactions)) {
  # 获取当前显著结果的数据
  current_slope_data <- significant_interactions$slope_data[[k]]
  current_task <- significant_interactions$dataname[k]
  current_variable <- significant_interactions$dependentvar[k]
  
  modified_task_name <- ifelse(current_task %in% names(task_mapping), task_mapping[current_task], current_task)
  modified_variable_name <- ifelse(current_variable %in% names(variable_mapping), variable_mapping[current_variable], current_variable)
  
  # 计算中位数和置信区间
  slope_summary <- current_slope_data %>%
    group_by(Age_year) %>%
    summarise(
      median_slope = median(slope, na.rm = TRUE),
      lower_95CI = quantile(slope, probs = 0.025, na.rm = TRUE),
      upper_95CI = quantile(slope, probs = 0.975, na.rm = TRUE)
    )
  
  # 绘制斜率随年龄变化的趋势
  slope_plot <- ggplot(slope_summary, aes(x = Age_year)) +
    geom_line(aes(y = median_slope), color = "black", size = 0.5) + 
    geom_ribbon(aes(ymin = lower_95CI, ymax = upper_95CI), alpha = 0.3, fill = "grey80") + 
    labs(title = paste0(modified_variable_name, " ~ ", modified_task_name), 
         x = NULL, 
         y = "Slope") +
    scale_x_continuous(name="", limits = c(11, 18),breaks = seq(11, 18, by = 2))+
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.25, color = "black"),  
      axis.text = element_text(size = 9,color = "black"), 
      axis.title.y = element_text(size = 9), 
      axis.title.x = element_blank(), 
      plot.title = element_text(size = 9, hjust = 0.5),
      panel.grid.major = element_blank(),    # 删除主网格线
      panel.grid.minor = element_blank(),    # 删除次网格线
      panel.border = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.25),  # 设置小线段颜色和大小
      axis.ticks.length = unit(0.05, "cm") 
    )
  
  print(slope_plot)
  # 保存图像
  ggsave(paste0(FigureFolder, "/significant_slope_change_", current_task, "_", current_variable, ".pdf"), plot = slope_plot, width = 6, height = 5.2, units = "cm")
}
# 遍历所有不显著交互效应
for (k in 1:nrow(nosignificant_interactions)) {
  # 获取当前显著结果的数据
  current_slope_data <- nosignificant_interactions$slope_data[[k]]
  current_task <- nosignificant_interactions$dataname[k]
  current_variable <- nosignificant_interactions$dependentvar[k]
  
  modified_task_name <- ifelse(current_task %in% names(task_mapping), task_mapping[current_task], current_task)
  modified_variable_name <- ifelse(current_variable %in% names(variable_mapping), variable_mapping[current_variable], current_variable)
  
  # 计算中位数和置信区间
  slope_summary <- current_slope_data %>%
    group_by(Age_year) %>%
    summarise(
      median_slope = median(slope, na.rm = TRUE),
      lower_95CI = quantile(slope, probs = 0.025, na.rm = TRUE),
      upper_95CI = quantile(slope, probs = 0.975, na.rm = TRUE)
    )
  
  # 绘制斜率随年龄变化的趋势
  slope_plot <- ggplot(slope_summary, aes(x = Age_year)) +
    geom_line(aes(y = median_slope), color = "black", size = 0.5) + 
    geom_ribbon(aes(ymin = lower_95CI, ymax = upper_95CI), alpha = 0.3, fill = "grey80") + 
    labs(title = paste0(modified_variable_name, " ~ ", modified_task_name), 
         x = NULL, 
         y = "Slope") +
    scale_x_continuous(name="", limits = c(11, 18),breaks = seq(11, 18, by = 2))+
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.25, color = "black"),  
      axis.text = element_text(size = 9,color = "black"), 
      axis.title.y = element_text(size = 9), 
      axis.title.x = element_blank(), 
      plot.title = element_text(size = 9, hjust = 0.5),
      panel.grid.major = element_blank(),    # 删除主网格线
      panel.grid.minor = element_blank(),    # 删除次网格线
      panel.border = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.25),  # 设置小线段颜色和大小
      axis.ticks.length = unit(0.05, "cm") 
    )
  
  print(slope_plot)
  # 保存图像
  ggsave(paste0(FigureFolder, "/Nosignificant_slope_change_", current_task, "_", current_variable, ".pdf"), plot = slope_plot, width = 6, height = 5.2, units = "cm")
}
