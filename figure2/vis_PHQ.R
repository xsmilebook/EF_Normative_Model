rm(list=ls())
library(ggplot2)
library(tidyverse)
library(mgcv)
library(readxl)
library(parallel)
library(gamlss)
library(scales)
library(tableone)
library(openxlsx)
# input directory
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- '/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/EF_Yunfu/interfileFolder_back12before'
  datapath_GNG <- '/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_all/interfileFolder'
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/EF_Yunfu/interfileFolder_back12before"
  functionFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/EF_Yunfu/code/Rcode/functions"
  resultFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/EF_Yunfu/results/corr"
}else{
  datapath <- 'D:/datasets/yunfu/raw_data'
  FigureFolder <- 'D:/datasets/yunfu/figures/fig1'
  interfileFolder <- "D:/datasets/yunfu/interfile_folder/temp"
  functionFolder <- "D:/code/EF_Normative_Model/functions"
  resultFolder <- "D:/datasets/yunfu/results/psy_corr"
}

# source functions
#source(paste0(functionFolder, "/gam_varyingcoefficients.R"))
source(paste0(functionFolder, "/gamcog_withsmoothvar_deviation.R"))
#source(paste0(functionFolder, "/ordinalcorr_new.R"))
# import dataset
#switch_data <- read_xlsx(paste0(datapath, '/Q_switch.xlsx'))
#head(switch_data)
PHQ_data <- read_xlsx(paste0(interfileFolder, '/PHQ.xlsx'))
GNGd_data <- read_rds(paste0(interfileFolder, '/GNGd_prime.deviations.rds'))
back1_data <- read_rds(paste0(interfileFolder, '/back1Acc.deviations.rds'))
back2_data <- read_rds(paste0(interfileFolder, '/back2Acc.deviations.rds'))

PHQ_data <- PHQ_data %>%
  select(用户ID, PHQ_y09) %>%  
  mutate(用户ID = as.character(用户ID))

GNGd_data <- GNGd_data %>% mutate(x__ID = as.character(x__ID))
back1_data <- back1_data %>% mutate(x__ID = as.character(x__ID))
back2_data <- back2_data %>% mutate(x__ID = as.character(x__ID))

GNGd_data <- inner_join(
  GNGd_data,
  PHQ_data,
  by = c("x__ID" = "用户ID")
)

back1_data <- inner_join(
  back1_data,
  PHQ_data,
  by = c("x__ID" = "用户ID")
)

back2_data <- inner_join(
  back2_data,
  PHQ_data,
  by = c("x__ID" = "用户ID")
)

plot_data <- bind_rows(
  GNGd_data %>% 
    select(x__ID, PHQ_y09, d_prime_deviationZ) %>%
    rename(acc_deviationZ = d_prime_deviationZ) %>%
    mutate(task = "GNGd"),
  
  back1_data %>% 
    select(x__ID, PHQ_y09, Oneback_acc_deviationZ) %>%
    rename(acc_deviationZ = Oneback_acc_deviationZ) %>%
    mutate(task = "Back1"),
  
  back2_data %>% 
    select(x__ID, PHQ_y09, Twoback_acc_deviationZ) %>%
    rename(acc_deviationZ = Twoback_acc_deviationZ) %>%
    mutate(task = "Back2")
) %>%
  filter(!is.na(PHQ_y09))  # 过滤掉PHQ_y09为NA的行

# 转换 PHQ_y09 为因子以用于分组
plot_data$PHQ_y09 <- as.factor(plot_data$PHQ_y09)

# 定义绘图函数 - 修正版本
plot_task <- function(data, task_name) {
  # 创建数据副本用于散点图，将PHQ_y09稍微向右偏移
  data_scatter <- data %>% 
    filter(task == task_name) %>%
    mutate(PHQ_y09_numeric = as.numeric(PHQ_y09) + 0.2)  # 向右偏移0.2
  
  # 箱线图层
  p_box <- ggplot(data %>% filter(task == task_name), 
                  aes(x = PHQ_y09, y = acc_deviationZ)) +
    geom_boxplot(color = "black", outlier.shape = NA, alpha = 0.6) +
    theme_minimal() +
    labs(title = task_name, x = "PHQ_y09", y = "Deviation Z-score") +
    ylim(-3, 3) +
    theme(axis.title.x = element_blank())
  
  # 散点图层
  p_scatter <- ggplot(data_scatter, 
                      aes(x = PHQ_y09_numeric, y = acc_deviationZ)) +
    geom_point(color = "lightblue", alpha = 0.6, size = 1) +
    theme_minimal() +
    labs(title = task_name, x = "PHQ_y09", y = "Deviation Z-score") +
    ylim(-3, 3) +
    theme(axis.title.x = element_blank())
  
  # 组合两个图层
  p_combined <- ggplot(data %>% filter(task == task_name), 
                       aes(x = PHQ_y09, y = acc_deviationZ)) +
    geom_boxplot(color = "black", outlier.shape = NA, alpha = 0.6) +
    geom_jitter(data = data %>% filter(task == task_name), 
                aes(x = as.numeric(PHQ_y09) + 0.2, y = acc_deviationZ),
                color = "lightblue", alpha = 0.6, size = 1, width = 0.1) +
    theme_minimal() +
    labs(title = task_name, x = "PHQ_y09", y = "Deviation Z-score") +
    ylim(-3, 3)
  
  return(p_combined)
}

plot_task <- function(data, task_name, jitter_width = 0.15, scatter_offset = 0.5) {
  ggplot(data %>% filter(task == task_name), 
         aes(x = PHQ_y09, y = acc_deviationZ)) +
    geom_boxplot(color = "black", 
                 outlier.shape = NA, 
                 alpha = 0.6,
                 width = 0.6) +
    geom_jitter(aes(x = as.numeric(PHQ_y09) + scatter_offset),
                color = "lightblue", 
                alpha = 0.6, 
                size = 1, 
                width = jitter_width,
                height = 0) +
    theme_minimal() +
    labs(title = task_name, x = "PHQ_y09", y = "Deviation Z-score") +
    ylim(-3, 3) +

    scale_x_discrete(expand = expansion(add = c(0.8, 0.8)))  # 增加左右边距
}

# 绘制三张图
p1 <- plot_task(plot_data, "GNGd")
p2 <- plot_task(plot_data, "Back1")
p3 <- plot_task(plot_data, "Back2")

# 组合显示三张图
combined_plot <- p1 + p2 + p3 +
  plot_layout(ncol = 3) +
  plot_annotation(title = "Deviation Z-scores across PHQ scores for different tasks")

# 显示图形
print(combined_plot)


simple_pairwise_test <- function(data, task_name) {
  cat("\n=== ", task_name, "任务的两两比较 ===\n")
  
  # 过滤数据
  task_data <- data %>% filter(task == task_name)
  
  # 获取PHQ_y09的所有水平
  phq_levels <- sort(unique(task_data$PHQ_y09))
  
  # 进行所有可能的两两比较
  results <- data.frame()
  
  for (i in 1:(length(phq_levels)-1)) {
    for (j in (i+1):length(phq_levels)) {
      group1 <- phq_levels[i]
      group2 <- phq_levels[j]
      
      # 提取两组数据
      data1 <- task_data %>% filter(PHQ_y09 == group1) %>% pull(acc_deviationZ)
      data2 <- task_data %>% filter(PHQ_y09 == group2) %>% pull(acc_deviationZ)
      
      # 进行t检验
      t_test_result <- t.test(data1, data2)
      
      # 记录结果
      temp_result <- data.frame(
        Task = task_name,
        Group1 = group1,
        Group2 = group2,
        Mean1 = round(mean(data1, na.rm = TRUE), 3),
        Mean2 = round(mean(data2, na.rm = TRUE), 3),
        P_value = round(t_test_result$p.value, 4),
        Significant = ifelse(t_test_result$p.value < 0.05, "Yes", "No"),
        stringsAsFactors = FALSE
      )
      
      results <- rbind(results, temp_result)
    }
  }
  
  # 显示结果
  print(results)
  
  # 显示显著的结果
  significant_results <- results %>% filter(Significant == "Yes")
  if (nrow(significant_results) > 0) {
    cat("\n显著差异 (p < 0.05):\n")
    print(significant_results)
  } else {
    cat("\n没有发现显著差异 (p < 0.05)\n")
  }
  
  return(results)
}

# 对每个任务进行分析
cat("开始简化统计分析...\n")

results_GNGd <- simple_pairwise_test(plot_data, "GNGd")
results_Back1 <- simple_pairwise_test(plot_data, "Back1")
results_Back2 <- simple_pairwise_test(plot_data, "Back2")

# 合并所有结果并保存
all_results <- rbind(results_GNGd, results_Back1, results_Back2)

# 保存结果
output_file <- paste0(resultFolder, "/pairwise_comparison_results.xlsx")
write.xlsx(all_results, file = output_file, rowNames = FALSE)
cat("\n结果已保存到:", output_file, "\n")

# 打印汇总结果
cat("\n=== 汇总结果 ===\n")
for (task in c("GNGd", "Back1", "Back2")) {
  task_results <- all_results %>% filter(Task == task & Significant == "Yes")
  cat("\n", task, "任务中的显著差异:\n")
  if (nrow(task_results) > 0) {
    print(task_results[, c("Group1", "Group2", "P_value")])
  } else {
    cat("  无显著差异\n")
  }
}














# 安装和加载必要的包
if (!require(ggsignif)) install.packages("ggsignif")
library(ggsignif)

# 简化的统计分析函数，返回显著差异结果
get_significant_pairs <- function(data, task_name) {
  task_data <- data %>% filter(task == task_name)
  phq_levels <- sort(unique(task_data$PHQ_y09))
  significant_pairs <- list()
  
  for (i in 1:(length(phq_levels)-1)) {
    for (j in (i+1):length(phq_levels)) {
      group1 <- phq_levels[i]
      group2 <- phq_levels[j]
      
      # 提取两组数据
      data1 <- task_data %>% filter(PHQ_y09 == group1) %>% pull(acc_deviationZ)
      data2 <- task_data %>% filter(PHQ_y09 == group2) %>% pull(acc_deviationZ)
      
      # 进行t检验
      t_test_result <- t.test(data1, data2)
      
      # 如果显著，记录
      if (t_test_result$p.value < 0.05) {
        significant_pairs[[length(significant_pairs) + 1]] <- c(as.character(group1), as.character(group2))
      }
    }
  }
  
  return(significant_pairs)
}

# 修正的绘图函数，包含显著性标注
plot_task_with_significance <- function(data, task_name, jitter_width = 0.15, scatter_offset = 0.25) {
  # 获取显著差异对
  sig_pairs <- get_significant_pairs(data, task_name)
  
  # 基础图形
  p <- ggplot(data %>% filter(task == task_name), 
              aes(x = PHQ_y09, y = acc_deviationZ)) +
    geom_boxplot(color = "black", 
                 outlier.shape = NA, 
                 alpha = 0.6,
                 width = 0.6) +
    geom_jitter(aes(x = as.numeric(PHQ_y09) + scatter_offset),
                color = "lightblue", 
                alpha = 0.6, 
                size = 1, 
                width = jitter_width,
                height = 0) +
    theme_minimal() +
    labs(title = task_name, x = "PHQ_y09", y = "Deviation Z-score") +
    ylim(-3, 3)
  
  # 如果有显著差异，添加标注
  if (length(sig_pairs) > 0) {
    # 转换显著对为ggsignif需要的格式
    comparisons <- lapply(sig_pairs, function(pair) {
      as.numeric(pair)
    })
    
    p <- p + ggsignif::geom_signif(
      comparisons = comparisons,
      map_signif_level = TRUE,
      textsize = 3,
      tip_length = 0.02,
      y_position = seq(max(data$acc_deviationZ, na.rm = TRUE) + 0.2, 
                       max(data$acc_deviationZ, na.rm = TRUE) + 0.2 * length(comparisons), 
                       0.2),
      vjust = 0.5
    )
  }
  
  return(p)
}

# 处理y轴位置避免重叠的函数
plot_task_with_significance <- function(data, task_name, jitter_width = 0.15, scatter_offset = 0.25) {
  # 获取显著差异对
  sig_pairs <- get_significant_pairs(data, task_name)
  
  # 基础图形
  p <- ggplot(data %>% filter(task == task_name), 
              aes(x = PHQ_y09, y = acc_deviationZ)) +
    geom_boxplot(color = "black", 
                 outlier.shape = NA, 
                 alpha = 0.6,
                 width = 0.6) +
    geom_jitter(aes(x = as.numeric(PHQ_y09) + scatter_offset),
                color = "lightblue", 
                alpha = 0.6, 
                size = 1, 
                width = jitter_width,
                height = 0) +
    theme_minimal() +
    labs(title = task_name, x = "PHQ_y09", y = "Deviation Z-score") +
    ylim(-3, 3)
  
  # 如果有显著差异，添加标注
  if (length(sig_pairs) > 0) {
    # 计算y轴位置，避免重叠
    max_y <- max(data$acc_deviationZ, na.rm = TRUE)
    y_positions <- seq(max_y + 0.3, max_y + 0.3 + 0.2 * (length(sig_pairs) - 1), 0.2)
    
    # 转换显著对为字符向量
    comparisons_list <- lapply(sig_pairs, as.character)
    
    p <- p + ggsignif::geom_signif(
      comparisons = comparisons_list,
      map_signif_level = TRUE,
      textsize = 3,
      tip_length = 0.02,
      y_position = y_positions,
      vjust = 0.5
    )
  }
  
  return(p)
}

# 绘制三张图
p1 <- plot_task_with_significance(plot_data, "GNGd")
p2 <- plot_task_with_significance(plot_data, "Back1")
p3 <- plot_task_with_significance(plot_data, "Back2")

# 组合显示三张图
combined_plot <- p1 + p2 + p3 +
  plot_layout(ncol = 3) +
  plot_annotation(title = "Deviation Z-scores across PHQ scores for different tasks")

# 显示图形
print(combined_plot)

# 保存图形
ggsave(filename = paste0(resultFolder, "/PHQ_deviation_plots_with_significance.png"), 
       plot = combined_plot, 
       width = 15, height = 5, dpi = 300)

# 同时输出显著性检验结果到控制台
cat("显著性检验结果:\n")
for (task in c("GNGd", "Back1", "Back2")) {
  cat("\n=== ", task, " ===\n")
  sig_pairs <- get_significant_pairs(plot_data, task)
  if (length(sig_pairs) > 0) {
    for (i in 1:length(sig_pairs)) {
      cat("PHQ", sig_pairs[[i]][1], " vs PHQ", sig_pairs[[i]][2], "\n")
    }
  } else {
    cat("无显著差异\n")
  }
}



