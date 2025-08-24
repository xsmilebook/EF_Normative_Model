# 1. 加载必要的R包
library(tidyverse)
library(ggplot2)
library(tools) # 用于处理文件名

# 2. 设置路径 (请确保这里的路径和您生成文件的路径一致)
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  # 服务器路径
  dataFolder <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/EF_Normative_Model/anova/results" # .rds文件所在文件夹
  FigureFolder <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/EF_Normative_Model/anova/figures" # 图片保存文件夹
} else {
  # 本地路径 (请根据您的实际情况修改)
  dataFolder <- "D:/datasets/yunfu/results/anovaPB"
  FigureFolder <- 'D:/datasets/yunfu/figures/fig2/anovaPB'
}

# 3. 检查并创建图片输出文件夹
if (!dir.exists(FigureFolder)) {
  dir.create(FigureFolder, recursive = TRUE)
}

# 4. 获取所有.rds文件的列表
files_to_process <- list.files(dataFolder, pattern = "\\.rds$", full.names = TRUE)

# 5. 循环处理每个文件并绘图
for (file_path in files_to_process) {
  # 读取RDS文件
  data <- readRDS(file_path)
  
  # 提取所需变量
  sim_df <- data.frame(simulated_stats = data$simulation$simulated_stats)
  observed_value <- data$simulation$observed_stat
  n_sim <- data$simulation$n_sim
  
  # 格式化p值
  p_value <- mean(sim_df$simulated_stats >= observed_value)
  if (p_value == 0) {
    p_label <- paste0("p < ", formatC(1/n_sim, format = "e", digits = 0), " ")
  } else {
    # 使用nsmall=3确保p值总是有3位小数，便于对齐
    p_label <- paste0(" p = ", format(p_value, digits = 3, nsmall = 3), " ")
  }
  
  # 自动判断文本注释的左右位置
  median_sim <- median(sim_df$simulated_stats)
  if (observed_value > median_sim) {
    hjust_value <- 1.1 # 观测值在右，文本在左
  } else {
    hjust_value <- -0.1 # 观测值在左，文本在右
  }
  
  # 预计算直方图高度
  hist_data <- hist(sim_df$simulated_stats, breaks = 50, plot = FALSE)
  max_y_height <- max(hist_data$counts)
  
  # 使用ggplot2绘图
  p <- ggplot(sim_df, aes(x = simulated_stats)) +
    geom_histogram(bins = 50, fill = "skyblue", alpha = 0.7, color = "white") +
    geom_vline(xintercept = observed_value, color = "black", linetype = "solid", linewidth = 1) +
    geom_point(aes(x = !!observed_value, y = 0), color = "black", size = 5, shape = 18) +
    
    labs(
      title = "Permutation Distribution vs. Observed Statistic",
      subtitle = paste("EFvar:", data$EFvar, " |  psyvar:", data$simulation$psyvar),
      x = "Statistic Value",
      y = "Frequency (Count)"
    ) +
    
    # 【修改点 2: 使用等宽字体对齐】
    annotate(
      "text",
      x = observed_value,
      y = max_y_height * 0.8,
      label = paste0("Observed = ", format(round(observed_value, 2), nsmall = 2), "\n", p_label),
      hjust = hjust_value,
      vjust = 0.5,
      color = "black",
      fontface = "bold",
      size = 4.5,
      # family = "mono" # 指定使用等宽字体
    ) +
    
    theme_minimal(base_size = 14) +
    # 【修改点 1: 美化背景网格线】
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      # 将主要网格线改为更淡的灰色虚线
      panel.grid.major = element_line(color = "grey90", linetype = "dashed", linewidth = 0.4),
      # 移除次要网格线
      panel.grid.minor = element_blank()
    )
  
  # 生成图片文件名
  output_filename <- paste0(file_path_sans_ext(basename(file_path)), ".png")
  output_path <- file.path(FigureFolder, output_filename)
  
  # 保存图片
  ggsave(output_path, plot = p, width = 8, height = 6, dpi = 300)
}

# 提示所有任务已完成
print(paste("所有", length(files_to_process), "个文件的可视化图片已生成并保存至:", FigureFolder))