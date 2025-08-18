# --- 0. 加载所有必要的 R 包 ---
library(ggplot2)
library(tidyverse)
library(mgcv)
library(readxl)
library(openxlsx)
library(cowplot) # Using cowplot for plot_grid, similar to the reference script

# --- 1. 文件路径设置 ---
# PLEASE UPDATE THESE PATHS TO MATCH YOUR FOLDER STRUCTURE
datapath <- 'D:/datasets/yunfu/raw_data'
interfileFolder <- "D:/datasets/yunfu/interfile_folder/temp"
functionFolder <- "D:/code/EF_Normative_Model/functions"
resultFolder <- "D:/datasets/yunfu/results/phq_sum_corr_modified"
FigureFolder <- 'D:/datasets/yunfu/figures/phq_sum_fig_modified'
# This is the new file you specified for the beta values
beta_file_path <- "D:/datasets/yunfu/results/psy_corr/PHQ_all_deviation_correlations_new.xlsx"

dir.create(resultFolder, showWarnings = FALSE)
dir.create(FigureFolder, showWarnings = FALSE)

# --- 2. 数据读取与预处理 ---

# 读取行为任务数据
GNGd_data <- readRDS(paste0(interfileFolder, '/GNGd_prime.deviations.rds'))
back1_data <- readRDS(paste0(interfileFolder, '/back1Acc.deviations.rds'))
back2_data <- readRDS(paste0(interfileFolder, '/back2Acc.deviations.rds'))

# 读取 PHQ 数据
PHQ_data <- read_xlsx(paste0(interfileFolder, '/PHQ.xlsx')) %>%
  select(用户ID, PHQ_sum) %>%
  mutate(用户ID = as.character(用户ID))

# 读取包含 Beta 值的相关性结果文件
beta_table <- read_xlsx(beta_file_path)

# 数据类型统一
GNGd_data <- GNGd_data %>% mutate(x__ID = as.character(x__ID), Sex = as.factor(Sex))
back1_data <- back1_data %>% mutate(x__ID = as.character(x__ID), Sex = as.factor(Sex))
back2_data <- back2_data %>% mutate(x__ID = as.character(x__ID), Sex = as.factor(Sex))

# 合并 PHQ 数据
GNGd_data <- GNGd_data %>% inner_join(PHQ_data, by = c("x__ID" = "用户ID"))
back1_data <- back1_data %>% inner_join(PHQ_data, by = c("x__ID" = "用户ID"))
back2_data <- back2_data %>% inner_join(PHQ_data, by = c("x__ID" = "用户ID"))

# 检查合并
if(nrow(GNGd_data) == 0) stop("ERROR: GNGd_data is empty after merging. Check ID columns.")

# 加载自定义函数
source(paste0(functionFolder, '/gam_varyingcoefficients.R'))

# --- 3. 定义变量与参数 ---
psyc_variables_continous <- c("PHQ_sum_z")

EFvars.set <- as.data.frame(matrix(c(
  "d_prime_deviationZ", "GNGd",
  "Oneback_acc_deviationZ", "back1",
  "Twoback_acc_deviationZ", "back2"
), byrow = TRUE, ncol = 2, dimnames = list(NULL, c("varname", "dataname"))))

# 对因变量进行 Z-score 标准化
GNGd_data[, "PHQ_sum_z"] <- scale(GNGd_data[, "PHQ_sum"])
back1_data[, "PHQ_sum_z"] <- scale(back1_data[, "PHQ_sum"])
back2_data[, "PHQ_sum_z"] <- scale(back2_data[, "PHQ_sum"])

psyc_variables_continous_z <- "PHQ_sum_z"

# 设置模型参数
knots <- 3
set_fx <- FALSE
increments <- 1000
draws <- 1000
return_posterior_coefficients <- TRUE

# --- 4. 初始化结果存储 ---
# MODIFICATION: Store all results in a list, like in the reference script
interaction_results <- list()

# --- 5. 主循环：执行 GAM 分析 ---
for (i in 1:nrow(EFvars.set)) {
  int_var <- EFvars.set$varname[i]
  dataname_str <- EFvars.set$dataname[i]
  dataname <- paste0(dataname_str, "_data")
  
  for (j in 1:length(psyc_variables_continous_z)) {
    dependentvar <- psyc_variables_continous_z[j]
    smooth_var <- "Age_year"
    covariates <- "Sex"
    
    cat(paste("Processing:", dataname, "~", int_var, "by", smooth_var, "for dependent var:", dependentvar, "\n"))
    
    # 调用 GAM 变系数模型函数
    result <- gam.varyingcoefficients(
      dependentvar = dependentvar,
      dataname = dataname,
      smooth_var = smooth_var,
      int_var = int_var,
      covariates = covariates,
      knots = knots,
      set_fx = set_fx,
      increments = increments,
      draws = draws,
      return_posterior_coefficients = return_posterior_coefficients
    )
    
    # 保存结果到列表
    model_name <- paste(dataname_str, dependentvar, sep = "_")
    interaction_results[[model_name]] <- result
  }
}

# 保存所有结果为RDS文件，方便后续快速加载
saveRDS(interaction_results, file = paste0(resultFolder, "/all_interaction_results_PHQ.rds"))

# --- 6. 新的可视化流程 ---
# This entire section is replaced with the logic from EF_psy_int_post_slope_compare.R
# with corrections based on your feedback.

# (可选) 如果已经运行过上面的分析，可以直接加载结果
# interaction_results <- readRDS(paste0(resultFolder, "/all_interaction_results_PHQ.rds"))

# 定义绘图用的标签映射
task_mapping <- c("GNGd" = "Go/No-Go", "back1" = "1-back", "back2" = "2-back")
variable_mapping <- c("PHQ_sum_z" = "PHQ Total Score")

# 定义统一的绘图主题
custom_theme <- theme_minimal() +
  theme(
    axis.line = element_line(size = 0.25, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.25),
    axis.ticks.length = unit(0.05, "cm")
  )

# 遍历所有模型结果并生成图表
for (model_name in names(interaction_results)) {
  result <- interaction_results[[model_name]]
  
  # 跳过没有有效结果的模型
  if (is.null(result) || length(result) < 2) {
    message(paste("Skipping empty result for model:", model_name))
    next
  }
  
  # 从模型名称中提取任务和变量信息
  parts <- strsplit(model_name, "_")[[1]]
  current_task <- parts[1]
  current_variable_raw <- paste(parts[2:length(parts)], collapse = "_")
  
  # 获取绘图所需的数据
  current_slope_data <- result[[2]]
  
  # --- *** KEY CORRECTION STARTS HERE *** ---
  
  # 1. 移除 "_z" 后缀以匹配 beta table 中的 'parcel' 列
  lookup_parcel_name <- current_variable_raw
  # lookup_parcel_name <- sub("_z$", "", current_variable_raw)
  
  # 2. 从beta_table中获取对应的总体beta值
  #    - 使用 'dataset' 列来匹配任务 (current_task)
  #    - 使用 'parcel' 列来匹配处理后的因变量名 (lookup_parcel_name)
  beta_value <- beta_table %>%
    filter(dataset == current_task, parcel == lookup_parcel_name) %>%
    pull(beta)
  
  # --- *** KEY CORRECTION ENDS HERE *** ---
  
  # 如果找不到beta值，跳过当前循环并打印更详细的错误信息
  if (length(beta_value) == 0) {
    message(paste("Could not find a beta value for dataset =", current_task, "and parcel =", lookup_parcel_name, ". Skipping plot."))
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
      is_significant = (beta_value < lower_95CI | beta_value > upper_95CI),
      abs_diff = abs(median_slope - beta_value)
    )
  
  # 绘制主图：斜率随年龄的变化趋势
  slope_plot_main <- ggplot(slope_summary, aes(x = Age_year)) +
    geom_line(aes(y = median_slope), color = "black", size = 0.5) +
    geom_ribbon(aes(ymin = lower_95CI, ymax = upper_95CI), alpha = 0.3, fill = "grey80") +
    geom_hline(yintercept = beta_value, linetype = "dashed", color = "red", size = 0.4) +
    labs(
      title = paste0(variable_mapping[current_variable_raw], " ~ ", task_mapping[current_task]),
      x = "",
      y = "Slope"
    ) +
    scale_x_continuous(name = "", limits = c(11, 18), breaks = seq(11, 18, by = 2)) +
    custom_theme
  
  # 绘制显著性条形图（子图）
  slope_plot_bar <- ggplot(slope_summary, aes(x = Age_year, y = 1)) +
    geom_col(aes(fill = abs_diff), data = filter(slope_summary, is_significant == TRUE),
             width = 1, color = "transparent") +
    geom_col(data = filter(slope_summary, is_significant == FALSE), fill = "white",
             width = 1, color = "transparent") +
    scale_fill_gradient(low = "#d9e0f0", high = "#A4C5DF", na.value = "white", guide = "none") +
    scale_x_continuous(name = "Age (years)", limits = c(11, 18), breaks = seq(11, 18, by = 2)) +
    theme_void() +
    theme(
      # --- 修改开始 ---
      axis.title.x = element_text(size = 14), # 恢复显示X轴标题 "Age (years)"
      axis.text.x = element_blank(),         # 保持移除X轴的刻度数字
      # --- 修改结束 ---
      plot.margin = margin(t = 0, r = 15, b = 5, l = 30, unit = "pt"),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
    )
  
  # 使用 cowplot 合并主图和子图
  final_plot <- plot_grid(slope_plot_main, slope_plot_bar, ncol = 1, rel_heights = c(0.90, 0.10), align = "v")
  
  # 打印并保存
  print(final_plot)
  ggsave(
    filename = paste0(FigureFolder, "/slope_comparison_", current_task, "_", lookup_parcel_name, ".pdf"),
    plot = final_plot,
    width = 6, height = 6.5, units = "cm"
  )
}

cat("Script finished. All plots have been generated and saved to:", FigureFolder, "\n")


# hist figure



# 假设 beta_table 已经从 PHQ_all_deviation_correlations_new.xlsx 文件中加载

# --- 1. 准备用于绘图的新数据框 ---
# 我们从 beta_table 中筛选出关于 PHQ_sum_z 的结果
plot_data_beta <- beta_table %>%
  filter(parcel == "PHQ_sum_z") %>% # 筛选出我们关心的连续变量结果
  mutate(
    # 创建显著性标签，使用置换检验的p值 (p_value_anova)
    sig_label = ifelse(p_value_anova < 0.05, "*", ""),
    # 将任务变量(dataset)转换为因子，以确保绘图顺序和标签正确
    task_factor = factor(dataset, 
                         levels = c("GNGd", "back1", "back2"),
                         labels = c("Go/No-Go", "1-back", "2-back"))
  )

# --- 2. 创建图表 ---
# 沿用您提供的代码框架，修改Y轴变量和数据源
plot_beta_values <- ggplot(plot_data_beta, aes(x = task_factor, y = beta, fill = task_factor)) +
  geom_col() + # 使用 geom_col() 直接绘制数值
  # 在条形图上方添加显著性标记 '*'
  geom_text(aes(label = sig_label), vjust = -0.5, size = 8, color = "black") +
  # 手动设置填充颜色
  scale_fill_manual(values = c("Go/No-Go" = "#dbe6eb", "1-back" = "#a1c1d5", "2-back" = "#3a7ca5")) +
  # 设置Y轴标签为 β (beta)，并让ggplot自动调整范围
  # 注意：移除了 limits 参数，让ggplot根据beta值的范围自动设定
  scale_y_continuous(expression(beta), limits = c(0, 0.035), expand = c(0, 0)) +
  # 设置X轴标签
  labs(x = "EF task", y = expression(beta)) +
  # 使用经典主题，并进行精细调整
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none", # 隐藏图例
    axis.line = element_line(size = 0.5),
    axis.ticks = element_line(size = 0.5),
    # 确保Y轴从0开始
    axis.line.y = element_line(arrow = arrow(length = unit(0.3, "cm"))),
    axis.ticks.y = element_line(size=0.5)
  )

# --- 3. 打印图表 ---
print(plot_beta_values)


