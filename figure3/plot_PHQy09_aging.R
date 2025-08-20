# 清空环境
rm(list = ls())

# --- 0. 加载所有必要的 R 包 ---
library(ggplot2)
library(tidyverse)
library(mgcv)
library(readxl)
library(openxlsx)
library(gamm4) # <-- 需要这个包
library(patchwork) # <-- 用于合并图形

# --- 1. 文件路径设置 ---
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  # ... 
} else {
  interfileFolder <- "D:/datasets/yunfu/interfile_folder/temp"
  functionFolder <- "D:/code/EF_Normative_Model/functions"
  resultFolder <- "D:/datasets/yunfu/results/phq_y09_interaction"
  FigureFolder <- 'D:/datasets/yunfu/figures/phq_y09_fig'
}
dir.create(resultFolder, showWarnings = FALSE)
dir.create(FigureFolder, showWarnings = FALSE)

# --- 2. 加载自定义函数 ---
source(paste0(functionFolder, '/gam_ordinal_interaction.R')) 

# --- 3. 数据读取与准备 ---

# 读取行为任务数据
# 【重要】确保您的数据中有一个唯一的被试ID列，这里我们假设是 'x__ID'
GNGd_data <- readRDS(paste0(interfileFolder, '/GNGd_prime.deviations.rds'))
back1_data <- readRDS(paste0(interfileFolder, '/back1Acc.deviations.rds'))
back2_data <- readRDS(paste0(interfileFolder, '/back2Acc.deviations.rds'))

# 读取 PHQ 数据
PHQ_data <- read_xlsx(paste0(interfileFolder, '/PHQ.xlsx')) %>%
  select(用户ID, PHQ_y09) %>%  
  mutate(用户ID = as.character(用户ID))

# 数据类型统一和重命名ID
# 【重要】我们将 'x__ID' 重命名为 'ID' 以便函数内部使用
GNGd_data <- GNGd_data %>% mutate(ID = as.character(x__ID), Gender = as.factor(Sex))
back1_data <- back1_data %>% mutate(ID = as.character(x__ID), Gender = as.factor(Sex))
back2_data <- back2_data %>% mutate(ID = as.character(x__ID), Gender = as.factor(Sex))

# 合并 PHQ 数据
GNGd_data <- GNGd_data %>% inner_join(PHQ_data, by = c("ID" = "用户ID"))
back1_data <- back1_data %>% inner_join(PHQ_data, by = c("ID" = "用户ID"))
back2_data <- back2_data %>% inner_join(PHQ_data, by = c("ID" = "用户ID"))

# --- 4. 定义变量与参数 ---
EFvars.set <- as.data.frame(matrix(c(
  "d_prime_deviationZ", "GNGd",
  "Oneback_acc_deviationZ", "back1",
  "Twoback_acc_deviationZ", "back2"
), byrow = TRUE, ncol = 2, dimnames = list(NULL, c("varname", "dataname"))))

# 设置模型参数
knots <- 3

# --- 5. 初始化结果存储 ---
int.results.df <- data.frame()
significant_interactions <- list() # 使用列表存储更灵活

# --- 6. 主循环：执行 GAM 分析 ---
for (i in 1:nrow(EFvars.set)) {
  int_var <- EFvars.set$varname[i]
  dataname_str <- EFvars.set$dataname[i]
  dataname <- paste0(dataname_str, "_data")
  
  dependentvar <- "PHQ_y09"
  smooth_var <- "Age_year"
  covariates <- "Gender" # 性别作为协变量
  
  # 调用新的交互效应模型函数
  result <- gam.ordinal.interaction(
    dependentvar = "PHQ_y09", 
    dataname = dataname, 
    smooth_var = "Age_year", 
    int_var = int_var, 
    covariates = covariates,
    knots = knots,
    draws = 1000,
    bootstrap_sims = 1000 # 可以调整自举次数
  )
  
  if (!is.null(result)) {
    stats_df <- result$stats
    stats_df$dataname <- dataname_str
    int.results.df <- rbind(int.results.df, stats_df)
    
    # 【重要】现在直接使用 p_value_final 列来判断显著性
    if (stats_df$p_value_final < 0.05) {
      significant_interactions[[dataname_str]] <- list(
        stats = stats_df,
        effect_data = result$effect_data
      )
    }
  }
}

# --- 7. 保存结果 ---
if(nrow(int.results.df) > 0) {
  write.xlsx(int.results.df, file = paste0(resultFolder, "/PHQ_y09_interaction_results.xlsx"), row.names = FALSE)
}
if(length(significant_interactions) > 0) {
  saveRDS(significant_interactions, file = paste0(resultFolder, "/PHQ_y09_significant_interactions.rds"))
}

# --- 8. 可视化 ---

# 定义任务映射，用于美化标题
task_mapping <- c("GNGd" = "Go/No-Go", "back1" = "1-back", "back2" = "2-back")

# (在您的主分析脚本中)

# --- 8. 可视化 (新版本：生成四宫格图) ---

# --- 加载必要的包 ---
library(patchwork)

# 首先，确保结果数据框存在且不为空
if (exists("int.results.df") && nrow(int.results.df) > 0) {
  
  # --- 定义任务的美化标签 ---
  task_labels <- c(
    "GNGd" = "Go/No-Go",
    "back1" = "1-back",
    "back2" = "2-back"
  )
  
  # --- PART 1: 创建左上角的条形图 (图A) ---
  
  plot_data_A <- int.results.df %>%
    mutate(
      sig_label = ifelse(p_value_final < 0.05, "*", ""),
      task_label = factor(dataname, levels = c("GNGd", "back1", "back2"), labels = task_labels)
    )
  
  plot_A <- ggplot(plot_data_A, aes(x = task_label, y = pseudo_rsq, fill = task_label)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sig_label), vjust = -0.5, size = 8) +
    scale_fill_manual(values = c("Go/No-Go" = "#dbe6eb", "1-back" = "#a1c1d5", "2-back" = "#3a7ca5")) +
    scale_y_continuous(
      expression(paste("Interaction Effect (Pseudo R"^"2", ")")),
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.15))
    ) +
    labs(x = "EF Task", y = expression(paste("Effect Size (", Delta, "R"^"2", ")"))) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none", axis.title.x = element_text(size=14))
  
  # --- PART 2: 创建三个平滑的效应曲线图 (图B, C, D) ---
  
  # 初始化一个列表来存储三个曲线图
  plot_list_curves <- list()
  
  # 循环遍历每个任务
  for (task_name_short in c("GNGd", "back1", "back2")) {
    
    # 从结果中找到对应的行
    current_task_row <- int.results.df %>% filter(dataname == task_name_short)
    
    # 重新运行函数，只为了获取平滑的 effect_data_latent
    # (确保这里的参数与您的主分析循环一致)
    temp_result <- gam.ordinal.interaction(
      dependentvar = "PHQ_y09",
      dataname = paste0(task_name_short, "_data"),
      smooth_var = "Age_year",
      int_var = current_task_row$int_var,
      covariates = "Gender",
      knots = 3,
      draws = 1000
    )
    
    if (is.null(temp_result)) next # 如果出错则跳过
    
    plot_data_curve <- temp_result$effect_data_latent
    
    # 查找Y轴的对称范围，使图形更美观
    y_max <- max(abs(c(plot_data_curve$lower_ci, plot_data_curve$upper_ci))) * 1.1
    
    # 创建曲线图
    p_curve <- ggplot(plot_data_curve, aes(x = Age_year)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "#d62828") +
      geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, fill = "gray50") +
      geom_line(aes(y = median_effect), color = "black", linewidth = 1) +
      scale_y_continuous(limits = c(-y_max, y_max)) +
      labs(
        title = paste("PHQ-y09 (Latent) ~", task_labels[task_name_short]),
        x = "Age (years)",
        y = "Latent Effect (Slope)"
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(size = 12, hjust = 0.5),
        axis.title = element_text(size = 12)
      )
    
    # 将创建好的图存入列表
    plot_list_curves[[task_name_short]] <- p_curve
  }
  
  # --- PART 3: 使用 patchwork 组合四宫格图 ---
  
  if (length(plot_list_curves) == 3) {
    # 定义布局
    # (A | B)
    # (C | D)
    # 我们用一个空白图 plot_spacer() 来占位
    final_plot <- (plot_A | plot_list_curves[["GNGd"]]) / 
      (plot_list_curves[["back1"]] | plot_list_curves[["back2"]]) +
      plot_annotation(tag_levels = 'A')
    
    print(final_plot)
    
    # 保存最终的图片
    ggsave(
      filename = paste0(FigureFolder, "/PHQ_y09_Four_Panel_Latent_Effect_Plot.pdf"),
      plot = final_plot,
      width = 10, # 宽度 (英寸)
      height = 8, # 高度 (英寸)
      dpi = 300
    )
    cat("\n--- Four-panel plot saved successfully. ---\n")
    
  } else {
    print("--- WARNING: Could not generate all three curve plots. Final panel plot not created. ---")
  }
  
} else {
  print("--- ERROR: The result dataframe 'int.results.df' is empty. Cannot create plots. ---")
}