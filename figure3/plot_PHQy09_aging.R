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

GNGd_data <- readRDS(paste0(interfileFolder, '/GNGd_prime.deviations.rds'))
back1_data <- readRDS(paste0(interfileFolder, '/back1Acc.deviations.rds'))
back2_data <- readRDS(paste0(interfileFolder, '/back2Acc.deviations.rds'))

PHQ_data <- read_xlsx(paste0(interfileFolder, '/PHQ.xlsx')) %>%
  dplyr::select(用户ID, PHQ_y09) %>%  
  dplyr::mutate(用户ID = as.character(用户ID))


GNGd_data <- GNGd_data %>% dplyr::mutate(ID = as.character(x__ID), Gender = as.factor(Sex))
back1_data <- back1_data %>% dplyr::mutate(ID = as.character(x__ID), Gender = as.factor(Sex))
back2_data <- back2_data %>% dplyr::mutate(ID = as.character(x__ID), Gender = as.factor(Sex))

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
  
  # 调用我们更新后的函数
  result <- gam.ordinal.interaction(
    dependentvar = "PHQ_y09", 
    dataname = dataname, 
    smooth_var = "Age_year", 
    int_var = int_var, 
    covariates = covariates,
    knots = knots,
    draws = 1000,
    bootstrap_sims = 1000 
  )
  
  if (!is.null(result)) {
    stats_df <- result$stats
    stats_df$dataname <- dataname_str
    
    # 这一步总会执行，为绘图准备数据
    int.results.df <- rbind(int.results.df, stats_df)
    
    # --- 【修改点】---
    # 1. 将 p_value_final 改为 interaction_p
    # 2. 检查 p 值是否存在且不为 NA，这更稳健
    if (!is.na(stats_df$interaction_p) && stats_df$interaction_p < 0.05) {
      
      # 额外保存显著的交互作用结果
      significant_interactions[[dataname_str]] <- list(
        stats = stats_df,
        # 明确指定保存 latent space 的曲线数据
        effect_data_latent = result$effect_data_latent 
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

# --- 8. 可视化 (最终修改版) ---

# --- 加载必要的包 ---
library(patchwork)

if (exists("int.results.df") && nrow(int.results.df) > 0) {
  task_labels <- c(
    "GNGd" = "Go/No-Go",
    "back1" = "1-back",
    "back2" = "2-back"
  )
  
  # --- PART 1: 创建左上角的条形图 (图A) - 显示主效应的 Beta 值和其显著性 ---
  
  plot_data_A <- int.results.df %>%
    mutate(
      # --- 修改：星号现在基于 main_effect_p ---
      sig_label = ifelse(main_effect_p < 0.05, "*", ""), 
      task_label = factor(dataname, levels = c("GNGd", "back1", "back2"), labels = task_labels)
    )
  
  plot_A <- ggplot(plot_data_A, aes(x = task_label, y = main_effect_beta, fill = task_label)) +
    geom_col(width = 0.7) +
    # --- 星号代表主效应的显著性 (即柱子高度的显著性) ---
    geom_text(aes(label = sig_label), vjust = -0.5, size = 8, color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
    scale_fill_manual(values = c("Go/No-Go" = "#dbe6eb", "1-back" = "#a1c1d5", "2-back" = "#3a7ca5")) +
    scale_y_continuous(
      expand = expansion(mult = c(0.1, 0.15))
    ) +
    labs(
      x = "EF Task", 
      y = expression(paste("Main Effect (", beta, " Coefficient)")) 
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none", axis.title.x = element_text(size=14))
  
  # --- PART 2: 创建三个平滑的效应曲线图 (图B, C, D) ---
  
  # 初始化一个列表来存储三个曲线图
  plot_list_curves <- list()
  
  # 循环遍历每个任务 (这部分代码无需修改)
  for (task_name_short in c("GNGd", "back1", "back2")) {
    
    current_task_row <- int.results.df %>% filter(dataname == task_name_short)
    
    if (nrow(current_task_row) == 0 || is.na(current_task_row$main_effect_beta)) next
    
    beta_value <- current_task_row$main_effect_beta
    
    # 我们需要重新运行模型以获取绘图数据
    temp_result <- gam.ordinal.interaction(
      dependentvar = "PHQ_y09",
      dataname = paste0(task_name_short, "_data"),
      smooth_var = "Age_year",
      int_var = current_task_row$int_var,
      covariates = "Gender",
      knots = 3,
      draws = 1000
    )
    
    if (is.null(temp_result)) next
    
    plot_data_curve <- temp_result$effect_data_latent
    
    y_range <- range(c(plot_data_curve$lower_ci, plot_data_curve$upper_ci, beta_value), na.rm = TRUE)
    y_buffer <- (y_range[2] - y_range[1]) * 0.1
    y_lims <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)
    
    p_curve <- ggplot(plot_data_curve, aes(x = Age_year)) +
      geom_hline(yintercept = beta_value, linetype = "dashed", color = "#d62828", linewidth = 1) +
      geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, fill = "gray50") +
      geom_line(aes(y = median_effect), color = "black", linewidth = 1) +
      scale_y_continuous(limits = y_lims) +
      labs(
        title = task_labels[task_name_short],
        x = "Age (years)",
        y = "Interaction Effect (Slope)"
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(size = 12, hjust = 0.5),
        axis.title = element_text(size = 12)
      )
    
    plot_list_curves[[task_name_short]] <- p_curve
  }
  
  # --- PART 3: 使用 patchwork 组合四宫格图 (无需修改) ---
  
  if (length(plot_list_curves) == 3) {
    final_plot <- (plot_A | plot_list_curves[["GNGd"]]) / 
      (plot_list_curves[["back1"]] | plot_list_curves[["back2"]]) +
      plot_annotation(tag_levels = 'A')
    
    print(final_plot)
    
    ggsave(
      filename = paste0(FigureFolder, "/PHQ_y09_Four_Panel_Main_Effect_and_Interaction(Main eff)_Plot.pdf"),
      plot = final_plot,
      width = 10,
      height = 8,
      dpi = 300
    )
    cat("\n--- Four-panel plot saved successfully. ---\n")
    
  } else {
    print("--- WARNING: Could not generate all three curve plots. Final panel plot not created. ---")
  }
  
} else {
  print("--- ERROR: The result dataframe 'int.results.df' is empty. Cannot create plots. ---")
}