# ===================================================================
#           完整且可直接运行的 R 脚本 (合并六条曲线于一图)
# ===================================================================

# --- 0. 加载所有必要的 R 包 ---
# 如果您尚未安装这些包，请先运行 install.packages("...")
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(forcats)
  library(ggplot2)
  library(mgcv)
  library(purrr)
  library(patchwork)
  library(readxl)
})

# --- 1. 文件路径设置 ---
# !!! 请确保这里的路径是正确的 !!!
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  # 如果在服务器上运行，请取消注释并修改以下路径
  # interfileFolder <- "/path/to/your/interfile_folder"
  # functionFolder <- "/path/to/your/functions"
  # resultFolder <- "/path/to/your/results/phq_y09_interaction_with_main_effect"
  # FigureFolder <- '/path/to/your/figures/phq_y09_fig_with_main_effect'
} else {
  # 本地计算机的路径
  interfileFolder <- "D:/datasets/yunfu/interfile_folder/temp"
  functionFolder <- "D:/code/EF_Normative_Model/functions"
  resultFolder <- "D:/datasets/yunfu/results/phq_y09_interaction_with_main_effect"
  FigureFolder <- 'D:/datasets/yunfu/figures/phq_y09_fig_with_main_effect'
}
# 创建输出文件夹（如果尚不存在）
dir.create(resultFolder, showWarnings = FALSE)
dir.create(FigureFolder, showWarnings = FALSE)


# --- 2. 数据读取与准备 ---
GNGd_data <- readRDS(paste0(interfileFolder, '/GNGd_prime.deviations.rds'))
back1_data <- readRDS(paste0(interfileFolder, '/back1Acc.deviations.rds'))
back2_data <- readRDS(paste0(interfileFolder, '/back2Acc.deviations.rds'))

PHQ_data <- read_xlsx(paste0(interfileFolder, '/PHQ.xlsx')) %>%
  dplyr::select(用户ID, PHQ_y09) %>%
  dplyr::mutate(用户ID = as.character(用户ID))

GNGd_data <- GNGd_data %>% dplyr::mutate(ID = as.character(x__ID), Sex = as.factor(Sex)) %>% dplyr::inner_join(PHQ_data, by = c("ID" = "用户ID"))
back1_data <- back1_data %>% dplyr::mutate(ID = as.character(x__ID), Sex = as.factor(Sex)) %>% dplyr::inner_join(PHQ_data, by = c("ID" = "用户ID"))
back2_data <- back2_data %>% dplyr::mutate(ID = as.character(x__ID), Sex = as.factor(Sex)) %>% dplyr::inner_join(PHQ_data, by = c("ID" = "用户ID"))

df_all <- dplyr::bind_rows(
  GNGd_data %>% dplyr::select(ID, Age_year, Sex, PHQ_y09, deviationZ = d_prime_deviationZ) %>% mutate(Task = "GNGd"),
  back1_data %>% dplyr::select(ID, Age_year, Sex, PHQ_y09, deviationZ = Oneback_acc_deviationZ) %>% mutate(Task = "back1"),
  back2_data %>% dplyr::select(ID, Age_year, Sex, PHQ_y09, deviationZ = Twoback_acc_deviationZ) %>% mutate(Task = "back2")
)

df_all$Task <- factor(df_all$Task, levels = c("GNGd", "back1", "back2"))
df_all$PHQ_y09 <- factor(df_all$PHQ_y09, levels = c("0", "1", "2", "3"), ordered = FALSE)


# --- 3. 分析函数 (无需修改) ---

# 函数 1: 计算与基线 "0" 的对比
analyze_task_effects <- function(task_name, data) {
  cat(sprintf("\n--- Analyzing Baseline Contrasts for Task: %s ---\n", task_name))
  dft <- data %>% filter(Task == task_name, !is.na(deviationZ), !is.na(Age_year), !is.na(PHQ_y09), !is.na(Sex))
  if(nrow(dft) < 50 || length(unique(dft$Age_year)) < 4) return(NULL)
  k_use <- 3
  formula_interaction <- deviationZ ~ Sex + s(Age_year, k = k_use, fx = F) + s(Age_year, by = PHQ_y09, k = k_use, fx = F) + PHQ_y09
  fit_interaction <- try(gam(formula_interaction, data = dft, method = "REML"), silent = TRUE)
  if (inherits(fit_interaction, "try-error")) return(NULL)
  ages <- seq(min(dft$Age_year), max(dft$Age_year), length.out = 100)
  sex_proportions <- table(dft$Sex) / nrow(dft)
  make_newdata <- function(phq_level, sex_prop) { map_dfr(names(sex_prop), ~data.frame(Age_year = ages, PHQ_y09 = factor(phq_level, levels = levels(dft$PHQ_y09)), Sex = factor(.x, levels = levels(dft$Sex)), weight = sex_prop[.x])) }
  nd0 <- make_newdata("0", sex_proportions); Xp0 <- predict(fit_interaction, newdata = nd0, type = "lpmatrix")
  effect_curves <- map_dfr(c("1", "2", "3"), function(level) {
    if (!level %in% levels(dft$PHQ_y09)) return(NULL)
    nd_level <- make_newdata(level, sex_proportions); Xp_level <- predict(fit_interaction, newdata = nd_level, type = "lpmatrix")
    D <- (Xp_level - Xp0) * nd_level$weight; D_with_age <- as.data.frame(D); D_with_age$Age_year <- nd_level$Age_year
    D_agg_df <- D_with_age %>% group_by(Age_year) %>% summarise(across(everything(), sum)) %>% ungroup()
    D_agg_matrix <- as.matrix(D_agg_df[, -1]); beta <- coef(fit_interaction); V <- vcov(fit_interaction, unconditional = TRUE)
    effect <- as.numeric(D_agg_matrix %*% beta); se <- sqrt(rowSums((D_agg_matrix %*% V) * D_agg_matrix))
    data.frame(Age_year = D_agg_df$Age_year, Contrast = paste("PHQ", level, "vs 0"), effect_size = effect, lower_ci = effect - 1.96 * se, upper_ci = effect + 1.96 * se)
  })
  formula_main <- deviationZ ~ Sex + s(Age_year, k = k_use, fx = F) + PHQ_y09
  fit_main <- try(gam(formula_main, data = dft, method = "REML"), silent = TRUE)
  main_effects_df <- NULL
  if (!inherits(fit_main, "try-error")) {
    p_table <- summary(fit_main)$p.table; phq_rows <- rownames(p_table)[str_detect(rownames(p_table), "PHQ_y09")]
    if (length(phq_rows) > 0) { main_effects_df <- tibble(Contrast = paste("PHQ", str_extract(phq_rows, "\\d$"), "vs 0"), effect_size = p_table[phq_rows, "Estimate"]) }
  }
  return(list(effect_curves = effect_curves, main_effects = main_effects_df))
}

# 函数 2: 计算两两对比
analyze_pairwise_effects <- function(task_name, data) {
  cat(sprintf("\n--- Analyzing Pairwise Contrasts for Task: %s ---\n", task_name))
  dft <- data %>% filter(Task == task_name, !is.na(deviationZ), !is.na(Age_year), !is.na(PHQ_y09), !is.na(Sex))
  if(nrow(dft) < 50 || length(unique(dft$Age_year)) < 4) return(NULL)
  k_use <- 3
  formula_interaction <- deviationZ ~ Sex + s(Age_year, k = k_use, fx = F) + s(Age_year, by = PHQ_y09, k = k_use, fx = F)
  fit_interaction <- try(gam(formula_interaction, data = dft, method = "REML"), silent = TRUE)
  if (inherits(fit_interaction, "try-error")) return(NULL)
  ages <- seq(min(dft$Age_year), max(dft$Age_year), length.out = 100)
  sex_proportions <- table(dft$Sex) / nrow(dft)
  make_newdata <- function(phq_level, sex_prop) { map_dfr(names(sex_prop), ~data.frame(Age_year = ages, PHQ_y09 = factor(phq_level, levels = levels(dft$PHQ_y09)), Sex = factor(.x, levels = levels(dft$Sex)), weight = sex_prop[.x])) }
  contrasts_list <- list(c("2", "1"), c("3", "1"), c("3", "2"))
  effect_curves_pairwise <- map_dfr(contrasts_list, function(pair) {
    level1 <- pair[1]; level2 <- pair[2]
    if (!level1 %in% levels(dft$PHQ_y09) || !level2 %in% levels(dft$PHQ_y09)) return(NULL)
    nd1 <- make_newdata(level1, sex_proportions); Xp1 <- predict(fit_interaction, newdata = nd1, type = "lpmatrix")
    nd2 <- make_newdata(level2, sex_proportions); Xp2 <- predict(fit_interaction, newdata = nd2, type = "lpmatrix")
    D <- (Xp1 - Xp2) * nd1$weight; D_with_age <- as.data.frame(D); D_with_age$Age_year <- nd1$Age_year
    D_agg_df <- D_with_age %>% group_by(Age_year) %>% summarise(across(everything(), sum)) %>% ungroup()
    D_agg_matrix <- as.matrix(D_agg_df[, -1]); beta <- coef(fit_interaction); V <- vcov(fit_interaction, unconditional = TRUE)
    effect <- as.numeric(D_agg_matrix %*% beta); se <- sqrt(rowSums((D_agg_matrix %*% V) * D_agg_matrix))
    data.frame(Age_year = D_agg_df$Age_year, Contrast = paste("PHQ", level1, "vs", level2), effect_size = effect, lower_ci = effect - 1.96 * se, upper_ci = effect + 1.96 * se)
  })
  formula_main <- deviationZ ~ Sex + s(Age_year, k = k_use, fx = F) + PHQ_y09
  fit_main <- try(gam(formula_main, data = dft, method = "REML"), silent = TRUE)
  main_effects_pairwise_df <- NULL
  if (!inherits(fit_main, "try-error")) {
    p_table <- summary(fit_main)$p.table
    required_coeffs <- c("PHQ_y091", "PHQ_y092", "PHQ_y093")
    if (all(required_coeffs %in% rownames(p_table))) {
      est <- p_table[, "Estimate"]
      main_effects_pairwise_df <- tibble(
        Contrast = c("PHQ 2 vs 1", "PHQ 3 vs 1", "PHQ 3 vs 2"),
        effect_size = c(est["PHQ_y092"] - est["PHQ_y091"], est["PHQ_y093"] - est["PHQ_y091"], est["PHQ_y093"] - est["PHQ_y092"])
      )
    }
  }
  return(list(effect_curves = effect_curves_pairwise, main_effects = main_effects_pairwise_df))
}


# --- 4. 运行所有分析 ---
tasks <- c("GNGd", "back1", "back2")
# 分析 1: vs Baseline
baseline_results <- map(tasks, ~analyze_task_effects(.x, data = df_all)) %>% set_names(tasks)
# 分析 2: Pairwise
pairwise_results <- map(tasks, ~analyze_pairwise_effects(.x, data = df_all)) %>% set_names(tasks)


# --- 5. 合并数据并进行最终的可视化 ---
# 清理空结果
baseline_results_clean <- baseline_results[!sapply(baseline_results, is.null)]
pairwise_results_clean <- pairwise_results[!sapply(pairwise_results, is.null)]

if (length(baseline_results_clean) > 0 && length(pairwise_results_clean) > 0) {
  
  # 准备数据框
  plot_data_curves_baseline <- map_dfr(names(baseline_results_clean), ~baseline_results_clean[[.x]]$effect_curves %>% mutate(Task = .x))
  plot_data_main_baseline <- map_dfr(names(baseline_results_clean), ~baseline_results_clean[[.x]]$main_effects %>% mutate(Task = .x))
  
  plot_data_curves_pairwise <- map_dfr(names(pairwise_results_clean), ~pairwise_results_clean[[.x]]$effect_curves %>% mutate(Task = .x))
  plot_data_main_pairwise <- map_dfr(names(pairwise_results_clean), ~pairwise_results_clean[[.x]]$main_effects %>% mutate(Task = .x))
  
  # *** 合并所有绘图数据 ***
  plot_data_curves_combined <- bind_rows(plot_data_curves_baseline, plot_data_curves_pairwise)
  plot_data_main_combined <- bind_rows(plot_data_main_baseline, plot_data_main_pairwise)
  
  # *** 定义包含所有6个对比的标签和颜色 ***
  task_labels <- c("GNGd" = "Go/No-Go", "back1" = "1-back", "back2" = "2-back")
  
  combined_contrast_colors <- c(
    # 基线对比 (蓝色系)
    "PHQ 1 vs 0" = "#6baed6", 
    "PHQ 2 vs 0" = "#3182bd", 
    "PHQ 3 vs 0" = "#08519c", 
    # 两两对比 (红色/紫色系)
    "PHQ 2 vs 1" = "#2ca02c", # 绿色
    "PHQ 3 vs 1" = "#d62728", # 红色
    "PHQ 3 vs 2" = "#9467bd"  # 紫色
  )
  
  # 确保所有绘图数据的因子顺序和标签正确
  plot_data_curves_combined$Task <- factor(plot_data_curves_combined$Task, levels = names(task_labels), labels = task_labels)
  plot_data_curves_combined$Contrast <- factor(plot_data_curves_combined$Contrast, levels = names(combined_contrast_colors))
  
  if(nrow(plot_data_main_combined) > 0){
    plot_data_main_combined$Task <- factor(plot_data_main_combined$Task, levels = names(task_labels), labels = task_labels)
    plot_data_main_combined$Contrast <- factor(plot_data_main_combined$Contrast, levels = names(combined_contrast_colors))
  }
  
  # *** 创建最终的组合图 ***
  combined_plot <- ggplot(plot_data_curves_combined, aes(x = Age_year, y = effect_size, color = Contrast, fill = Contrast)) +
    # 0线
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    # 交互效应置信区间
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.15, linetype = 0) +
    # 交互效应曲线
    geom_line(linewidth = 1) +
    # 主效应基准线
    geom_hline(
      data = plot_data_main_combined,
      aes(yintercept = effect_size, color = Contrast),
      linetype = "dotted",
      linewidth = 1
    ) +
    # 按任务分面
    facet_wrap(~ Task, ncol = 3) +
    # 视觉元素设置
    scale_color_manual(values = combined_contrast_colors, name = "Contrast", drop = FALSE) +
    scale_fill_manual(values = combined_contrast_colors, name = "Contrast", drop = FALSE) +
    labs(
      title = "All PHQ Contrasts on Cognitive Deviation Across Age",
      subtitle = "Solid lines are interaction effects; Dotted lines are main effects.",
      x = "Age (years)",
      y = "Effect Size (Difference in Deviation Z-score)"
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face="bold", size=16),
      plot.subtitle = element_text(size=12, color="grey30"),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold")
    ) +
    guides(color = guide_legend(nrow = 2)) # 将图例分为两行，使其更紧凑
  
  # 在R环境中显示图片
  print(combined_plot)
  
  # 保存图片
  ggsave(
    filename = file.path(FigureFolder, "PHQ_y09_All_Six_Contrasts_Combined.pdf"),
    plot = combined_plot,
    width = 13, # 稍微加宽以容纳图例
    height = 6,
    dpi = 300
  )
  
  cat(sprintf("\n--- Combined analysis complete. Plot saved to %s ---\n", FigureFolder))
  
} else {
  cat("\n--- One or more analysis sets failed. Could not generate combined plot. ---\n")
}