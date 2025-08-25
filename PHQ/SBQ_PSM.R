# 1. 加载所需的库
# 确保已安装所有必要的包
library(dplyr)
library(openxlsx)
library(ggplot2)
library(patchwork) 
library(MatchIt)   # 用于倾向性评分匹配
library(cobalt)    # 用于评估匹配后的平衡性

# 2. 读取和准备数据
datapath <- "D:/datasets/yunfu/interfile_folder/deviationZ_temp"
sbq_path <- "D:/datasets/yunfu/raw_data/2.0_云浮队列基线学生问卷_自杀相关变量.xlsx"

oneback_data <- readRDS(file.path(datapath, "back1Acc.deviations.rds"))
twoback_data <- readRDS(file.path(datapath, "back2Acc.deviations.rds"))
GNG_data <- readRDS(file.path(datapath, "GNGd_prime.deviations.rds"))
SBQ_data <- read_xlsx(sbq_path)

# 3. 数据预处理和合并
SBQ_data_prepared <- SBQ_data %>%
  rename(x__ID = "用户ID") %>%
  mutate(x__ID = as.numeric(x__ID))

# --- 修改点在这里 ---
process_and_merge_data <- function(task_data, sbq_data) {
  task_data %>%
    mutate(
      Gender = as.factor(Gender),
      Age_year = as.numeric(Age_year),
      School = as.factor(School)
    ) %>%
    inner_join(sbq_data, by = "x__ID") %>%
    # 修改筛选逻辑：只要任意一个自杀未遂史变量有效，就保留该行
    filter(
      (SuicideAttempt_Lifetime %in% c(0, 1)) | (SuicideAttempt_Pastyear %in% c(0, 1))
    ) %>%
    # 协变量的NA值仍然需要被过滤
    filter(
      !is.na(Gender) & !is.na(Age_year) & !is.na(School)
    )
}

oneback_final <- process_and_merge_data(oneback_data, SBQ_data_prepared)
twoback_final <- process_and_merge_data(twoback_data, SBQ_data_prepared)
gng_final     <- process_and_merge_data(GNG_data, SBQ_data_prepared)


# 4. 定义PSM分析与绘图主函数 (无需修改)
run_psm_analysis_and_plot <- function(data, test_var, group_var) {
  
  # 注意：这里的filter会为当前分析移除在group_var上为NA的行
  current_data <- data %>%
    filter(!is.na(.data[[test_var]]) & !is.na(.data[[group_var]])) %>%
    mutate(!!sym(group_var) := as.numeric(as.character(.data[[group_var]])))
  
  # ... 函数其余部分与之前相同 ...
  min_group_size <- min(table(current_data[[group_var]]), na.rm = TRUE)
  if (min_group_size < 5) {
    cat(paste("Skipping PSM for", group_var, "on", test_var, "due to small minority group size (n =", min_group_size, ").\n"))
    return(NULL)
  }
  
  match_formula <- as.formula(paste(group_var, "~ Gender + Age_year + School"))
  
  match_obj <- matchit(
    match_formula,
    data = current_data,
    method = "nearest", 
    ratio = 1,          
    caliper = 0.2
  )
  
  matched_data <- match.data(match_obj)
  
  bal_plot <- love.plot(
    match_obj, 
    binary = "std", 
    thresholds = c(m = .1),
    title = paste("Covariate Balance for", group_var),
  )
  
  stats_label <- tryCatch({
    treatment_group <- matched_data %>% filter(.data[[group_var]] == 1) %>% arrange(subclass)
    control_group <- matched_data %>% filter(.data[[group_var]] == 0) %>% arrange(subclass)
    
    if (nrow(treatment_group) != nrow(control_group)) { stop("Matched groups have different sizes.") }
    
    t_test_result <- t.test(treatment_group[[test_var]], control_group[[test_var]], paired = TRUE)
    
    t_val_formatted <- sprintf("%.2f", t_test_result$statistic)
    p_val_formatted <- if (t_test_result$p.value < 0.001) "p < 0.001" else paste("p =", sprintf("%.3f", t_test_result$p.value))
    
    paste0("Paired t = ", t_val_formatted, ", ", p_val_formatted)
    
  }, error = function(e) { 
    print(e$message)
    "Paired t-test failed" 
  })
  
  matched_data[[group_var]] <- as.factor(matched_data[[group_var]])
  plot_title <- paste(test_var, "by", group_var, "\n(After 1:1 PSM)")
  
  sample_sizes <- matched_data %>% group_by(.data[[group_var]]) %>% summarise(n = n())
  axis_labels <- paste0(levels(as.factor(sample_sizes[[group_var]])), "\n(n=", sample_sizes$n, ")")
  
  result_plot <- ggplot(matched_data, aes_string(x = group_var, y = test_var)) +
    geom_jitter(width = 0.2, height = 0, color = "lightblue", alpha = 0.5) +
    geom_boxplot(width = 0.5, alpha = 0.6, outlier.shape = NA, aes_string(fill = group_var)) +
    geom_line(aes(group = subclass), color = "grey", alpha = 0.4) + 
    annotate("text", x = 1.5, y = Inf, vjust = 1.2, label = stats_label, size = 4) +
    scale_x_discrete(labels = axis_labels) +
    labs(title = plot_title, x = group_var, y = test_var) +
    theme_bw(base_size = 12) + 
    theme(legend.position = "none", plot.title = element_text(size = 10))
  
  return(list(result_plot = result_plot, balance_plot = bal_plot))
}


# 5. 定义要执行的所有分析 (无需修改)
analyses_to_run <- list(
  list(data = oneback_final, test = "Oneback_acc_deviationZ", group = "SuicideAttempt_Lifetime"),
  list(data = oneback_final, test = "Oneback_acc_deviationZ", group = "SuicideAttempt_Pastyear"),
  list(data = twoback_final, test = "Twoback_acc_deviationZ", group = "SuicideAttempt_Lifetime"),
  list(data = twoback_final, test = "Twoback_acc_deviationZ", group = "SuicideAttempt_Pastyear"),
  list(data = gng_final,     test = "d_prime_deviationZ",   group = "SuicideAttempt_Lifetime"),
  list(data = gng_final,     test = "d_prime_deviationZ",   group = "SuicideAttempt_Pastyear")
)

# 6. 循环执行所有分析并收集图表 (无需修改)
all_result_plots <- list()
for (i in 1:length(analyses_to_run)) {
  params <- analyses_to_run[[i]]
  plot_name <- paste(params$test, params$group, sep = "_by_")
  
  cat("=====================================================\n")
  cat("Running PSM analysis for:", params$test, "by", params$group, "...\n")
  
  psm_results <- run_psm_analysis_and_plot(params$data, params$test, params$group)
  
  if (!is.null(psm_results)) {
    all_result_plots[[plot_name]] <- psm_results$result_plot
  }
}

# 7. 按照 GNG, Oneback, Twoback 的顺序重新排列图表并拼接 (无需修改)
if (length(all_result_plots) > 0) {
  gng_lifetime_plot <- all_result_plots[["d_prime_deviationZ_by_SuicideAttempt_Lifetime"]]
  oneback_lifetime_plot <- all_result_plots[["Oneback_acc_deviationZ_by_SuicideAttempt_Lifetime"]]
  twoback_lifetime_plot <- all_result_plots[["Twoback_acc_deviationZ_by_SuicideAttempt_Lifetime"]]
  
  gng_pastyear_plot <- all_result_plots[["d_prime_deviationZ_by_SuicideAttempt_Pastyear"]]
  oneback_pastyear_plot <- all_result_plots[["Oneback_acc_deviationZ_by_SuicideAttempt_Pastyear"]]
  twoback_pastyear_plot <- all_result_plots[["Twoback_acc_deviationZ_by_SuicideAttempt_Pastyear"]]
  
  ordered_plots <- list(
    gng_lifetime_plot, oneback_lifetime_plot, twoback_lifetime_plot,
    gng_pastyear_plot, oneback_pastyear_plot, twoback_pastyear_plot
  )
  
  final_plot <- patchwork::wrap_plots(ordered_plots, ncol = 3, nrow = 2) +
    plot_annotation(
      title = "Cognitive Performance by Suicide Attempt History (After PSM)",
      subtitle = "Row 1: Lifetime Suicide Attempt | Row 2: Past-Year Suicide Attempt",
      tag_levels = 'A'
    )
  
  output_dir <- "D:/datasets/yunfu/results/SBQ"
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  output_pdf_path <- file.path(output_dir, "cognitive_performance_PSM_2x3_layout.pdf")
  ggsave(output_pdf_path, final_plot, width = 18, height = 10, units = "in", device = "pdf", limitsize = FALSE)
  
  cat(paste("\n所有分析图表已按 2x3 布局成功拼接并保存到PDF文件:\n", output_pdf_path, "\n"))
  
} else {
  cat("\n没有生成任何图表，无法创建PDF。\n")
}