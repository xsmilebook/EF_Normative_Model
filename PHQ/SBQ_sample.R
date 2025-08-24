# 1. 加载所需的库
# 确保已安装所有必要的包: install.packages(c("dplyr", "openxlsx", "ggplot2", "patchwork", "scales"))
library(dplyr)
library(openxlsx)
library(ggplot2)
library(patchwork) # 用于拼接和嵌入图表
library(scales)    # 用于格式化百分比

# 2. 读取和准备数据 (与您之前的代码相同)
datapath <- "D:/datasets/yunfu/interfile_folder/deviationZ_temp"
sbq_path <- "D:/datasets/yunfu/raw_data/2.0_云浮队列基线学生问卷_自杀相关变量.xlsx"

if (!dir.exists(datapath) || !file.exists(sbq_path)) {
  stop("一个或多个文件路径不正确，请检查datapath和sbq_path变量。")
}

oneback_data <- readRDS(file.path(datapath, "back1Acc.deviations.rds"))
twoback_data <- readRDS(file.path(datapath, "back2Acc.deviations.rds"))
GNG_data <- readRDS(file.path(datapath, "GNGd_prime.deviations.rds"))
SBQ_data <- read_xlsx(sbq_path)

SBQ_data$"用户ID" <- as.numeric(SBQ_data$"用户ID")
SBQ_data_renamed <- SBQ_data %>% rename(x__ID = "用户ID")

process_data <- function(task_data, sbq_data) {
  task_data %>%
    inner_join(sbq_data, by = "x__ID") %>%
    filter(
      SuicideAttempt_Lifetime %in% c(0, 1),
      SuicideAttempt_Pastyear %in% c(0, 1)
    )
}

oneback_final <- process_data(oneback_data, SBQ_data_renamed)
twoback_final <- process_data(twoback_data, SBQ_data_renamed)
gng_final     <- process_data(GNG_data, SBQ_data_renamed)


# 3. 定义核心分析函数：执行100次重采样t检验 (已修正)
run_resampling_ttest <- function(data, test_var, group_var, n_iterations = 100) {
  
  group_0_data <- data %>% filter(.data[[group_var]] == 0)
  group_1_data <- data %>% filter(.data[[group_var]] == 1)
  
  n_small <- nrow(group_1_data)
  
  if (n_small < 3) {
    return(NULL)
  }
  
  results <- replicate(n_iterations, {
    group_0_sample <- group_0_data %>% sample_n(n_small)
    resampled_data <- bind_rows(group_0_sample, group_1_data)
    
    formula <- as.formula(paste(test_var, "~", group_var))
    t_test_result <- t.test(formula, data = resampled_data)
    
    c(t_stat = t_test_result$statistic, p_value = t_test_result$p.value)
  })
  
  # --- 这是修正过的部分 ---
  # 将结果转换为整洁的数据框，并明确重命名，确保后续操作无误
  df_results <- as.data.frame(t(results))
  colnames(df_results) <- c("t_stat", "p_value")
  return(df_results)
}


# 4. 定义新的绘图函数：可视化重采样结果 (无需改动)
plot_resampling_results <- function(original_data, resampled_stats, test_var, group_var) {
  
  plot_title <- paste(test_var, "by", group_var, "\n(Results from 100 resamples)")
  
  mean_t <- mean(resampled_stats$t_stat)
  median_p <- median(resampled_stats$p_value)
  prop_sig <- mean(resampled_stats$p_value < 0.05)
  
  stats_label <- paste0(
    "Avg. t-statistic = ", sprintf("%.2f", mean_t),
    "\nMedian p-value = ", sprintf("%.3f", median_p),
    "\n% Significant (p<0.05) = ", percent(prop_sig, accuracy = 1)
  )
  
  original_data[[group_var]] <- as.factor(original_data[[group_var]])
  sample_sizes <- original_data %>% group_by(.data[[group_var]]) %>% summarise(n = n())
  axis_labels <- paste0(sample_sizes[[group_var]], "\n(n=", sample_sizes$n, ")")
  
  y_range <- range(original_data[[test_var]], na.rm = TRUE)
  y_pos <- y_range[1] + 0.15 * (y_range[2] - y_range[1])
  
  main_plot <- ggplot(original_data, aes_string(x = group_var, y = test_var)) +
    geom_jitter(width = 0.2, height = 0, color = "lightblue", alpha = 0.3) +
    geom_boxplot(width = 0.5, alpha = 0.6, outlier.shape = NA, aes_string(fill = group_var)) +
    annotate("text", x = 1.5, y = y_pos, label = stats_label, size = 3.5, lineheight = 1.2, hjust = 0.5) +
    scale_x_discrete(labels = axis_labels) +
    labs(title = plot_title, x = group_var, y = test_var) +
    theme_bw(base_size = 12) + 
    theme(legend.position = "none", plot.title = element_text(size = 10))
  
  hist_t <- ggplot(resampled_stats, aes(x = t_stat)) +
    geom_histogram(bins = 20, fill = "grey70", color = "black") +
    geom_vline(xintercept = mean_t, color = "red", linetype = "dashed") +
    labs(title = "t-statistic Distribution", x = NULL, y = NULL) +
    theme_minimal(base_size = 8)
  
  hist_p <- ggplot(resampled_stats, aes(x = p_value)) +
    geom_histogram(bins = 20, fill = "grey70", color = "black") +
    geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") +
    labs(title = "p-value Distribution", x = NULL, y = NULL) +
    theme_minimal(base_size = 8)
  
  final_composite_plot <- main_plot + 
    inset_element(hist_t, left = 0.02, bottom = 0.6, right = 0.45, top = 0.98) +
    inset_element(hist_p, left = 0.55, bottom = 0.6, right = 0.98, top = 0.98)
  
  return(final_composite_plot)
}

# 5. 执行所有分析并生成图表
all_plots <- list()
analyses_to_run <- list(
  list(data = oneback_final, test = "Oneback_acc_deviationZ", group = "SuicideAttempt_Lifetime"),
  list(data = oneback_final, test = "Oneback_acc_deviationZ", group = "SuicideAttempt_Pastyear"),
  list(data = twoback_final, test = "Twoback_acc_deviationZ", group = "SuicideAttempt_Lifetime"),
  list(data = twoback_final, test = "Twoback_acc_deviationZ", group = "SuicideAttempt_Pastyear"),
  list(data = gng_final,     test = "d_prime_deviationZ",   group = "SuicideAttempt_Lifetime"),
  list(data = gng_final,     test = "d_prime_deviationZ",   group = "SuicideAttempt_Pastyear")
)

for (i in 1:length(analyses_to_run)) {
  params <- analyses_to_run[[i]]
  cat("Running analysis for:", params$test, "by", params$group, "...\n")
  
  resampling_results <- run_resampling_ttest(params$data, params$test, params$group)
  
  if (!is.null(resampling_results)) {
    plot_name <- paste0("p", i)
    all_plots[[plot_name]] <- plot_resampling_results(params$data, resampling_results, params$test, params$group)
  }
}

# 6. 拼接并保存最终的PDF
if (length(all_plots) > 0) {
  final_plot <- patchwork::wrap_plots(all_plots, ncol = 2) +
    plot_annotation(
      title = "Cognitive Performance by Suicide History (Assessed with 100x Resampling)",
      tag_levels = 'A' 
    )
  
  output_dir <- "D:/datasets/yunfu/results/SBQ"
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  output_pdf_path <- file.path(output_dir, "cognitive_performance_resampling_analysis.pdf")
  ggsave(output_pdf_path, final_plot, width = 14, height = 15, units = "in", device = "pdf")
  
  cat(paste("\n所有图表已成功拼接并保存到PDF文件:\n", output_pdf_path, "\n"))
  
} else {
  cat("\n没有生成任何图表，无法创建PDF。\n")
}