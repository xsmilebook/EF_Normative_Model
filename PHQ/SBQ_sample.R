# 1. 加载所需的库 (新增 stringr 用于备用方案，但本方案未使用)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(patchwork)
library(scales)
# library(stringr) # 如果需要自动换行，则取消此行注释

# 2. 读取和准备数据 (无变化)
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
      SuicideAttempt_Pastyear %in% c(0, 1),
      !is.na(Age_year), !is.na(Gender) 
    )
}

oneback_final <- process_data(oneback_data, SBQ_data_renamed)
twoback_final <- process_data(twoback_data, SBQ_data_renamed)
gng_final     <- process_data(GNG_data, SBQ_data_renamed)


# 3. 新增：名称美化函数 (核心改动点 1)
# 此函数将代码中的长变量名转换为图表中显示的短名称
get_pretty_name <- function(var_name) {
  # 创建一个查找列表 (key = 原始名, value = 漂亮名)
  name_map <- list(
    "Oneback_acc_deviationZ" = "1-Back Deviation-Z",
    "Twoback_acc_deviationZ" = "2-Back Deviation-Z",
    "d_prime_deviationZ"     = "Go/No-Go Deviation-Z",
    "SuicideAttempt_Lifetime" = "Lifetime Suicide Attempt",
    "SuicideAttempt_Pastyear" = "Past-Year Suicide Attempt"
  )
  
  # 如果在列表中找到匹配项，则返回漂亮名；否则返回原始名
  return(name_map[[var_name]] %||% var_name)
}


# 4. 稳健分层抽样函数 (与上一版相同，无需修改)
perform_robust_stratified_sampling <- function(group_0_data, group_1_data) {
  n_to_sample <- nrow(group_1_data)
  samples_to_take <- group_1_data %>% count(Gender, Age_year, name = "n_needed")
  pool_group_0 <- group_0_data
  final_sample_ids <- c()
  
  exact_matches <- pool_group_0 %>%
    inner_join(samples_to_take, by = c("Gender", "Age_year")) %>%
    group_by(Gender, Age_year) %>%
    sample_n(size = min(n(), first(n_needed))) %>%
    ungroup()
  
  if(nrow(exact_matches) > 0) {
    final_sample_ids <- c(final_sample_ids, exact_matches$x__ID)
    pool_group_0 <- pool_group_0 %>% filter(!x__ID %in% final_sample_ids)
  }
  
  n_matched_so_far <- exact_matches %>% count(Gender, Age_year, name = "n_matched")
  deficit_list <- samples_to_take %>%
    left_join(n_matched_so_far, by = c("Gender", "Age_year")) %>%
    mutate(n_matched = ifelse(is.na(n_matched), 0, n_matched), deficit = n_needed - n_matched) %>%
    filter(deficit > 0)
  
  if(nrow(deficit_list) > 0) {
    for(i in 1:nrow(deficit_list)) {
      target_row <- deficit_list[i, ]; target_gender <- target_row$Gender
      target_age <- target_row$Age_year; num_to_find <- target_row$deficit
      for(j in 1:num_to_find){
        if(nrow(pool_group_0) == 0) break
        best_fallback <- pool_group_0 %>%
          filter(Gender == target_gender) %>%
          mutate(age_diff = abs(Age_year - target_age)) %>%
          arrange(age_diff) %>% slice(1)
        if(nrow(best_fallback) > 0){
          final_sample_ids <- c(final_sample_ids, best_fallback$x__ID)
          pool_group_0 <- pool_group_0 %>% filter(!x__ID %in% final_sample_ids)
        }
      }
    }
  }
  return(group_0_data %>% filter(x__ID %in% final_sample_ids))
}


# 5. 核心分析函数 (无变化)
run_resampling_ttest <- function(data, test_var, group_var, n_iterations = 100) {
  group_0_data <- data %>% filter(.data[[group_var]] == 0); group_1_data <- data %>% filter(.data[[group_var]] == 1)
  n_small <- nrow(group_1_data)
  if (n_small < 3) {
    warning(paste("组别 '1' 的样本量过小 (n < 3) for", test_var, "by", group_var, ". 跳过分析。")); return(NULL)
  }
  successful_t_stats <- c(); successful_p_values <- c()
  cat("开始重采样，目标是获得 ", n_iterations, " 次成功的均衡样本...\n")
  max_attempts <- n_iterations * 5; attempt <- 1
  while(length(successful_t_stats) < n_iterations && attempt <= max_attempts) {
    group_0_sample <- perform_robust_stratified_sampling(group_0_data, group_1_data)
    if (nrow(group_0_sample) != n_small) { attempt <- attempt + 1; next }
    resampled_data <- bind_rows(group_0_sample, group_1_data)
    gender_test_p <- chisq.test(table(resampled_data$Gender, resampled_data[[group_var]]))$p.value
    age_test_p <- t.test(Age_year ~ get(group_var), data = resampled_data)$p.value
    if (gender_test_p > 0.05 && age_test_p > 0.05) {
      formula <- as.formula(paste(test_var, "~", group_var))
      t_test_result <- t.test(formula, data = resampled_data)
      successful_t_stats <- c(successful_t_stats, t_test_result$statistic)
      successful_p_values <- c(successful_p_values, t_test_result$p.value)
      cat("\r成功获取样本: ", length(successful_t_stats), "/", n_iterations)
    }
    attempt <- attempt + 1
  }
  cat("\n")
  if(length(successful_t_stats) < n_iterations) { warning("在最大尝试次数内未能收集到100个成功的均衡样本。仅返回已成功的部分。") }
  cat("重采样完成。\n"); return(data.frame(t_stat = successful_t_stats, p_value = successful_p_values))
}

# 6. 绘图函数 (核心改动点 2)
plot_resampling_results <- function(resampled_stats, test_var, group_var) {
  resampled_stats <- na.omit(resampled_stats)
  if (nrow(resampled_stats) == 0) {
    return(ggplot() + labs(title = "没有有效数据用于绘图", subtitle = paste(test_var, "by", group_var)))
  }
  
  # --- 在这里使用美化函数来创建标题 ---
  pretty_test_var <- get_pretty_name(test_var)
  pretty_group_var <- get_pretty_name(group_var)
  plot_title <- paste(pretty_test_var, "by", pretty_group_var)
  
  mean_t <- mean(resampled_stats$t_stat); median_p <- median(resampled_stats$p_value)
  prop_sig <- mean(resampled_stats$p_value < 0.05)
  stats_label <- paste0(
    "avg t-statistic = ", sprintf("%.2f", mean_t), "\n",
    "p-value median = ", sprintf("%.3f", median_p), "\n",
    "ratio(p<0.05) = ", percent(prop_sig, accuracy = 1)
  )
  x_range <- range(resampled_stats$t_stat, na.rm = TRUE)
  y_range_info <- ggplot_build(ggplot(resampled_stats, aes(x = t_stat)) + geom_histogram(bins=20))$layout$panel_scales_y[[1]]$range$range
  x_pos <- x_range[1] + 0.05 * (x_range[2] - x_range[1]); y_pos <- y_range_info[2] * 0.9
  
  hist_plot <- ggplot(resampled_stats, aes(x = t_stat)) +
    geom_histogram(bins = 20, fill = "lightblue", color = "black", alpha = 0.8) +
    geom_vline(xintercept = mean_t, color = "red", linetype = "dashed", size = 1) +
    annotate("text", x = x_pos, y = y_pos, label = stats_label, hjust = 0, vjust = 1, size = 3.5, lineheight = 1.2, fontface = "bold") +
    labs(
      title = plot_title, # 使用美化后的标题
      subtitle = paste0("t-statistic Distribution (from ", nrow(resampled_stats), " samples)"),
      x = "t-statistic",
      y = "Frequency"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 12, face = "bold"), # 也可以适当调整字体大小
      plot.subtitle = element_text(size = 9)
    )
  return(hist_plot)
}

# 7. 执行所有分析并生成图表 (无变化)
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
  params <- analyses_to_run[[i]]; cat("----------------------------------------------------------\n")
  cat("正在运行分析: ", params$test, " by ", params$group, "...\n")
  resampling_results <- run_resampling_ttest(params$data, params$test, params$group)
  if (!is.null(resampling_results) && nrow(resampling_results) > 0) {
    plot_name <- paste0("p", i)
    all_plots[[plot_name]] <- plot_resampling_results(resampling_results, params$test, params$group)
  }
}

# 8. 拼接并保存最终的PDF (无变化)
if (length(all_plots) > 0) {
  final_plot <- patchwork::wrap_plots(all_plots, ncol = 2) +
    plot_annotation(
      title = "100 bootstrap results",
      tag_levels = 'A' 
    ) & theme(plot.title = element_text(size = 16, face = "bold"))
  
  output_dir <- "D:/datasets/yunfu/results/SBQ"
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # 更新输出文件名以作区分
  output_pdf_path <- file.path(output_dir, "cognitive_performance_resampling_analysis_V7_pretty_titles.pdf")
  ggsave(output_pdf_path, final_plot, width = 14, height = 15, units = "in", device = "pdf")
  
  cat(paste("\n所有图表已成功拼接并保存到PDF文件:\n", output_pdf_path, "\n"))
  
} else {
  cat("\n没有生成任何图表，无法创建PDF。\n")
}