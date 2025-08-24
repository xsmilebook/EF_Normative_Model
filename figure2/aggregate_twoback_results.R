# 1. 加载必要的包
library(tidyverse)
library(ggplot2)
library(patchwork)

# ==============================================================================
# 2. 设置路径和关键变量
# +++ 唯一的修改点: 将任务名称改为 "twoback" +++
task_name <- "twoback" 

wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  resultFolder_individual <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/datasets/yunfu/results/EF_psy"
  resultFolder_final <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/EF_Normative_Model/anova/results"
  FigureFolder_final <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/EF_Normative_Model/anova/figures"
} else {
  resultFolder_individual <- "D:/datasets/yunfu/results/EF_psy"
  resultFolder_final <- "D:/datasets/yunfu/results/EF_psy_anovaPB"
  FigureFolder_final <- 'D:/datasets/yunfu/figures/fig2/anovaPB'
}

# 确保输出文件夹存在
dir.create(resultFolder_final, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureFolder_final, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 3. 识别所有独特的心理学变量 (psyvar)
# (无需修改，代码会自动适应 task_name 的变化)
rds_pattern <- paste0("^anova_simulation_", task_name, "_.*\\.rds$")
all_rds_files <- list.files(path = resultFolder_individual, pattern = rds_pattern, full.names = FALSE)

if (length(all_rds_files) == 0) {
  stop(paste("错误: 在指定文件夹中找不到任何与", task_name, "相关的 'anova_simulation' RDS 文件。"))
}

regex_pattern <- paste0("anova_simulation_", task_name, "_(.*)_Time_\\d+\\.rds$")
unique_psy_vars <- unique(str_match(all_rds_files, regex_pattern)[,2])

cat(paste("找到以下与", task_name, "任务相关的独特心理学变量进行合并:\n"))
print(unique_psy_vars)

# ==============================================================================
# 4. 循环合并每个变量的模拟结果
final_results_list <- list()

for (psyvar in unique_psy_vars) {
  
  cat(paste("\n===== 正在合并变量:", psyvar, "=====\n"))
  
  pattern_rds <- paste0("anova_simulation_", task_name, "_", psyvar, "_Time_.*\\.rds$")
  rds_files_for_var <- list.files(path = resultFolder_individual, pattern = pattern_rds, full.names = TRUE)
  
  if (length(rds_files_for_var) == 0) {
    cat(paste("警告: 找不到变量", psyvar, "的RDS文件，跳过。\n"))
    next
  }
  
  all_sim_data <- map(rds_files_for_var, readRDS)
  combined_sim_stats <- unlist(map(all_sim_data, ~ .$simulation$simulated_stats))
  observed_stat <- all_sim_data[[1]]$simulation$observed_stat
  final_p_value <- mean(combined_sim_stats >= observed_stat, na.rm = TRUE)
  
  cat(paste("合并了", length(combined_sim_stats), "次模拟。\n"))
  cat(paste("观测值:", round(observed_stat, 4), "\n"))
  cat(paste("重新计算的 P-value:", format.pval(final_p_value, digits = 4, eps = 0.0001), "\n"))
  
  pattern_csv <- paste0("corr_", task_name, "_", psyvar, "_Time_.*\\.csv$")
  csv_file_template <- list.files(path = resultFolder_individual, pattern = pattern_csv, full.names = TRUE)[1]
  
  if (is.na(csv_file_template)) {
    cat(paste("警告: 找不到变量", psyvar, "的CSV文件模板，无法生成统计摘要。\n"))
    next
  }
  
  final_summary_row <- read.csv(csv_file_template)
  final_summary_row$anova.pvalues <- final_p_value
  final_summary_row$Time_id <- "merged_10000"
  final_summary_row$n_sim <- length(combined_sim_stats)
  final_results_list[[psyvar]] <- final_summary_row
  
  bootstrap_folder_final <- file.path(FigureFolder_final, "bootstrap_distributions_final")
  dir.create(bootstrap_folder_final, showWarnings = FALSE, recursive = TRUE)
  
  clean_psyvar <- gsub("[/:*?\"<>|]", "_", psyvar)
  file_path <- file.path(bootstrap_folder_final, paste0("sim_dist_final_", task_name, "_", clean_psyvar, ".png"))
  
  png(filename = file_path, width = 800, height = 600, res = 150)
  hist_title <- paste("Combined Bootstrap Distribution (N=", length(combined_sim_stats), ")\n", task_name, "vs", psyvar)
  hist(combined_sim_stats, 
       main = hist_title,
       xlab = "Deviance Difference", col = "#56B1F7", border = "white", breaks = 100)
  abline(v = observed_stat, col = "#D55E00", lwd = 2.5)
  legend("topright", 
         legend = paste("Observed =", round(observed_stat, 3), 
                        "\nP-value =", format.pval(final_p_value, digits = 4, eps = 0.0001)), 
         bty = "n", cex = 1.2)
  dev.off()
  
  cat(paste("最终的分布图已保存至:", file_path, "\n"))
}

# ==============================================================================
# 5. 将所有变量的最终结果合并成一个数据框
corr.result.df <- bind_rows(final_results_list)

# ==============================================================================
# 6. 进行多重比较校正 (Bonferroni)
corr.result.df <- corr.result.df %>%
  mutate(
    correstimate = as.numeric(correstimate),
    anova.pvalues = as.numeric(anova.pvalues),
    corrp = as.numeric(corrp),
    anovap.bonf = p.adjust(anova.pvalues, method = "bonferroni"),
    sig = anovap.bonf < 0.05
  )

# ==============================================================================
# 7. 绘制最终的热图
# 提醒: 如果 twoback 任务分析的心理学变量不同，您需要修改下面的 y_levels 和 y_labels
y_levels <- c("SDQ_PP_sum_z", "SDQ_ES_sum_z") 
y_labels <- c("SDQ_ES_sum_z" = "Emotional Symptoms", "SDQ_PP_sum_z" = "Peer Problems")

corr.result.df$parcel <- factor(corr.result.df$psyvar, levels = rev(y_levels))

plot_title <- paste("Correlation between", task_name, "and Psychiatric Scores")

Fig <- ggplot(corr.result.df, aes(x = period, y = parcel, fill = correstimate)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(data = filter(corr.result.df, sig == TRUE), 
            aes(label = "*"), vjust = 0.75, size = 10, color = "black") +
  scale_fill_distiller(palette = "RdBu", direction = -1, name = "Estimate") +
  scale_y_discrete(labels = y_labels) +
  labs(
    title = plot_title,
    x = NULL, y = "Psychiatric Scores"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    panel.grid = element_blank()
  )

# ==============================================================================
# 8. 保存最终结果
final_plot_path <- file.path(FigureFolder_final, paste0(task_name, "_correlation_heatmap_final.pdf"))
ggsave(final_plot_path, plot = Fig, width = 8, height = 7)
cat("\n最终热图已保存至:", final_plot_path, "\n")

final_csv_path <- file.path(resultFolder_final, paste0(task_name, "_corr_results_with_bonf_aggregated.csv"))
write.csv(corr.result.df, file = final_csv_path, row.names = FALSE)
cat("最终Bonferroni校正结果已保存至:", final_csv_path, "\n")

cat("\n所有合并和分析已完成！\n")