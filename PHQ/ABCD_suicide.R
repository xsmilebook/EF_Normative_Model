# 1. 加载所需的库
# 确保已安装所有必要的包: install.packages(c("dplyr", "openxlsx", "ggplot2", "patchwork", "car"))
library(dplyr)
library(openxlsx)
library(ggplot2)
library(patchwork) # 用于拼接多张图表
library(car)       # 用于执行Levene's test

# 2. 读取、准备并合并您的新数据
datapath  <- "D:/datasets/yunfu/interfile_folder/ABCD"
flanker_data <- read.csv(file.path(datapath, "Flanker.deviations.csv"))
# ksads_data <- read.xlsx(file.path(datapath, "mh_p_ksads_SuicideAttempt.xlsx"))
ksads_data <- read.csv(file.path(datapath, "new", "mh_y_ksads_ss.csv"))
ksads_prepared <- ksads_data %>%
  mutate(ksads_import_id_t = paste(src_subject_id, eventname, sep = "_"))

flanker_prepared <- flanker_data %>%
  rename(ksads_import_id_t = ID)

merged_data <- inner_join(
  ksads_prepared,
  flanker_prepared,
  by = "ksads_import_id_t"
)

# 3. 数据过滤与最终准备 (新增分组变量)
analysis_ready_data <- merged_data %>%
  mutate(
    ksads_23_954_t = as.numeric(as.character(ksads_23_954_t)),
    ksads_23_965_t = as.numeric(as.character(ksads_23_965_t))
  ) %>%
  # 过滤掉两个变量都为NA的行，以及认知数据为NA的行
  filter(
    (!is.na(ksads_23_954_t) | !is.na(ksads_23_965_t)) &
      !is.na(nihtbx_flanker_uncorrected_deviationZ)
  ) %>%
  # --- 新增步骤：创建新的分组变量 'any_suicide_item' ---
  mutate(
    any_suicide_item = case_when(
      # 条件1: 任意一个为1，则新变量为1 (风险组)
      ksads_23_954_t == 1 | ksads_23_965_t == 1 ~ 1,
      # 条件2: 两个都为0，则新变量为0 (对照组)
      ksads_23_954_t == 0 & ksads_23_965_t == 0 ~ 0,
      # 其他情况 (例如一个为0一个为NA) 标记为NA，后续会被自动过滤
      TRUE ~ NA_real_
    )
  ) %>%
  # 确保新变量存在有效值 (0 或 1)
  filter(!is.na(any_suicide_item))


# 4. 复用并修正分析与绘图函数 (无需改动，可直接复用)
analyze_and_plot <- function(data, test_var, group_var) {
  
  plot_data <- data %>%
    filter(!is.na(.data[[group_var]]), .data[[group_var]] %in% c(0, 1)) %>%
    mutate(!!sym(group_var) := as.factor(.data[[group_var]]))
  
  plot_title <- paste(test_var, "\nby", group_var)
  
  sample_sizes <- plot_data %>%
    group_by(.data[[group_var]]) %>%
    summarise(n = n(), .groups = 'drop')
  
  if (nrow(sample_sizes) < 2 || any(sample_sizes$n < 3)) {
    cat(paste0("Skipping plot '", plot_title, "' due to insufficient data for comparison.\n"))
    return(NULL)
  }
  
  axis_labels <- paste0(sample_sizes[[group_var]], "\n(n=", sample_sizes$n, ")")
  
  stats_label <- tryCatch({
    formula <- as.formula(paste(test_var, "~", group_var))
    
    levene_res <- leveneTest(formula, data = plot_data)
    levene_p_val <- levene_res$`Pr(>F)`[1]
    levene_f_val <- levene_res$`F value`[1]
    levene_label <- paste0("Levene's F = ", sprintf("%.2f", levene_f_val), 
                           ", p = ", sprintf("%.3f", levene_p_val))
    
    t_test_result <- t.test(formula, data = plot_data, var.equal = FALSE)
    
    t_val_formatted <- sprintf("%.2f", t_test_result$statistic)
    p_val_formatted <- if (t_test_result$p.value < 0.001) {
      "p < 0.001"
    } else {
      paste("p =", sprintf("%.3f", t_test_result$p.value))
    }
    
    paste0(levene_label, "\nWelch's t = ", t_val_formatted, ", ", p_val_formatted)
    
  }, error = function(e) { "" })
  
  y_range <- range(plot_data[[test_var]], na.rm = TRUE)
  y_pos <- y_range[1] + 0.95 * (y_range[2] - y_range[1])
  
  p <- ggplot(plot_data, aes_string(x = group_var, y = test_var)) +
    geom_jitter(width = 0.2, height = 0, color = "lightblue", alpha = 0.3) +
    geom_boxplot(width = 0.5, alpha = 0.6, outlier.shape = NA, aes_string(fill = group_var)) +
    annotate("text", x = 1.5, y = y_pos, label = stats_label, size = 3.5, lineheight = 1.2) +
    scale_x_discrete(labels = axis_labels) +
    labs(title = plot_title, x = group_var, y = test_var) +
    theme_bw(base_size = 12) + 
    theme(legend.position = "none",
          plot.title = element_text(size = 10))
  
  return(p)
}


# 5. 生成所有图表并存储在一个列表中 (新增第三个图)
plot_list <- list()

plot_list$p1 <- analyze_and_plot(analysis_ready_data, "nihtbx_flanker_uncorrected_deviationZ", "ksads_23_954_t")
plot_list$p2 <- analyze_and_plot(analysis_ready_data, "nihtbx_flanker_uncorrected_deviationZ", "ksads_23_965_t")
# 新增的分析: 使用新创建的 'any_suicide_item' 变量
plot_list$p3 <- analyze_and_plot(analysis_ready_data, "nihtbx_flanker_uncorrected_deviationZ", "any_suicide_item")


# 6. 移除列表中可能存在的NULL元素
if (requireNamespace("purrr", quietly = TRUE)) {
  plot_list <- purrr::compact(plot_list)
} else {
  plot_list <- plot_list[!sapply(plot_list, is.null)]
}


# 7. 使用 patchwork 拼接图表
if (length(plot_list) > 0) {
  # 将3个图拼接成 2x2 的网格，最后一个位置留空
  final_plot <- patchwork::wrap_plots(plot_list, ncol = 2) +
    plot_annotation(
      title = "ABCD Flanker Task Performance by KSADS Suicide Items",
      tag_levels = 'A' 
    )
  
  # 8. 将拼接好的图表保存为 PDF 文件
  output_dir <- "D:/datasets/yunfu/results/ABCD"
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  output_pdf_path <- file.path(output_dir, "flanker_ksads_analysis_combined.pdf")
  # 调整高度以适应两行图表
  ggsave(output_pdf_path, final_plot, width = 12, height = 10, units = "in", device = "pdf")
  
  cat(paste("\nABCD数据集分析图表（含合并组比较）已成功保存到PDF文件:\n", output_pdf_path, "\n"))
  
} else {
  cat("\n没有生成任何图表，无法创建PDF。\n")
}