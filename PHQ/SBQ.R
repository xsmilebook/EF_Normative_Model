# 1. 加载所需的库
# 确保已安装所有必要的包: install.packages(c("dplyr", "openxlsx", "ggplot2", "patchwork"))
library(dplyr)
library(openxlsx)
library(ggplot2)
library(patchwork) # 用于拼接多张图表
# 注意：我们不再需要 ggpubr 包

# 2. 定义文件路径并加载数据 (请确保文件路径正确)
datapath <- "D:/datasets/yunfu/interfile_folder/deviationZ_temp"
sbq_path <- "D:/datasets/yunfu/raw_data/2.0_云浮队列基线学生问卷_自杀相关变量.xlsx"

# 检查文件是否存在，避免出错
if (!dir.exists(datapath) || !file.exists(sbq_path)) {
  stop("一个或多个文件路径不正确，请检查datapath和sbq_path变量。")
}

oneback_data <- readRDS(file.path(datapath, "back1Acc.deviations.rds"))
twoback_data <- readRDS(file.path(datapath, "back2Acc.deviations.rds"))
GNG_data <- readRDS(file.path(datapath, "GNGd_prime.deviations.rds"))
SBQ_data <- read_xlsx(sbq_path)

# 3. 数据预处理
SBQ_data$"用户ID" <- as.numeric(SBQ_data$"用户ID")
SBQ_data_renamed <- SBQ_data %>%
  rename(x__ID = "用户ID")

# 4. 定义一个统一的数据处理函数
process_data <- function(task_data, sbq_data) {
  task_data %>%
    inner_join(sbq_data, by = "x__ID") %>%
    filter(
      SuicideAttempt_Lifetime %in% c(0, 1),
      SuicideAttempt_Pastyear %in% c(0, 1)
    ) %>%
    mutate(
      SuicideAttempt_Lifetime = as.factor(SuicideAttempt_Lifetime),
      SuicideAttempt_Pastyear = as.factor(SuicideAttempt_Pastyear)
    )
}

# 分别处理三个任务的数据
oneback_final <- process_data(oneback_data, SBQ_data_renamed)
twoback_final <- process_data(twoback_data, SBQ_data_renamed)
gng_final     <- process_data(GNG_data, SBQ_data_renamed)


# 5. 创建一个用于分析和绘图的函数 (最终稳健版本)
analyze_and_plot <- function(data, test_var, group_var) {
  
  # 创建图表标题
  plot_title <- paste(test_var, "by", group_var)
  
  # 计算样本量并创建X轴标签
  sample_sizes <- data %>%
    group_by(.data[[group_var]]) %>%
    summarise(n = n(), .groups = 'drop')
  
  # 检查是否有足够的组别和数据进行比较
  if (nrow(sample_sizes) < 2 || any(sample_sizes$n < 2)) {
    cat(paste0("Skipping plot '", plot_title, "' due to insufficient data in one or both groups.\n"))
    return(NULL) # 如果不足两组或任一组样本小于2，则跳过
  }
  
  axis_labels <- paste0(sample_sizes[[group_var]], "\n(n=", sample_sizes$n, ")")
  
  # --- 新增: 手动执行t检验并创建标签 ---
  stats_label <- tryCatch({
    formula <- as.formula(paste(test_var, "~", group_var))
    t_test_result <- t.test(formula, data = data)
    
    t_val_formatted <- sprintf("%.2f", t_test_result$statistic)
    # 格式化p值，使其更易读
    p_val_formatted <- if (t_test_result$p.value < 0.001) {
      "p < 0.001"
    } else {
      paste("p =", sprintf("%.3f", t_test_result$p.value))
    }
    
    paste0("t = ", t_val_formatted, "\n", p_val_formatted)
  }, error = function(e) {
    # 如果t检验因任何原因失败，返回一个空标签
    ""
  })
  
  # --- 确定标签的Y轴位置 ---
  y_range <- range(data[[test_var]], na.rm = TRUE)
  # 将标签放置在Y轴范围顶部90%的位置
  y_pos <- y_range[1] + 0.95 * (y_range[2] - y_range[1])
  
  # 使用 ggplot2 绘制箱线散点图
  p <- ggplot(data, aes_string(x = group_var, y = test_var)) +
    geom_jitter(width = 0.2, height = 0, color = "lightblue", alpha = 0.3) +
    geom_boxplot(width = 0.5, alpha = 0.6, outlier.shape = NA, aes_string(fill = group_var)) +

    # *** 使用 annotate 添加手动创建的标签 ***
    annotate(
      geom = "text",
      x = 1.5, # 将标签放置在两组 (组1和组2) 的中间
      y = y_pos,
      label = stats_label,
      size = 4,
      lineheight = 1 # 调整行距
    ) +
    
    # 应用带有样本量的新标签
    scale_x_discrete(labels = axis_labels) +
    
    labs(
      title = plot_title,
      x = group_var,
      y = test_var
    ) +
    theme_bw(base_size = 12) + 
    theme(legend.position = "none")
  
  return(p)
}

# 6. 生成所有图表并存储在一个列表中
plot_list <- list()

# --- One-Back 分析 ---
plot_list$p1 <- analyze_and_plot(oneback_final, "Oneback_acc_deviationZ", "SuicideAttempt_Lifetime")
plot_list$p2 <- analyze_and_plot(oneback_final, "Oneback_acc_deviationZ", "SuicideAttempt_Pastyear")

# --- Two-Back 分析 ---
plot_list$p3 <- analyze_and_plot(twoback_final, "Twoback_acc_deviationZ", "SuicideAttempt_Lifetime")
plot_list$p4 <- analyze_and_plot(twoback_final, "Twoback_acc_deviationZ", "SuicideAttempt_Pastyear")

# --- Go/No-Go 分析 ---
plot_list$p5 <- analyze_and_plot(gng_final, "d_prime_deviationZ", "SuicideAttempt_Lifetime")
plot_list$p6 <- analyze_and_plot(gng_final, "d_prime_deviationZ", "SuicideAttempt_Pastyear")

# 7. 移除列表中可能存在的NULL元素（由跳过的绘图产生）
plot_list <- purrr::compact(plot_list)

# 8. 使用 patchwork 拼接图表
if (length(plot_list) > 0) {
  # 根据图表数量自动调整PDF高度
  pdf_height <- 5 * ceiling(length(plot_list) / 2)
  final_plot <- patchwork::wrap_plots(plot_list, ncol = 2) +
    plot_annotation(
      title = "Cognitive Task Performance by Suicide Attempt History",
      tag_levels = 'A' 
    )
  
  # 9. 将拼接好的图表保存为 PDF 文件
  output_pdf_path <- "D:/datasets/yunfu/results/SBQ/cognitive_performance_analysis_with_N.pdf"
  ggsave(output_pdf_path, final_plot, width = 12, height = pdf_height, units = "in")
  
  # 打印消息，告知用户文件已保存
  cat(paste("\n所有图表已成功拼接并保存到PDF文件:\n", file.path(getwd(), output_pdf_path), "\n"))
  
} else {
  cat("\n没有生成任何图表，无法创建PDF。\n")
}