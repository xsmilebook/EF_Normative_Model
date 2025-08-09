library(ggplot2)
library(dplyr)
library(readxl)
library(gridExtra)

# --- 1. 数据读取 ---
switch <- read_excel("D:/datasets/yunfu/raw_data/switch_results/switch_dependent_var.xlsx")

# --- 2. 确保 Age_year 是数值型 ---
# (不进行四舍五入，保持其连续性)
# switch <- switch %>% mutate(Age_year = as.numeric(Age_year)) # 如果需要强制转换
# 检查 Age_year 的范围，用于固定 x 轴
age_range <- range(switch$Age_year, na.rm = TRUE)
cat("Age_year 的范围是:", age_range[1], "到", age_range[2], "\n")

# --- 3. 定义要绘制的指标 ---
# 注意：这里 y_column 是原始数据 switch 中的列名
plot_info <- tibble::tribble(
  ~y_column, ~title, ~y_label,
  # --- 原始指标 ---
  # sequence
  "M_acc", "M Block Accuracy", "Accuracy",
  "P_acc", "P Block Accuracy", "Accuracy",
  "B_acc", "Mix Trials Accuracy", "Accuracy",
  "M_trials_RT", "M Block RT", "Reaction Time (ms)",
  "P_trials_RT", "P Block RT", "Reaction Time (ms)",
  "B_rt", "Mix Trials RT", "Reaction Time (ms)",
  
  
  "Repeat_Trials_Acc", "Repeat Trials Accuracy", "Accuracy",
  "Switch_Trials_Acc", "Switch Trials Accuracy", "Accuracy",
  "Pure_trials_acc", "Pure Trials Accuracy", "Accuracy",
  "Repeat_Trials_RT", "Repeat Trials RT", "Reaction Time (ms)",
  "Switch_Trials_RT", "Switch Trials RT", "Reaction Time (ms)",
  "Pure_RT", "Pure Trials RT", "Reaction Time (ms)",
  
  "SC_SR_acc", "SC: Repeat - Switch Acc", "Accuracy Difference",
  "SC_SP_acc", "SC: Pure - Switch Acc", "Accuracy Difference",
  "SC_BP_acc", "SC: Pure - Mix Acc", "Accuracy Difference",
  "SC_SR_rt", "SC: Repeat - Switch RT", "RT Difference (ms)",
  "SC_SP_rt", "SC: Pure - Switch RT", "RT Difference (ms)",
  "SC_BP_rt", "SC: Pure - Mix RT", "RT Difference (ms)",
  
  # --- 计算出的 IES 指标 ---
  # 同样，假设这些列已存在
  "IES_all", "IES: Overall", "IES (RT/Acc)",
  "IES_pure", "IES: Pure Trials", "IES (RT/Acc)",
  "IES_mix", "IES: Mix Trials", "IES (RT/Acc)",
  "IES_re", "IES: Repeat Trials", "IES (RT/Acc)",
  "IES_switch", "IES: Switch Trials", "IES (RT/Acc)",
  
  # --- 计算出的 Proportion 指标 ---
  "SC_SR_rt_proportion", "SC RT Proportion (R-S)/R", "Proportion",
  "SC_SR_acc_proportion", "SC Acc Proportion (R-S)/R", "Proportion",
  "SC_SP_rt_proportion", "SC RT Proportion (P-S)/P", "Proportion",
  "SC_SP_acc_proportion", "SC Acc Proportion (P-S)/P", "Proportion",
  "Mean_RT", "Overall Mean RT", "Reaction Time (ms)"
)

# --- 在绘图前计算需要的复合指标 ---
# switch <- switch %>%
#   mutate(
#     # Switch Cost
#     SC_SR_acc = Repeat_Trials_Acc - Switch_Trials_Acc,
#     SC_SR_rt = Repeat_Trials_RT - Switch_Trials_RT,
#     SC_SP_acc = Pure_trials_acc - Switch_Trials_Acc,
#     SC_SP_rt = Pure_RT - Switch_Trials_RT,
#     SC_BP_acc = Pure_trials_acc - B_acc,
#     SC_BP_rt = Pure_RT - B_rt,
#     # IES (注意处理 Acc 为 0 的情况)
#     IES_all = Mean_RT / Switch_acc,
#     IES_pure = Pure_RT / Pure_trials_acc,
#     IES_mix = B_rt / B_acc,
#     IES_re = Repeat_Trials_RT / Repeat_Trials_Acc,
#     IES_switch = Switch_Trials_RT / Switch_Trials_Acc,
#     # Proportion
#     SC_SR_rt_proportion = (Repeat_Trials_RT - Switch_Trials_RT) / Repeat_Trials_RT,
#     SC_SR_acc_proportion = (Repeat_Trials_Acc - Switch_Trials_Acc) / Repeat_Trials_Acc,
#     SC_SP_rt_proportion = (Pure_RT - Switch_Trials_RT) / Pure_RT,
#     SC_SP_acc_proportion = (Pure_trials_acc - Switch_Trials_Acc) / Pure_trials_acc
#   )

# --- 4. 生成图表 ---
plot_list <- list()

png_output_dir <- "D:/datasets/yunfu/raw_data/switch_results/figures"

# for (i in 1:2) {
for (i in 1:nrow(plot_info)) {
  y_col <- plot_info$y_column[i]
  p_title <- plot_info$title[i]
  y_lab <- plot_info$y_label[i]
  
  # 检查列是否存在
  if (!y_col %in% names(switch)) {
    warning(paste("列", y_col, "在 switch 数据中不存在，跳过绘图。"))
    plot_list[[i]] <- ggplot() + annotate("text", x=0, y=0, label=paste("Missing:", y_col)) + theme_void()
    next
  }
  
  # 过滤掉 y 值或 x 值为非有限的行，防止绘图错误
  # 同时移除 IES 中可能的 Inf 值
  plot_data <- switch[is.finite(switch[[y_col]]) & is.finite(switch$Age_year) & 
                        (!startsWith(y_col, "IES") | abs(switch[[y_col]]) != Inf), ]
  
  if (nrow(plot_data) == 0) {
    warning(paste("列", y_col, "没有有效的有限数值用于绘图，跳过。"))
    plot_list[[i]] <- ggplot() + annotate("text", x=0, y=0, label=paste("No valid ", y_col)) + theme_void()
    next
  }
  
  # 构建 ggplot 对象
  # 使用所有原始数据点进行拟合
  if (grepl("RT", y_col, ignore.case = TRUE)) {
    y_scale <- scale_y_continuous(expand = c(0.05, 0))
  } else if (grepl("Acc", y_col, ignore.case = TRUE) & !grepl("SC", y_col, ignore.case = TRUE)) {
    # y_scale <- scale_y_continuous(limits = c(0.5, 0.95))
    y_scale <- scale_y_continuous(expand = c(0.05, 0))
  } else {
    # 对于其他指标 (如 SC, IES, Proportion) 使用默认扩展
    y_scale <- scale_y_continuous(expand = c(0.05, 0))
  }
  
  p <- ggplot(plot_data, aes(x = Age_year, y = !!sym(y_col))) +
    # 使用轻微透明度的点，因为点会很多
    # geom_point(color = "#FFB6C1", size = 1, alpha = 0.3) + 
    # 使用 GAM 进行平滑拟合，置信区间基于模型
    geom_smooth(method = "gam", formula = y ~ s(x, k = 3), 
                se = TRUE, 
                color = "black", 
                fill = "grey",   # 改变置信区间的颜色
                alpha = 0.3) +   # 调整置信区间的透明度
    theme_minimal() +
    theme(
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 13, hjust = 0.5),
      panel.grid.minor = element_blank()
    ) +
    # 根据数据的实际范围设置 x 轴
    # 在循环中根据 y_col 设置不同的范围
    y_scale + 
    scale_x_continuous(limits = age_range, expand = c(0.01, 0.01)) +
    labs(title = p_title, x = "Age (years)", y = y_lab)
  
  plot_list[[i]] <- p
  # --- 新增：保存当前图表为 PNG ---

  # gsub("[^a-zA-Z0-9_]", "_", ...) 将所有非字母、数字、下划线的字符替换为下划线
  safe_title <- gsub("[^a-zA-Z0-9_]", "_", p_title) 
  # 定义文件名，例如 "plot_01_SC_Pure__Mix_RT.png"
  filename <- file.path(png_output_dir, paste0("plot_", sprintf("%02d", i), "_", safe_title, ".jpg"))
  
  # 使用 ggsave 保存图表
  # width, height: 图片尺寸 (英寸)
  # dpi: 图片分辨率 (dots per inch)，300 是常用的打印质量
  # 可以根据需要调整这些参数
  ggsave(filename, plot = p, width = 8, height = 6, dpi = 300) 
  
  cat("已保存图表:", filename, "\n") # 打印保存信息
}

# --- 5. 保存到 PDF ---
output_pdf_path <- "D:/datasets/yunfu/raw_data/switch_results/all_switch_analysis_plots_continuous_age.pdf"

# 计算需要多少页 (每页16图)
num_plots <- length(plot_list)
plots_per_page <- 12 # 4x4
num_pages <- ceiling(num_plots / plots_per_page)

if (num_plots == 0) {
  stop("没有生成任何图表，无法保存 PDF。")
}

# 打开PDF设备
pdf(output_pdf_path, width = 16, height = 12)

# 分页绘制
for (page in 1:num_pages) {
  start_idx <- (page - 1) * plots_per_page + 1
  end_idx <- min(page * plots_per_page, num_plots)
  current_plots <- plot_list[start_idx:end_idx]
  
  grid.arrange(grobs = current_plots, nrow = 4, ncol = 3, top = paste("Page", page))
}

# 关闭PDF设备
dev.off()

cat("成功生成", num_plots, "个图表并保存至:", output_pdf_path, "\n")
cat("共", num_pages, "页。\n")

