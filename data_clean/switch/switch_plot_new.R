library(ggplot2)
library(dplyr)
library(readxl)
library(gridExtra)

# --- 1. 数据读取 ---
switch <- read_excel("D:/datasets/yunfu/raw_data/switch_results/switch_dependent_var_new.xlsx")

# --- 2. 数据处理：四舍五入年龄 ---
switch <- switch %>%
  mutate(Age_year = round(Age_year, 0))

# --- 3. 数据汇总 ---
# 包含所有需要绘制的原始指标和计算出的指标
switch_ana <- switch %>%
  group_by(Age_year) %>%
  summarise(
    # --- 原始指标 (Mean, SD, SE) ---
    Mean_RT = mean(Mean_RT, na.rm = TRUE),
    sd_Mean_RT = sd(Mean_RT, na.rm = TRUE),
    se_Mean_RT = sd_Mean_RT / sqrt(n()),
    Switch_acc = mean(Switch_acc, na.rm = TRUE),
    sd_Switch_acc = sd(Switch_acc, na.rm = TRUE),
    se_Switch_acc = sd_Switch_acc / sqrt(n()),
    Repeat_Trials_Acc = mean(Repeat_Trials_Acc, na.rm = TRUE),
    sd_Repeat_Trials_Acc = sd(Repeat_Trials_Acc, na.rm = TRUE),
    se_Repeat_Trials_Acc = sd_Repeat_Trials_Acc / sqrt(n()),
    Repeat_Trials_RT = mean(Repeat_Trials_RT, na.rm = TRUE),
    sd_Repeat_Trials_RT = sd(Repeat_Trials_RT, na.rm = TRUE),
    se_Repeat_Trials_RT = sd_Repeat_Trials_RT / sqrt(n()),
    Switch_Trials_Acc = mean(Switch_Trials_Acc, na.rm = TRUE),
    sd_Switch_Trials_Acc = sd(Switch_Trials_Acc, na.rm = TRUE),
    se_Switch_Trials_Acc = sd_Switch_Trials_Acc / sqrt(n()),
    Switch_Trials_RT = mean(Switch_Trials_RT, na.rm = TRUE),
    sd_Switch_Trials_RT = sd(Switch_Trials_RT, na.rm = TRUE),
    se_Switch_Trials_RT = sd_Switch_Trials_RT / sqrt(n()),
    Pure_trials_acc = mean(Pure_trials_acc, na.rm = TRUE),
    sd_Pure_trials_acc = sd(Pure_trials_acc, na.rm = TRUE),
    se_Pure_trials_acc = sd_Pure_trials_acc / sqrt(n()),
    Pure_RT = mean(Pure_RT, na.rm = TRUE),
    sd_Pure_RT = sd(Pure_RT, na.rm = TRUE),
    se_Pure_RT = sd_Pure_RT / sqrt(n()),
    B_acc = mean(B_acc, na.rm = TRUE),
    sd_B_acc = sd(B_acc, na.rm = TRUE),
    se_B_acc = sd_B_acc / sqrt(n()),
    B_rt = mean(B_rt, na.rm = TRUE),
    sd_B_rt = sd(B_rt, na.rm = TRUE),
    se_B_rt = sd_B_rt / sqrt(n()),
    M_acc = mean(M_acc, na.rm = TRUE),
    sd_M_acc = sd(M_acc, na.rm = TRUE),
    se_M_acc = sd_M_acc / sqrt(n()),
    M_trials_RT = mean(M_trials_RT, na.rm = TRUE),
    sd_M_trials_RT = sd(M_trials_RT, na.rm = TRUE),
    se_M_trials_RT = sd_M_trials_RT / sqrt(n()),
    P_acc = mean(P_acc, na.rm = TRUE),
    sd_P_acc = sd(P_acc, na.rm = TRUE),
    se_P_acc = sd_P_acc / sqrt(n()),
    P_trials_RT = mean(P_trials_RT, na.rm = TRUE),
    sd_P_trials_RT = sd(P_trials_RT, na.rm = TRUE),
    se_P_trials_RT = sd_P_trials_RT / sqrt(n()),
    # --- 计算出的 Switch Cost 指标 (Mean, SD, SE) ---
    SC_SR_acc = mean(Repeat_Trials_Acc - Switch_Trials_Acc, na.rm = TRUE),
    sd_SC_SR_acc = sd(Repeat_Trials_Acc - Switch_Trials_Acc, na.rm = TRUE),
    se_SC_SR_acc = sd_SC_SR_acc / sqrt(n()),
    SC_SR_rt = mean(Repeat_Trials_RT - Switch_Trials_RT, na.rm = TRUE),
    sd_SC_SR_rt = sd(Repeat_Trials_RT - Switch_Trials_RT, na.rm = TRUE),
    se_SC_SR_rt = sd_SC_SR_rt / sqrt(n()),
    SC_SP_acc = mean(Pure_trials_acc - Switch_Trials_Acc, na.rm = TRUE),
    sd_SC_SP_acc = sd(Pure_trials_acc - Switch_Trials_Acc, na.rm = TRUE),
    se_SC_SP_acc = sd_SC_SP_acc / sqrt(n()),
    SC_SP_rt = mean(Pure_RT - Switch_Trials_RT, na.rm = TRUE),
    sd_SC_SP_rt = sd(Pure_RT - Switch_Trials_RT, na.rm = TRUE),
    se_SC_SP_rt = sd_SC_SP_rt / sqrt(n()),
    SC_BP_acc = mean(Pure_trials_acc - B_acc, na.rm = TRUE),
    sd_SC_BP_acc = sd(Pure_trials_acc - B_acc, na.rm = TRUE),
    se_SC_BP_acc = sd_SC_BP_acc / sqrt(n()),
    SC_BP_rt = mean(Pure_RT - B_rt, na.rm = TRUE),
    sd_SC_BP_rt = sd(Pure_RT - B_rt, na.rm = TRUE),
    se_SC_BP_rt = sd_SC_BP_rt / sqrt(n()),
    # --- 计算出的 IES 指标 (Mean) ---
    IES_all = mean(Mean_RT / Switch_acc, na.rm = TRUE),
    IES_pure = mean(Pure_RT / Pure_trials_acc, na.rm = TRUE),
    IES_mix = mean(B_rt / B_acc, na.rm = TRUE),
    IES_re = mean(Repeat_Trials_RT / Repeat_Trials_Acc, na.rm = TRUE),
    IES_switch = mean(Switch_Trials_RT / Switch_Trials_Acc, na.rm = TRUE),
    # --- 计算出的 Proportion 指标 (Mean) ---
    SC_SR_rt_proportion = mean((Repeat_Trials_RT - Switch_Trials_RT) / Repeat_Trials_RT, na.rm = TRUE),
    SC_SR_acc_proportion = mean((Repeat_Trials_Acc - Switch_Trials_Acc) / Repeat_Trials_Acc, na.rm = TRUE),
    SC_SP_rt_proportion = mean((Pure_RT - Switch_Trials_RT) / Pure_RT, na.rm = TRUE),
    SC_SP_acc_proportion = mean((Pure_trials_acc - Switch_Trials_Acc) / Pure_trials_acc, na.rm = TRUE),
    .groups = 'drop'
  )

# --- 4. 定义要绘制的指标 ---
plot_info <- tibble::tribble(
  ~y_column, ~title, ~y_label,
  # --- 原始指标 ---
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
  
  "SC_SR_acc", "SC: Switch - Repeat Acc", "Accuracy Difference",
  "SC_SP_acc", "SC: Switch - Pure Acc", "Accuracy Difference",
  "SC_BP_acc", "SC: Mixed - Pure Acc", "Accuracy Difference",
  "SC_SR_rt", "SC: Switch - Repeat RT", "RT Difference (ms)",
  "SC_SP_rt", "SC: Switch - Pure RT", "RT Difference (ms)",
  "SC_BP_rt", "SC: Mix - Pure RT", "RT Difference (ms)",
  # --- IES 指标 ---
  "IES_all", "IES: Overall", "IES (RT/Acc)",
  "IES_pure", "IES: Pure Trials", "IES (RT/Acc)",
  "IES_mix", "IES: Mix Trials", "IES (RT/Acc)",
  "IES_re", "IES: Repeat Trials", "IES (RT/Acc)",
  "IES_switch", "IES: Switch Trials", "IES (RT/Acc)",
  # --- Proportion 指标 ---
  "SC_SR_rt_proportion", "SC RT Proportion (R-S)/R", "Proportion",
  "SC_SR_acc_proportion", "SC Acc Proportion (R-S)/R", "Proportion",
  "SC_SP_rt_proportion", "SC RT Proportion (P-S)/P", "Proportion",
  "SC_SP_acc_proportion", "SC Acc Proportion (P-S)/P", "Proportion"
)

# --- 5. 生成图表 ---
plot_list <- list()

# --- 新增：定义保存单个 PNG 图片的目录 ---
png_output_dir <- "D:/datasets/yunfu/raw_data/switch_new_results/figures_no_points"
if (!dir.exists(png_output_dir)) {
  dir.create(png_output_dir, recursive = TRUE)
}

for (i in 1:nrow(plot_info)) {
  y_col <- plot_info$y_column[i]
  p_title <- plot_info$title[i]
  y_lab <- plot_info$y_label[i]
  
  # 检查列是否存在
  if (!y_col %in% names(switch_ana)) {
    warning(paste("列", y_col, "在 switch_ana 中不存在，跳过绘图。"))
    plot_list[[i]] <- ggplot() + annotate("text", x=0, y=0, label=paste("Missing:", y_col)) + theme_void()
    
    # 保存占位图
    safe_title <- gsub("[^a-zA-Z0-9_]", "_", p_title)
    filename <- file.path(png_output_dir, paste0("plot_", sprintf("%02d", i), "_", safe_title, ".png"))
    ggsave(filename, plot = plot_list[[i]], width = 8, height = 6, dpi = 300)
    cat("已保存占位图:", filename, "\n")
    next
  }
  
  # --- 新增：数据清理 ---
  plot_data <- switch_ana[is.finite(switch_ana[[y_col]]) & is.finite(switch_ana$Age_year), ]
  
  if (nrow(plot_data) == 0) {
    warning(paste("列", y_col, "没有有效的有限数值用于绘图，跳过。"))
    plot_list[[i]] <- ggplot() + annotate("text", x=0, y=0, label=paste("No valid data for", y_col)) + theme_void()
    
    # 保存无数据图
    safe_title <- gsub("[^a-zA-Z0-9_]", "_", p_title)
    filename <- file.path(png_output_dir, paste0("plot_", sprintf("%02d", i), "_", safe_title, ".png"))
    ggsave(filename, plot = plot_list[[i]], width = 8, height = 6, dpi = 300)
    cat("已保存无数据图:", filename, "\n")
    next
  }
  # --- 新增结束 ---
  
  # --- 根据 y 列类型设置 y 轴范围 ---
  if (grepl("RT", y_col, ignore.case = TRUE)) {
    y_scale_obj <- scale_y_continuous(expand = c(0.05, 0))
  } else if (grepl("Acc", y_col, ignore.case = TRUE) & !grepl("SC", y_col, ignore.case = TRUE)) {
    # 为原始准确率指标设置固定范围
    y_scale_obj <- scale_y_continuous(limits = c(0.65, 1), expand = c(0.01, 0))
  } else {
    # 对于其他指标 (如 SC, IES, Proportion) 使用默认扩展
    y_scale_obj <- scale_y_continuous(expand = c(0.05, 0))
  }
  
  # 构建 ggplot 对象 (不包含 geom_point)
  p <- ggplot(plot_data, aes(x = Age_year, y = !!sym(y_col))) +
    # 使用 GAM 进行平滑拟合，置信区间基于模型
    geom_smooth(method = "gam", formula = y ~ s(x, k = 3), 
                se = TRUE, 
                color = "black", 
                fill = "grey",   # 使用你喜欢的粉色
                alpha = 0.2) +      # 调整置信区间的透明度
    theme_minimal() +
    theme(
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 13, hjust = 0.5),
      panel.grid.minor = element_blank()
    ) +
    # 根据数据的实际范围设置 x 轴
    scale_x_continuous(breaks = seq(11, 18, by = 1), limits = c(11, 18)) +
    # 添加动态 y 轴 scale
    y_scale_obj +
    labs(title = p_title, x = "Age (years)", y = y_lab)
  
  plot_list[[i]] <- p
  
  # --- 新增：保存当前图表为 PNG ---
  safe_title <- gsub("[^a-zA-Z0-9_]", "_", p_title) 
  filename <- file.path(png_output_dir, paste0("plot_", sprintf("%02d", i), "_", safe_title, ".png"))
  ggsave(filename, plot = p, width = 8, height = 6, dpi = 300) 
  cat("已保存图表:", filename, "\n")
}

# --- 6. 保存到 PDF ---
output_pdf_path <- "D:/datasets/yunfu/raw_data/switch_results/all_switch_new_analysis_plots_no_points.pdf"

# 计算需要多少页 (每页12图 - 4行3列)
num_plots <- length(plot_list)
plots_per_page <- 12 # 4x3
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
cat("所有图表也已单独保存至目录:", png_output_dir, "\n")
