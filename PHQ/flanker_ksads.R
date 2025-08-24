# ===================================================================
# 1. 设置与数据加载 (Setup & Data Loading)
# ===================================================================
library(openxlsx)
library(dplyr)
library(ggplot2)
library(effectsize)

library(ggpubr)

path  <- "D:/datasets/yunfu/interfile_folder/ABCD"
ksads_data <- read.csv(file.path(path, "new", "mh_y_ksads_ss.csv"))
# ksads_data <- read.csv("Z:/abcd-data-release-5.0/abcd-data-release-5.1/core/mental-health/mh_p_ksads_ss.csv")
flanker_data <- read.csv(file.path(path, "Flanker.deviations.csv"))


# ===================================================================
# 2. ID 创建与数据合并 (ID Creation & Merging)
# ===================================================================
ksads_prepared <- ksads_data %>%
  mutate(ksads_import_id_t = paste(src_subject_id, eventname, sep = "_"))

flanker_prepared <- flanker_data %>%
  rename(ksads_import_id_t = ID)

merged_data <- inner_join(
  ksads_prepared,
  flanker_prepared,
  by = "ksads_import_id_t"
)


# ===================================================================
# 3. 数据清洗与预分析 (Data Cleaning & Pre-analysis)
# ===================================================================
symptom_cols <- c("ksads_23_954_t", "ksads_23_965_t")
response_variable <- "nihtbx_flanker_uncorrected_deviationZ"

required_cols <- c(symptom_cols, response_variable)
if (!all(required_cols %in% names(merged_data))) {
  stop("错误：一个或多个指定的症状/响应列在合并后的数据中不存在。")
}

cat("--- 各症状列'1'的数量统计 (在删除前) ---\n")
counts_of_ones <- colSums(merged_data[symptom_cols] == 1, na.rm = TRUE)
print(counts_of_ones)

# 计算并报告将要删除的行数
cat("\n--- 缺失值分析 ---\n")
na_all_symptoms_count <- merged_data %>%
  filter(if_all(all_of(symptom_cols), is.na)) %>%
  nrow()
total_rows <- nrow(merged_data)
cat(sprintf("在 %d 总行数中，有 %d 行的全部5个症状列均为 NA。\n", 
            total_rows, 
            na_all_symptoms_count))

# 【修改】创建最终分析数据集：删除全NA行，并转换类型
analysis_data <- merged_data %>%
  # 步骤 1: 删除所有症状列均为 NA 的行
  filter(!if_all(all_of(symptom_cols), is.na)) %>%
  # 步骤 2: 将症状列转换为 numeric 类型
  mutate(across(all_of(symptom_cols), as.numeric))

# 【新增】报告删除操作的结果
cat(sprintf("\n已删除 %d 行。用于分析的数据集现在包含 %d 行。\n",
            na_all_symptoms_count,
            nrow(analysis_data)))


# ===================================================================
# 4. 特征工程：创建分组变量 (Feature Engineering)
# ===================================================================
analysis_data <- analysis_data %>%
  rowwise() %>%
  mutate(
    ksads = as.integer(1 %in% c_across(all_of(symptom_cols)))
  ) %>%
  ungroup() %>%
  mutate(
    ksads_group = factor(ksads, 
                         levels = c(0, 1), 
                         labels = c("Control Group (ksads=0)", "Symptom Group (ksads=1)"))
  )

cat("\n--- 创建的分组计数 (删除后) ---\n")
print(table(analysis_data$ksads_group))


# ===================================================================
# 5. 统计分析：两样本 T 检验 (Statistical Analysis)
# ===================================================================
ttest_result <- t.test(as.formula(paste(response_variable, "~ ksads")), 
                       data = analysis_data)

cat("\n\n--- 两样本 t 检验结果 (删除后) ---\n")
print(ttest_result)



# ===================================================================
# 1-4. (与之前相同) 数据加载、合并、清洗和特征工程
# 假设您已经运行了之前的代码，并得到了 'analysis_data' 数据框
# ===================================================================
# ... [前面的所有代码步骤] ...
# 最终得到的 analysis_data 包含 ksads, ksads_group 和 response_variable


# ===================================================================
# 5. 统计分析：T检验 与 科恩'd'效应量
# ===================================================================
# --- 5a. 两样本 T 检验 (同前) ---
response_variable <- "nihtbx_flanker_uncorrected_deviationZ"
ttest_result <- t.test(as.formula(paste(response_variable, "~ ksads")), 
                       data = analysis_data)

cat("\n\n--- 两样本 t 检验结果 ---\n")
print(ttest_result)

# --- 5b. 【新增】计算科恩'd' (Cohen's d) ---
cat("\n\n--- 科恩'd' 效应量计算结果 ---\n")
d_estimate <- cohens_d(as.formula(paste(response_variable, "~ ksads")), 
                       data = analysis_data)
print(d_estimate)


# ===================================================================
# 6. 【修改】结果可视化：箱线-散点组合图
# ===================================================================
box_jitter_plot <- ggplot(
  data = analysis_data, 
  aes(x = ksads_group, y = .data[[response_variable]], color = ksads_group)
) +
  # --- 添加图层 ---
  
  # 添加抖动散点图层，展示每个数据点
  # alpha设置透明度以显示数据密度
  # size设置点的大小，width控制水平抖动幅度
  geom_jitter(width = 0.25, alpha = 0.1, size = 0.8) +
  
  # 在散点图上层叠加箱线图
  # alpha=0.7使其半透明，可以看到后面的点
  # outlier.shape = NA 避免箱线图重复绘制离群点
  geom_boxplot(width = 0.2, alpha = 0.5, outlier.shape = NA) +
  
  # --- 添加统计检验结果 ---
  stat_compare_means(
    comparisons = list(c("Control Group (ksads=0)", "Symptom Group (ksads=1)")),
    method = "t.test",
    label = "p.format"
  ) +
  
  # --- 美化与自定义 ---
  labs(
    title = "Flanker Task Performance: Box Plot with Jitter",
    subtitle = "Comparison of Deviation Z-Scores (after removing NA rows)",
    x = "Group Status",
    y = "Deviation Z-Score"
  ) +
  
  # 使用 scale_color_manual 匹配颜色
  scale_color_manual(values = c("Control Group (ksads=0)" = "#0072B2", 
                                "Symptom Group (ksads=1)" = "#D55E00")) +
  
  # 使用简洁主题并进行微调
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major.x = element_blank(), # 移除垂直网格线
    panel.grid.minor.x = element_blank()
  )

# 打印新的图表
print(box_jitter_plot)

# ===================================================================
# 6. 结果可视化 (Visualization)
# ===================================================================
comparison_plot <- ggplot(
  data = analysis_data, 
  aes(x = ksads_group, y = .data[[response_variable]], fill = ksads_group)
) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  stat_compare_means(
    comparisons = list(c("Control Group (ksads=0)", "Symptom Group (ksads=1)")),
    method = "t.test",
    label = "p.format"
  ) +
  labs(
    title = "Flanker Task Performance by KSADS Symptom Group",
    subtitle = "Comparison of Deviation Z-Scores (after removing NA rows)",
    x = "Group Status",
    y = "Deviation Z-Score"
  ) +
  scale_fill_manual(values = c("Control Group (ksads=0)" = "#0072B2", 
                               "Symptom Group (ksads=1)" = "#D55E00")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(comparison_plot)