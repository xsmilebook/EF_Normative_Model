# ===================================================================
# 1. 设置与库加载 (Setup & Library Loading)
# ===================================================================
# 确保已安装所有必要的包
library(openxlsx)
library(dplyr)
library(ggplot2)
library(effectsize)
library(MatchIt)   # 用于倾向性评分匹配
library(cobalt)    # 用于评估匹配后的平衡性

# 2. 数据加载与合并 (与您的代码相同)
path  <- "D:/datasets/yunfu/interfile_folder/ABCD"
ksads_data <- read.csv(file.path(path, "new", "mh_y_ksads_ss.csv"))
flanker_data <- read.csv(file.path(path, "Flanker.deviations.csv"))

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
# 3. 数据清洗与预分析 (与您的代码相同)
# ===================================================================
symptom_cols <- c("ksads_23_954_t", "ksads_23_965_t")
response_variable <- "nihtbx_flanker_uncorrected_deviationZ"
covariate_cols <- c("Age_year", "Sex") # 定义协变量

# 检查所需列是否存在
required_cols <- c(symptom_cols, response_variable, covariate_cols)
if (!all(required_cols %in% names(merged_data))) {
  stop("错误：一个或多个指定的症状/响应/协变量列在数据中不存在。")
}

# 删除全NA行，并转换类型
analysis_data <- merged_data %>%
  filter(!if_all(all_of(symptom_cols), is.na)) %>%
  mutate(across(all_of(symptom_cols), as.numeric)) %>%
  # 确保协变量没有缺失值，且类型正确
  filter(!is.na(Age_year) & !is.na(Sex)) %>%
  mutate(Sex = as.factor(Sex))

# ===================================================================
# 4. 特征工程：创建分组变量 (与您的代码相同)
# ===================================================================
analysis_data <- analysis_data %>%
  rowwise() %>%
  mutate(ksads = as.integer(1 %in% c_across(all_of(symptom_cols)))) %>%
  ungroup()

cat("\n--- 原始数据分组计数 ---\n")
print(table(analysis_data$ksads))

# ===================================================================
# 5. 【新增】倾向性评分匹配 (Propensity Score Matching)
# ===================================================================
cat("\n\n--- 正在执行倾向性评分匹配... ---\n")
match_formula <- as.formula("ksads ~ Age_year + Sex")

# 使用 MatchIt 包进行1:1最近邻匹配，并设置卡尺
match_obj <- matchit(
  match_formula,
  data = analysis_data,
  method = "nearest", 
  ratio = 1,          
  caliper = 0.2       
)

# 打印匹配摘要，显示匹配前后的样本量
cat("\n--- PSM 匹配摘要 ---\n")
print(summary(match_obj))

# ===================================================================
# 6. 【新增】评估匹配质量 (Assess Matching Quality)
# ===================================================================
cat("\n--- 正在生成协变量平衡图 (Love Plot)... ---\n")
bal_plot <- love.plot(
  match_obj, 
  binary = "std", 
  thresholds = c(m = .1), # 标出均衡性的阈值线
  title = "Covariate Balance Before and After Matching"
)
print(bal_plot)

# 提取匹配后的数据用于后续分析
matched_data <- match.data(match_obj)

# 7. 【修改】在匹配后的数据上进行统计分析
# --- 7a. 配对 t 检验 (Paired T-test) --- (最推荐的方法)
cat("\n\n--- 配对 t 检验结果 (在匹配数据上) ---\n")
treatment_group <- matched_data %>% filter(ksads == 1) %>% arrange(subclass)
control_group <- matched_data %>% filter(ksads == 0) %>% arrange(subclass)
ttest_result_psm <- t.test(treatment_group[[response_variable]], 
                           control_group[[response_variable]], 
                           paired = TRUE)
print(ttest_result_psm)


# --- 7c. 【新增演示】使用单样本t检验 (结果与配对t检验相同) ---
cat("\n\n--- (演示) 单样本 t 检验结果 (在配对差异上) ---\n")
# 1. 将数据重塑为宽格式，每行一个配对
paired_diff_data <- matched_data %>%
  select(subclass, ksads, !!sym(response_variable)) %>%
  tidyr::pivot_wider(
    id_cols = subclass,
    names_from = ksads,
    values_from = !!sym(response_variable),
    names_prefix = "group_"
  ) %>%
  # 2. 计算每对的差异
  mutate(difference = group_1 - group_0)

# 3. 对差异列执行单样本 t 检验，检验其均值是否为 0
ttest_one_sample <- t.test(paired_diff_data$difference, mu = 0)
print(ttest_one_sample)
cat("--- 注意：单样本t检验的t值和p值与上面的配对t检验完全相同。---\n")


# ===================================================================
# 8. 【修改】在匹配后的数据上进行结果可视化 (保持不变，仍使用配对t检验的结果)
# ===================================================================
matched_data <- matched_data %>%
  mutate(
    ksads_group = factor(ksads, 
                         levels = c(0, 1), 
                         labels = c("Control Group (ksads=0)", "Symptom Group (ksads=1)"))
  )

# 我们仍然使用更直接的配对t检验结果来标注图表
p_value_for_plot <- scales::pvalue(ttest_result_psm$p.value)

box_jitter_plot_psm <- ggplot(
  data = matched_data, 
  aes(x = ksads_group, y = .data[[response_variable]], color = ksads_group)
) +
  geom_jitter(width = 0.25, alpha = 0.2, size = 1) +
  geom_line(aes(group = subclass), color = "grey80", alpha = 0.5) +
  geom_boxplot(width = 0.2, alpha = 0.5, outlier.shape = NA) +
  annotate(
    geom = "text",
    x = 1.5, y = Inf, vjust = 1.5,
    label = paste("Paired t-test, p =", p_value_for_plot),
    size = 4
  ) +
  labs(
    title = "Flanker Task Performance After 1:1 PSM",
    subtitle = "Matched on Age and Sex",
    x = "Group Status",
    y = "Deviation Z-Score"
  ) +
  scale_color_manual(values = c("Control Group (ksads=0)" = "#0072B2", 
                                "Symptom Group (ksads=1)" = "#D55E00")) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

print(box_jitter_plot_psm)