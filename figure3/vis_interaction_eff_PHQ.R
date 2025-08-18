
# --- 0. 加载所有必要的 R 包 ---
library(ggplot2)
library(tidyverse)
library(mgcv)
library(readxl)
library(openxlsx)
library(cowplot) # Using cowplot for plot_grid, similar to the reference script

# --- 1. 文件路径设置 ---
# PLEASE UPDATE THESE PATHS TO MATCH YOUR FOLDER STRUCTURE
datapath <- 'D:/datasets/yunfu/raw_data'
interfileFolder <- "D:/datasets/yunfu/interfile_folder/temp"
functionFolder <- "D:/code/EF_Normative_Model/functions"
resultFolder <- "D:/datasets/yunfu/results/phq_sum_corr_modified"
FigureFolder <- 'D:/datasets/yunfu/figures/phq_sum_fig_modified'
# This is the new file you specified for the beta values
beta_file_path <- "D:/datasets/yunfu/results/psy_corr/PHQ_all_deviation_correlations_new.xlsx"

dir.create(resultFolder, showWarnings = FALSE)
dir.create(FigureFolder, showWarnings = FALSE)

# --- 2. 数据读取与预处理 ---

# 读取行为任务数据
GNGd_data <- readRDS(paste0(interfileFolder, '/GNGd_prime.deviations.rds'))
back1_data <- readRDS(paste0(interfileFolder, '/back1Acc.deviations.rds'))
back2_data <- readRDS(paste0(interfileFolder, '/back2Acc.deviations.rds'))

# 读取 PHQ 数据
PHQ_data <- read_xlsx(paste0(interfileFolder, '/PHQ.xlsx')) %>%
  select(用户ID, PHQ_sum) %>%
  mutate(用户ID = as.character(用户ID))

# 读取包含 Beta 值的相关性结果文件
beta_table <- read_xlsx(beta_file_path)

# 数据类型统一
GNGd_data <- GNGd_data %>% mutate(x__ID = as.character(x__ID), Sex = as.factor(Sex))
back1_data <- back1_data %>% mutate(x__ID = as.character(x__ID), Sex = as.factor(Sex))
back2_data <- back2_data %>% mutate(x__ID = as.character(x__ID), Sex = as.factor(Sex))

# 合并 PHQ 数据
GNGd_data <- GNGd_data %>% inner_join(PHQ_data, by = c("x__ID" = "用户ID"))
back1_data <- back1_data %>% inner_join(PHQ_data, by = c("x__ID" = "用户ID"))
back2_data <- back2_data %>% inner_join(PHQ_data, by = c("x__ID" = "用户ID"))

# 检查合并
if(nrow(GNGd_data) == 0) stop("ERROR: GNGd_data is empty after merging. Check ID columns.")

# 加载自定义函数
source(paste0(functionFolder, '/gam_varyingcoefficients.R'))

# --- 3. 定义变量与参数 ---
psyc_variables_continous <- c("PHQ_sum")

EFvars.set <- as.data.frame(matrix(c(
  "d_prime_deviationZ", "GNGd",
  "Oneback_acc_deviationZ", "back1",
  "Twoback_acc_deviationZ", "back2"
), byrow = TRUE, ncol = 2, dimnames = list(NULL, c("varname", "dataname"))))

# 对因变量进行 Z-score 标准化
GNGd_data[, "PHQ_sum_z"] <- scale(GNGd_data[, "PHQ_sum"])
back1_data[, "PHQ_sum_z"] <- scale(back1_data[, "PHQ_sum"])
back2_data[, "PHQ_sum_z"] <- scale(back2_data[, "PHQ_sum"])

psyc_variables_continous_z <- "PHQ_sum_z"


# --- 可视化验证代码 ---

# 确保 back1_data 已经加载并且包含了 PHQ_sum_z 和 Oneback_acc_deviationZ
# 如果没有 PHQ_sum_z, 您可以先运行之前的标准化代码，或者直接使用 PHQ_sum

# 1. 创建年龄分组变量
# 我们将年龄分为三组：11-13 (早期), 14-15 (中期), 16-18 (晚期)
back1_validation_data <- back1_data %>%
  # 剔除必要的缺失值以保证绘图成功
  filter(!is.na(PHQ_sum_z), !is.na(Oneback_acc_deviationZ), !is.na(Age_year)) %>%
  mutate(
    AgeGroup = case_when(
      Age_year <= 13 ~ "1. 青少年早期 (11-13岁)",
      Age_year > 13 & Age_year <= 15 ~ "2. 青少年中期 (14-15岁)",
      Age_year > 15 ~ "3. 青少年晚期 (16-18岁)"
    )
  )

# 2. 使用 ggplot2 创建分面散点图
validation_plot <- ggplot(back1_validation_data, 
                          aes(x = Oneback_acc_deviationZ, y = PHQ_sum_z)) +
  # 绘制散点
  geom_point(alpha = 0.4, color = "skyblue") +
  # 添加线性回归趋势线 (lm = linear model)
  geom_smooth(method = "lm", se = TRUE, color = "navy") +
  # 按照年龄组分面
  facet_wrap(~ AgeGroup, ncol = 3) +
  # 添加标签和主题
  labs(
    title = "按年龄组检验：PHQ分数与1-back表现的关系",
    x = "1-back 表现 (Z-score偏差)",
    y = "PHQ 总分 (Z-score)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    strip.text = element_text(size = 11, face = "bold") # 设置分面标题样式
  )

# 3. 打印图形
print(validation_plot)




# --- 全新可视化验证代码 ---

# 确保 back1_data 已经加载并且包含了必要的数据列

# 1. 创建认知表现分组
# 我们使用中位数将 1-back 表现分为“高分组”和“低分组”
back1_validation_data_new <- back1_data %>%
  filter(!is.na(PHQ_sum_z), !is.na(Oneback_acc_deviationZ), !is.na(Age_year)) %>%
  mutate(
    performance_median = median(Oneback_acc_deviationZ, na.rm = TRUE),
    Oneback_Group = ifelse(Oneback_acc_deviationZ >= performance_median,
                           "高表现组",
                           "低表现组")
  )

# 2. 使用 ggplot2 绘制两组的年龄趋势线
new_validation_plot <- ggplot(back1_validation_data_new, 
                              aes(x = Age_year, y = PHQ_sum_z, color = Oneback_Group)) +
  # 绘制平滑的年龄趋势线 (使用 GAM 平滑) 和置信区间
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE) +
  # 也可以选择性地加上散点
  # geom_point(alpha = 0.3) +
  # 设置颜色
  scale_color_manual(values = c("高表现组" = "darkgreen", "低表现组" = "darkred")) +
  # 添加标签和主题
  labs(
    title = "条件效应可视化：不同1-back表现组的PHQ-年龄趋势",
    subtitle = "如果曲线不平行，则证明存在交互效应",
    x = "年龄 (岁)",
    y = "PHQ 总分 (Z-score)",
    color = "1-back 表现分组"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "top"
  )

# 3. 打印图形
print(new_validation_plot)





# --- 验证主效应 Beta 值的可视化代码 ---

# 确保 back1_data 数据已加载
# 注意：请确保协变量与您在 plot_PHQ.R 中使用的模型一致
# 例如，如果原始数据中使用的是 Gender，请将下面代码中的 Sex 改为 Gender

# 1. 准备数据，剔除缺失值
validation_data_main_effect <- back1_data %>%
  filter(!is.na(PHQ_sum_z), !is.na(Oneback_acc_deviationZ), !is.na(Age_year), !is.na(Sex))

# 2. 计算部分残差
# a. 计算 PHQ 中无法被 Age 和 Sex 解释的部分
phq_residual_model <- gam(PHQ_sum_z ~ s(Age_year) + Sex, data = validation_data_main_effect)
validation_data_main_effect$phq_partial_residual <- residuals(phq_residual_model)

# b. 计算 1-back 中无法被 Age 和 Sex 解释的部分
oneback_residual_model <- gam(Oneback_acc_deviationZ ~ s(Age_year) + Sex, data = validation_data_main_effect)
validation_data_main_effect$oneback_partial_residual <- residuals(oneback_residual_model)

# 3. 绘制部分残差图
main_effect_validation_plot <- ggplot(validation_data_main_effect, 
                                      aes(x = oneback_partial_residual, y = phq_partial_residual)) +
  geom_point(alpha = 0.5, color = "gray50") +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick", fill = "pink", alpha = 0.3) +
  labs(
    title = "主效应 Beta 值可视化 (部分残差图)",
    subtitle = "趋势线的斜率 ≈ GAM模型计算出的总体Beta值",
    x = "1-back 表现的残差 (已控制年龄和性别)",
    y = "PHQ分数的残差 (已控制年龄和性别)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )

# 4. 打印图形
print(main_effect_validation_plot)

# (可选) 额外进行数值验证：
# 计算残差之间的线性模型，查看其系数
numerical_validation <- lm(phq_partial_residual ~ oneback_partial_residual, data = validation_data_main_effect)
cat("\n--- 数值验证 ---\n")
cat("部分残差图趋势线的斜率 (beta值) 约为:\n")
print(coef(numerical_validation)[2])
cat("您可以将此值与您在 PHQ_all_deviation_correlations_new.xlsx 文件中为 back1 找到的 beta 值进行比较。\n")



