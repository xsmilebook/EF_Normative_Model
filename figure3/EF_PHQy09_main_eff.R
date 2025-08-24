# ===================================================================
#               优化后的 R 代码 (整合您的要求)
# ===================================================================

# --- 0. 加载所有必要的 R 包 ---
suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(tidyr); library(forcats)
  library(ggplot2); library(mgcv); library(purrr); library(multcomp)
  library(patchwork); library(readxl)
})

# --- 1. 文件路径设置 ---
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  # interfileFolder <- "/path/to/your/interfile_folder"
  # functionFolder <- "/path/to/your/functions"
  # resultFolder <- "/path/to/your/results/phq_y09_interaction_with_main_effect"
  # FigureFolder <- '/path/to/your/figures/phq_y09_fig_with_main_effect'
} else {
  interfileFolder <- "D:/datasets/yunfu/interfile_folder/temp"
  functionFolder <- "D:/code/EF_Normative_Model/functions"
  resultFolder <- "D:/datasets/yunfu/results/phq_y09_interaction_with_main_effect"
  FigureFolder <- 'D:/datasets/yunfu/figures/phq_y09_fig_with_main_effect'
}
dir.create(resultFolder, showWarnings = FALSE)
dir.create(FigureFolder, showWarnings = FALSE)


# --- 2. 数据读取与准备 ---
GNGd_data <- readRDS(paste0(interfileFolder, '/GNGd_prime.deviations.rds'))
back1_data <- readRDS(paste0(interfileFolder, '/back1Acc.deviations.rds'))
back2_data <- readRDS(paste0(interfileFolder, '/back2Acc.deviations.rds'))

# 读取 PHQ 数据
PHQ_data <- read_xlsx(paste0(interfileFolder, '/PHQ.xlsx')) %>%
  dplyr::select(用户ID, PHQ_y09) %>%
  dplyr::mutate(用户ID = as.character(用户ID))

# 数据类型统一和重命名ID
GNGd_data <- GNGd_data %>% dplyr::mutate(ID = as.character(x__ID), Sex = as.factor(Sex))
back1_data <- back1_data %>% dplyr::mutate(ID = as.character(x__ID), Sex = as.factor(Sex))
back2_data <- back2_data %>% dplyr::mutate(ID = as.character(x__ID), Sex = as.factor(Sex))

# 合并 PHQ 数据
GNGd_data <- GNGd_data %>% dplyr::inner_join(PHQ_data, by = c("ID" = "用户ID"))
back1_data <- back1_data %>% dplyr::inner_join(PHQ_data, by = c("ID" = "用户ID"))
back2_data <- back2_data %>% dplyr::inner_join(PHQ_data, by = c("ID" = "用户ID"))



# 将所有数据合并到一个数据框中，并创建任务因子
df_all <- dplyr::bind_rows(
  GNGd_data %>%
    dplyr::select(ID, Age_year, Sex, PHQ_y09, deviationZ = d_prime_deviationZ) %>%
    mutate(Task = "GNGd"),
  
  back1_data %>%
    dplyr::select(ID, Age_year, Sex, PHQ_y09, deviationZ = Oneback_acc_deviationZ) %>%
    mutate(Task = "back1"),
  
  back2_data %>%
    dplyr::select(ID, Age_year, Sex, PHQ_y09, deviationZ = Twoback_acc_deviationZ) %>%
    mutate(Task = "back2")
)

# 确保Task和PHQ_y09是因子，并设置正确的顺序
df_all$Task <- factor(df_all$Task, levels = c("GNGd", "back1", "back2"))
# 对于模型拟合，将PHQ_y09作为标准（无序）因子更稳健，模型会自动处理对比
df_all$PHQ_y09 <- factor(df_all$PHQ_y09, levels = c("0", "1", "2", "3"), ordered = FALSE)


# --- 3. 核心建模与分析函数 (修改后，包含主效应模型) ---

analyze_task_effects <- function(task_name, data) {
  
  cat(sprintf("\n--- Analyzing Task: %s ---\n", task_name))
  
  # 准备分析所需数据
  dft <- data %>% filter(Task == task_name) %>%
    filter(!is.na(deviationZ), !is.na(Age_year), !is.na(PHQ_y09), !is.na(Sex))
  
  if(nrow(dft) < 50) {
    cat("Warning: Insufficient data points. Skipping.\n")
    return(NULL)
  }
  
  k_use <- 3
  if (length(unique(dft$Age_year)) < k_use + 1) {
    cat("Warning: Not enough unique age points to fit a smooth. Skipping.\n")
    return(NULL)
  }
  
  # --- 定义两个模型公式 ---
  # 模型1: 包含年龄与PHQ的交互作用
  formula_interaction <- deviationZ ~ Sex + s(Age_year, k = k_use, fx = F) + s(Age_year, by = PHQ_y09, k = k_use, fx = F)
  # 模型2: 仅包含PHQ的主效应 
  formula_main <- deviationZ ~ Sex + s(Age_year, k = k_use, fx = F) + PHQ_y09
  
  # 拟合交互模型
  fit_interaction <- try(gam(formula_interaction, data = dft, method = "REML"), silent = TRUE)
  
  if (inherits(fit_interaction, "try-error")) {
    cat("ERROR: Interaction GAM fitting failed. Details:", as.character(fit_interaction), "\n")
    return(NULL)
  }
  
  cat("Interaction model fitted. Calculating effect size curves...\n")
  
  # --- 1. 计算效应大小曲线 (来自交互模型) ---
  ages <- seq(min(dft$Age_year), max(dft$Age_year), length.out = 100)
  sex_proportions <- table(dft$Sex) / nrow(dft)
  
  # 辅助函数：创建用于预测的新数据集
  make_newdata <- function(phq_level, sex_prop) {
    map_dfr(names(sex_prop), ~data.frame(
      Age_year = ages,
      PHQ_y09 = factor(phq_level, levels = levels(dft$PHQ_y09)),
      Sex = factor(.x, levels = levels(dft$Sex)),
      weight = sex_prop[.x]
    ))
  }
  
  # 获取基线(PHQ=0)的预测矩阵
  nd0 <- make_newdata("0", sex_proportions)
  Xp0 <- predict(fit_interaction, newdata = nd0, type = "lpmatrix")
  
  # 循环计算每个对比组 (1, 2, 3 vs 0) 的差异曲线
  effect_curves <- map_dfr(c("1", "2", "3"), function(level) {
    if (!level %in% levels(dft$PHQ_y09)) return(NULL)
    
    nd_level <- make_newdata(level, sex_proportions)
    Xp_level <- predict(fit_interaction, newdata = nd_level, type = "lpmatrix")
    
    # 计算对比矩阵 D, 并按性别比例加权
    D <- (Xp_level - Xp0) * nd_level$weight
    D_with_age <- as.data.frame(D)
    D_with_age$Age_year <- nd_level$Age_year
    
    # 按年龄聚合，消除性别的影响
    D_agg_df <- D_with_age %>%
      group_by(Age_year) %>%
      summarise(across(everything(), sum)) %>%
      ungroup()
    
    D_agg_matrix <- as.matrix(D_agg_df[, -1])
    
    beta <- coef(fit_interaction)
    V <- vcov(fit_interaction, unconditional = TRUE)
    
    # 计算效应大小和95%置信区间
    effect <- as.numeric(D_agg_matrix %*% beta)
    se <- sqrt(rowSums((D_agg_matrix %*% V) * D_agg_matrix))
    
    data.frame(
      Age_year = D_agg_df$Age_year,
      Contrast = paste("PHQ", level, "vs 0"),
      effect_size = effect,
      lower_ci = effect - 1.96 * se,
      upper_ci = effect + 1.96 * se
    )
  })
  
  # --- 2. 拟合主效应模型并提取系数 ---
  cat("Main effect model fitting...\n")
  fit_main <- try(gam(formula_main, data = dft, method = "REML"), silent = TRUE)
  
  main_effects_df <- NULL
  if (!inherits(fit_main, "try-error")) {
    main_summary <- summary(fit_main)
    p_table <- main_summary$p.table
    
    # 动态、安全地找到 PHQ 相关的行
    phq_rows <- rownames(p_table)[str_detect(rownames(p_table), "PHQ_y09")]
    
    if (length(phq_rows) > 0) {
      main_effects_df <- tibble(
        # 从行名中提取数字 (e.g., "PHQ_y091" -> "1")
        Contrast = paste("PHQ", str_extract(phq_rows, "\\d$"), "vs 0"),
        effect_size = p_table[phq_rows, "Estimate"]
      )
    }
    cat("Main effect coefficients extracted.\n")
  } else {
    cat("WARNING: Main effect GAM fitting failed.\n")
  }
  
  # 返回包含两种结果的列表
  return(list(effect_curves = effect_curves, main_effects = main_effects_df, fit_main = fit_main))
}

# --- 4. 运行所有任务的分析 ---
all_results <- map(c("GNGd", "back1", "back2"), ~analyze_task_effects(.x, data = df_all)) %>%
  set_names(c("GNGd", "back1", "back2"))

# --- 5. 可视化 (整合主效应和交互效应) ---
# 移除分析失败的空结果
all_results <- all_results[!sapply(all_results, is.null)]
if (length(all_results) == 0) {
  stop("All analyses failed. No results to plot.")
}

# 准备交互效应曲线的绘图数据
plot_data_curves <- map_dfr(names(all_results), ~all_results[[.x]]$effect_curves %>% mutate(Task = .x))

# 准备主效应基准线的绘图数据
plot_data_main <- map_dfr(names(all_results), ~all_results[[.x]]$main_effects %>% mutate(Task = .x))

# 定义任务和对比的标签及颜色，以确保绘图一致性
task_labels <- c("GNGd" = "Go/No-Go", "back1" = "1-back", "back2" = "2-back")
contrast_colors <- c("PHQ 1 vs 0" = "#4E79A7", "PHQ 2 vs 0" = "#F28E2B", "PHQ 3 vs 0" = "#E15759")

# 确保所有绘图数据的因子顺序和标签正确
plot_data_curves$Task <- factor(plot_data_curves$Task, levels = names(task_labels), labels = task_labels)
plot_data_curves$Contrast <- factor(plot_data_curves$Contrast, levels = names(contrast_colors))
if(nrow(plot_data_main) > 0){
  plot_data_main$Task <- factor(plot_data_main$Task, levels = names(task_labels), labels = task_labels)
  plot_data_main$Contrast <- factor(plot_data_main$Contrast, levels = names(contrast_colors))
}

# 创建最终的组合图
final_plot <- ggplot(plot_data_curves, aes(x = Age_year, y = effect_size, color = Contrast, fill = Contrast)) +
  # 0线
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth=0.5) +
  # 交互效应置信区间
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, linetype = 0) +
  # 交互效应曲线
  geom_line(linewidth = 1) +
  
  # --- 新增: 添加主效应基准线 ---
  geom_hline(
    data = plot_data_main,
    aes(yintercept = effect_size, color = Contrast),
    linetype = "dotted", # 使用点线以作区分
    linewidth = 1       # 设置线宽
  ) +
  
  # 按任务分面
  facet_wrap(~ Task, ncol = 3) +
  
  # 视觉元素设置
  scale_color_manual(values = contrast_colors, name = "Contrast") +
  scale_fill_manual(values = contrast_colors, name = "Contrast") +
  labs(
    title = "Effect of PHQ Levels (vs Baseline) on Cognitive Deviation Across Age",
    subtitle = "Solid lines are interaction effects with age; Dotted lines are main effects.",
    x = "Age (years)",
    y = "Effect Size (Difference in Deviation Z-score)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face="bold", size=16),
    plot.subtitle = element_text(size=12, color="grey30"),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold")
  )

# 在R环境中显示图片
print(final_plot)

# 保存图片
ggsave(
  filename = file.path(FigureFolder, "PHQ_y09_Effect_Curves_with_Main_Effect.pdf"),
  plot = final_plot,
  width = 12,
  height = 5.5, # 稍微增加高度以容纳副标题
  dpi = 300
)

cat(sprintf("\n--- Analysis complete. Plot saved to %s ---\n", FigureFolder))