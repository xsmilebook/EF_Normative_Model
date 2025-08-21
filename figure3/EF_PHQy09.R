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
# (保持您的路径设置不变)
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  # ... 
} else {
  interfileFolder <- "D:/datasets/yunfu/interfile_folder/temp"
  functionFolder <- "D:/code/EF_Normative_Model/functions"
  resultFolder <- "D:/datasets/yunfu/results/phq_y09_interaction_optimized"
  FigureFolder <- 'D:/datasets/yunfu/figures/phq_y09_fig_optimized'
}
dir.create(resultFolder, showWarnings = TRUE) # showWarnings=TRUE 可以在创建时看到提示
dir.create(FigureFolder, showWarnings = TRUE)

# --- 2. 数据读取与准备 ---
# (保持您的数据读取和清洗代码不变)
# 读取行为任务数据
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
df_all$PHQ_y09 <- factor(df_all$PHQ_y09, levels = c("0", "1", "2", "3"), ordered = TRUE) # 使用有序因子更符合变量性质

# --- 3. 核心建模与分析函数 (无备用模型) ---

analyze_task_interaction <- function(task_name, data) {
  
  cat(sprintf("\n--- Analyzing Task: %s ---\n", task_name))
  
  dft <- data %>% filter(Task == task_name) %>%
    filter(!is.na(deviationZ), !is.na(Age_year), !is.na(PHQ_y09), !is.na(Sex))
  
  if(nrow(dft) < 50) {
    cat("Warning: Insufficient data points. Skipping.\n")
    return(NULL)
  }
  
  # 动态调整k值，防止k值过大导致模型无法拟合
  k_use <- 3
  if (k_use < 3) {
    cat("Warning: Not enough unique age points to fit a smooth. Skipping.\n")
    return(NULL)
  }
  
  # 定义模型公式，PHQ=0为参考水平
  formula_interaction <- deviationZ ~ Sex + s(Age_year, k = k_use, fx = F) + s(Age_year, by = PHQ_y09, k = k_use, fx = F)
  formula_main <- deviationZ ~ Sex + s(Age_year, k = k_use, fx = F) + PHQ_y09
  # 拟合GAM模型
  fit <- try(gam(formula_interaction, data = dft, method = "REML"), silent = TRUE)
  
  if (inherits(fit, "try-error")) {
    cat("ERROR: GAM fitting failed. Details:", as.character(fit), "\n")
    return(NULL)
  }
  
  cat("Model fitted successfully. Now calculating effect size curves...\n")
  
  # --- 计算效应大小曲线 (1v0, 2v0, 3v0) ---
  ages <- seq(min(dft$Age_year), max(dft$Age_year), length.out = 100)
  
  # 这个辅助函数用于生成预测所需的新数据
  make_newdata <- function(phq_level, sex_prop) {
    # 对性别进行边际化（加权平均）
    map_dfr(names(sex_prop), ~data.frame(
      Age_year = ages,
      PHQ_y09 = factor(phq_level, levels = levels(dft$PHQ_y09)),
      Sex = factor(.x, levels = levels(dft$Sex)),
      weight = sex_prop[.x]
    ))
  }
  
  sex_proportions <- table(dft$Sex) / nrow(dft)
  
  # 获取基线(PHQ=0)的预测矩阵
  nd0 <- make_newdata("0", sex_proportions)
  Xp0 <- predict(fit, newdata = nd0, type = "lpmatrix")
  
  # 循环计算每个对比组的差异
  effect_curves <- map_dfr(c("1", "2", "3"), function(level) {
    if (!level %in% levels(dft$PHQ_y09)) return(NULL) 
    
    nd_level <- make_newdata(level, sex_proportions)
    Xp_level <- predict(fit, newdata = nd_level, type = "lpmatrix")
    
    # 对比矩阵 D = X_level - X_0, 同时考虑性别权重
    D <- (Xp_level - Xp0) * nd_level$weight
    
    # --- *** MODIFICATION START *** ---
    # 使用 dplyr 进行聚合，这比 aggregate 更稳健
    D_with_age <- as.data.frame(D)
    D_with_age$Age_year <- nd_level$Age_year
    
    D_agg_df <- D_with_age %>%
      group_by(Age_year) %>%
      summarise(across(everything(), sum)) %>%
      ungroup()
    
    # 将聚合后的结果转换为纯数值矩阵用于乘法
    D_agg_matrix <- as.matrix(D_agg_df[, -1])
    # --- *** MODIFICATION END *** ---
    
    beta <- coef(fit)
    V <- vcov(fit, unconditional = TRUE)
    
    # 计算效应大小 (点估计) 和标准误
    effect <- as.numeric(D_agg_matrix %*% beta)
    se <- sqrt(rowSums((D_agg_matrix %*% V) * D_agg_matrix)) # 注意这里的矩阵乘法也需要用 D_agg_matrix
    
    data.frame(
      Age_year = D_agg_df$Age_year, # 从聚合后的数据框中获取年龄
      Contrast = paste("PHQ", level, "vs 0"),
      effect_size = effect,
      lower_ci = effect - 1.96 * se,
      upper_ci = effect + 1.96 * se
    )
  })
  
  return(list(fit = fit, effect_curves = effect_curves))
}

# --- 4. 运行所有任务的分析 ---
all_results <- map(c("GNGd", "back1", "back2"), ~analyze_task_interaction(.x, data = df_all)) %>%
  set_names(c("GNGd", "back1", "back2"))

# --- 5. 可视化 ---
# 移除空结果
all_results <- all_results[!sapply(all_results, is.null)] 
if (length(all_results) == 0) {
  stop("All analyses failed. No results to plot.")
}

# 准备绘图数据
plot_data <- map_dfr(names(all_results), ~all_results[[.x]]$effect_curves %>% mutate(Task = .x))

# 定义任务和对比的标签及颜色
task_labels <- c("GNGd" = "Go/No-Go", "back1" = "1-back", "back2" = "2-back")
contrast_colors <- c("PHQ 1 vs 0" = "#4E79A7", "PHQ 2 vs 0" = "#F28E2B", "PHQ 3 vs 0" = "#E15759")
plot_data$Task <- factor(plot_data$Task, levels = names(task_labels), labels = task_labels)
plot_data$Contrast <- factor(plot_data$Contrast, levels = names(contrast_colors))

# 创建最终的组合图
final_plot <- ggplot(plot_data, aes(x = Age_year, y = effect_size, color = Contrast, fill = Contrast)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, linetype = 0) + # 移除ribbon的边线
  geom_line(linewidth = 1) +
  facet_wrap(~ Task, ncol = 3) + # 每个任务一个面板
  scale_color_manual(values = contrast_colors, name = "Contrast") +
  scale_fill_manual(values = contrast_colors, name = "Contrast") +
  labs(
    title = "Effect Size of PHQ Levels (vs Baseline) on Cognitive Deviation Across Age",
    x = "Age (years)",
    y = "Effect Size (Difference in Deviation Z-score)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold")
  )

print(final_plot)

# 保存图片
ggsave(
  filename = file.path(FigureFolder, "PHQ_y09_Effect_Size_Curves_Optimized_no_main_eff.pdf"),
  plot = final_plot,
  width = 12,
  height = 5,
  dpi = 300
)

cat(sprintf("\n--- Analysis complete. Plot saved to %s ---\n", FigureFolder))