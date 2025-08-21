# 清空环境
rm(list = ls())

# --- 0. 加载所有必要的 R 包 ---
library(ggplot2)
library(tidyverse)
library(mgcv)
library(readxl)
library(openxlsx)
# library(gamm4) # <-- Not strictly needed for this model
# library(patchwork) # <-- Not used for this specific plot

# --- 1. 文件路径设置 ---
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  # ... 
} else {
  interfileFolder <- "D:/datasets/yunfu/interfile_folder/temp"
  functionFolder <- "D:/code/EF_Normative_Model/functions" # 
  resultFolder <- "D:/datasets/yunfu/results/phq_y09_interaction_contrasts"
  FigureFolder <- 'D:/datasets/yunfu/figures/phq_y09_fig'
}
dir.create(resultFolder, showWarnings = FALSE)
dir.create(FigureFolder, showWarnings = FALSE)

# --- 2. 加载自定义函数 ---
# Filename: gam.ordinal.contrasts.R
# Description: 
#   ADVANCED version. This function fits an ordinal GAM and then extracts the 
#   varying-coefficient (interaction) beta values for each of the R-1 contrasts 
#   (e.g., for logit(P(Y>1)), logit(P(Y>2)), etc.). It also performs
#   significance testing with Bonferroni correction.

# --- Required Libraries ---
library(mgcv)
library(tidyverse)
# library(MASS) # For mvrnorm

# --- Function Definition ---

gam.ordinal.contrasts <- function(
    dependentvar, 
    dataname, 
    smooth_var, 
    int_var, 
    covariates, 
    knots, 
    draws = 2000 # Increased draws for stable CIs
) {
  
  # --- 1. Model Fitting (using the no-main-effect version as requested) ---
  cat("\n--- Running Ordinal Contrasts GAM for:", dependentvar, "&", int_var, "---\n")
  gam.data <- get(dataname)
  gam.data <- gam.data %>% filter(!is.na(.data[[dependentvar]]) & !is.na(.data[[int_var]]))
  if (nrow(gam.data) < 50) { cat("!!! WARNING: Insufficient data. Skipping. !!!\n"); return(NULL) }
  
  # Ensure dependent var is integer and starts from 1
  gam.data[[dependentvar]] <- as.integer(gam.data[[dependentvar]])
  original_levels <- sort(unique(gam.data[[dependentvar]]))
  if (min(original_levels, na.rm = TRUE) == 0) {
    gam.data[[dependentvar]] <- gam.data[[dependentvar]] + 1
  }
  R <- length(unique(gam.data[[dependentvar]]))
  gam.family <- ocat(R = R)
  
  # Using the model formula without the main effect of int_var
  formula_full <- as.formula(sprintf("%s ~ s(%s, k=%d) + s(%s, by = %s, k=%d) + %s", 
                                     dependentvar, smooth_var, knots, smooth_var, int_var, knots, covariates))
  
  cat("Fitting model...\n")
  model_full <- gam(formula_full, data = gam.data, family = gam.family, method = "REML")
  
  # --- 2. Prepare for Post-processing ---
  cat("Extracting contrast-specific beta coefficients...\n")
  
  # Get model coefficients and variance-covariance matrix
  b <- coef(model_full)
  V <- vcov(model_full, unconditional = TRUE)
  
  # Posterior simulation of coefficients
  sim_b <- MASS::mvrnorm(draws, mu = b, Sigma = V)
  
  # Identify the interaction smooth term
  interaction_smooth_label <- sprintf("s(%s):%s", smooth_var, int_var)
  smooth_idx <- which(sapply(model_full$smooth, function(s) s$label) == interaction_smooth_label)
  if (length(smooth_idx) == 0) {
    cat("!!! ERROR: Could not find the interaction smooth term. Skipping. !!!\n")
    return(NULL)
  }
  
  # Get parameter indices for this smooth for the FIRST contrast
  para_indices <- model_full$smooth[[smooth_idx]]$first.para:model_full$smooth[[smooth_idx]]$last.para
  
  # Get number of coefficients per contrast
  p_per_contrast <- length(b) / (R - 1)
  
  # Create a prediction grid to evaluate the smooth function f(age)
  # We set the 'by' variable to 1 to get the basis functions for the smooth
  smooth_range <- range(gam.data[[smooth_var]])
  pred_grid <- data.frame(value = seq(smooth_range[1], smooth_range[2], length.out = 100))
  names(pred_grid) <- smooth_var
  pred_grid[[int_var]] <- 1 # Set 'by' variable to 1
  
  # Handle covariates by setting them to their mode
  if (!is.null(covariates) && covariates != "") {
    gam.data[[covariates]] <- as.factor(gam.data[[covariates]])
    mode_level <- names(which.max(table(gam.data[[covariates]])))
    pred_grid[[covariates]] <- factor(mode_level, levels = levels(gam.data[[covariates]]))
  }
  
  # Get the linear predictor matrix (Xp)
  Xp <- predict(model_full, newdata = pred_grid, type = "lpmatrix")
  
  # Extract the columns of Xp corresponding to the interaction smooth
  Xp_int <- Xp[, para_indices, drop = FALSE]
  
  # --- 3. Calculate Betas for Each Contrast ---
  all_contrasts_results <- list()
  
  for (j in 1:(R - 1)) {
    # The contrast label refers to the original data levels
    # ocat models P(Y > j). If original data is 0,1,2,3 -> new is 1,2,3,4.
    # j=1 -> P(Y>1) -> P(orig Y > 0) -> "PHQ 1 vs 0"
    # j=2 -> P(Y>2) -> P(orig Y > 1) -> "PHQ 2 vs 0/1"
    # j=3 -> P(Y>3) -> P(orig Y > 2) -> "PHQ 3 vs 0/1/2"
    # We will use the labels from the example image for simplicity
    contrast_label <- paste("PHQ", j, "vs 0") 
    
    # Calculate index offset for the j-th contrast
    offset <- (j - 1) * p_per_contrast
    current_indices <- para_indices + offset
    
    # Extract simulated coefficients for this specific contrast's interaction
    sim_b_int_j <- sim_b[, current_indices, drop = FALSE]
    
    # Calculate the posterior distribution of the beta coefficient (which is f(age))
    # This is a matrix of [n_age_points x n_draws]
    posterior_betas <- Xp_int %*% t(sim_b_int_j)
    
    # Summarize the posterior distribution for each age point
    contrast_df <- data.frame(
      Age = pred_grid[[smooth_var]],
      Contrast = contrast_label,
      beta = apply(posterior_betas, 1, median),
      lower_ci = apply(posterior_betas, 1, quantile, probs = 0.025),
      upper_ci = apply(posterior_betas, 1, quantile, probs = 0.975)
    )
    all_contrasts_results[[j]] <- contrast_df
  }
  
  # Combine all results into a single dataframe
  final_results <- do.call(rbind, all_contrasts_results)
  
  # --- 4. Perform Per-Task Bonferroni Correction for Significance ---
  n_tests <- nrow(final_results) # Number of age points * number of contrasts
  alpha <- 0.05 / n_tests
  
  # Recalculate CIs and determine significance based on the corrected alpha
  # This is computationally intensive, so we re-use the posterior samples
  
  corrected_results_list <- list()
  for (j in 1:(R - 1)) {
    # Re-extract the posterior betas for this contrast
    offset <- (j - 1) * p_per_contrast
    current_indices <- para_indices + offset
    sim_b_int_j <- sim_b[, current_indices, drop = FALSE]
    posterior_betas <- Xp_int %*% t(sim_b_int_j)
    
    # Find the corrected CIs
    lower_ci_corr <- apply(posterior_betas, 1, quantile, probs = alpha / 2)
    upper_ci_corr <- apply(posterior_betas, 1, quantile, probs = 1 - (alpha / 2))
    
    # Filter the original results and add the significance column
    temp_df <- final_results %>% 
      filter(Contrast == paste("PHQ", j, "vs 0")) %>%
      mutate(significant = (lower_ci_corr > 0) | (upper_ci_corr < 0))
    
    corrected_results_list[[j]] <- temp_df
  }
  
  final_results_corrected <- do.call(rbind, corrected_results_list)
  
  cat("--- Contrast analysis complete ---\n")
  return(final_results_corrected)
}


# --- 3. 数据读取与准备 (带诊断的版本) ---

# 读取行为任务数据
cat("--- Loading RDS files ---\n")
GNGd_data <- readRDS(paste0(interfileFolder, '/GNGd_prime.deviations.rds'))
back1_data <- readRDS(paste0(interfileFolder, '/back1Acc.deviations.rds'))
back2_data <- readRDS(paste0(interfileFolder, '/back2Acc.deviations.rds'))
cat(sprintf("Loaded GNGd_data with %d rows.\n", nrow(GNGd_data)))

# 读取 PHQ 数据
cat("\n--- Loading PHQ.xlsx file ---\n")
PHQ_data <- read_xlsx(paste0(interfileFolder, '/PHQ.xlsx')) %>%
  select(用户ID, PHQ_y09) %>%  
  mutate(用户ID = as.character(用户ID))
cat(sprintf("Loaded PHQ_data with %d rows.\n", nrow(PHQ_data)))


# --- 关键检查点：在合并前，检查两个数据框的 ID ---
cat("\n\n--- CRITICAL CHECK: Inspecting IDs before join ---\n")

# 准备行为数据，但【不】合并
GNGd_data_pre_join <- GNGd_data %>% mutate(ID = as.character(x__ID), Gender = as.factor(Sex))

# 打印双方ID的前6个，用肉眼检查格式！
cat("First 6 IDs from GNGd_data (column 'ID'):\n")
print(head(GNGd_data_pre_join$ID))

cat("\nFirst 6 IDs from PHQ_data (column '用户ID'):\n")
print(head(PHQ_data$用户ID))

# 你需要检查：
# 1. 是不是一个是纯数字，另一个是 "yunfu_xxxx" 这样的格式？
# 2. 是不是有隐藏的空格？
# 3. 是不是有 .0 这样的后缀？
# 4. 是不是一个有前导零 (e.g., "0123") 一个没有 (e.g., "123")？

# --- 执行合并，并立即检查结果 ---
cat("\n--- Performing the inner_join ---\n")
GNGd_data_after_join <- GNGd_data_pre_join %>% 
  inner_join(PHQ_data, by = c("ID" = "用户ID"))

# 【决定性证据】 打印合并后的数据框行数
cat("\n\n--- MOMENT OF TRUTH ---\n")
cat(sprintf("Number of rows in GNGd_data AFTER inner_join: %d\n", nrow(GNGd_data_after_join)))
cat("-----------------------\n\n")

# 如果上面的行数是 0，你就找到了问题的根源：ID 完全不匹配。
# 如果行数 > 0，问题可能在于列名，但这与你之前的诊断相矛盾。

# 为保证脚本能继续运行，我们重新赋值
GNGd_data <- GNGd_data %>% inner_join(PHQ_data, by = c("ID" = "用户ID"))
back1_data <- back1_data %>% inner_join(PHQ_data, by = c("ID" = "用户ID"))
back2_data <- back2_data %>% inner_join(PHQ_data, by = c("ID" = "用户ID"))
# --- 4. 定义变量与参数 ---
EFvars.set <- as.data.frame(matrix(c(
  "d_prime_deviationZ", "GNGd",
  "Oneback_acc_deviationZ", "back1",
  "Twoback_acc_deviationZ", "back2"
), byrow = TRUE, ncol = 2, dimnames = list(NULL, c("varname", "dataname"))))
knots <- 3

# --- 5. 初始化结果存储 ---
all_tasks_contrasts <- list()


cat("\n--- DIAGNOSTIC CHECK: Checking available data for each task ---\n")

# 检查 GNGd
valid_gngd <- GNGd_data %>% 
  filter(!is.na(PHQ_y09) & !is.na(d_prime_deviationZ))
cat(sprintf("Task: GNGd | Valid observations (PHQ_y09 & d_prime_deviationZ): %d\n", nrow(valid_gngd)))

# 检查 back1
valid_back1 <- back1_data %>%
  filter(!is.na(PHQ_y09) & !is.na(Oneback_acc_deviationZ))
cat(sprintf("Task: back1 | Valid observations (PHQ_y09 & Oneback_acc_deviationZ): %d\n", nrow(valid_back1)))

# 检查 back2
valid_back2 <- back2_data %>%
  filter(!is.na(PHQ_y09) & !is.na(Twoback_acc_deviationZ))
cat(sprintf("Task: back2 | Valid observations (PHQ_y09 & Twoback_acc_deviationZ): %d\n", nrow(valid_back2)))

cat("--- END OF DIAGNOSTIC CHECK ---\n\n")

# --- 6. 主循环：执行 GAM 对比度分析 ---
for (i in 1:nrow(EFvars.set)) {
  int_var <- EFvars.set$varname[i]
  dataname_str <- EFvars.set$dataname[i]
  dataname <- paste0(dataname_str, "_data")
  
  # 调用新的对比度分析函数
  result_df <- gam.ordinal.contrasts(
    dependentvar = "PHQ_y09", 
    dataname = dataname, 
    smooth_var = "Age_year", 
    int_var = int_var, 
    covariates = "Gender",
    knots = knots
  )
  
  if (!is.null(result_df)) {
    all_tasks_contrasts[[dataname_str]] <- result_df
  }
}

# --- 7. 合并并保存结果 ---
if (length(all_tasks_contrasts) > 0) {
  # 使用 bind_rows 并通过 .id 参数添加任务名称列
  final_contrast_df <- bind_rows(all_tasks_contrasts, .id = "Task")
  
  # 保存详细的数值结果
  write.xlsx(final_contrast_df, file = paste0(resultFolder, "/PHQ_y09_all_tasks_contrasts.xlsx"), row.names = FALSE)
} else {
  stop("--- ERROR: Analysis failed for all tasks. No results to plot. ---")
}


# --- 8. 可视化：生成 Beta 热力图 ---

# 准备绘图数据
plot_data <- final_contrast_df %>%
  mutate(
    # 将任务和对比度转换为因子，以控制绘图顺序
    Task = factor(Task, levels = c("back1", "back2", "GNGd")),
    Contrast = factor(Contrast, levels = c("PHQ 3 vs 0", "PHQ 2 vs 0", "PHQ 1 vs 0"))
  )

# 找到显著区域用于绘制黑框
significant_boxes <- plot_data %>%
  filter(significant == TRUE) %>%
  group_by(Task, Contrast) %>%
  # 对于每个有显著点的任务-对比度组合，计算其年龄范围
  summarise(
    xmin = min(Age),
    xmax = max(Age),
    .groups = 'drop'
  ) %>%
  # 将因子转换为数值以便计算矩形框的 y 范围
  mutate(
    y_center = as.numeric(Contrast),
    ymin = y_center - 0.5,
    ymax = y_center + 0.5
  )

# 确定颜色范围的对称极值
max_abs_beta <- max(abs(plot_data$beta), na.rm = TRUE)

# --- 开始绘图 ---
beta_heatmap <- ggplot(plot_data, aes(x = Age, y = Contrast)) +
  # 1. 绘制热力图瓦片
  geom_tile(aes(fill = beta)) +
  
  # 2. 添加显著区域的黑框
  geom_rect(
    data = significant_boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    color = "black",   # 框线颜色
    fill = NA,         # 透明填充
    linewidth = 0.6,
    inherit.aes = FALSE # 不继承上层的 aesthetics
  ) +
  
  # 3. 设置颜色映射
  scale_fill_gradient2(
    name = expression(beta),
    low = "#0571b0",      # 负值颜色 (蓝)
    mid = "white",        # 零值颜色
    high = "#ca0020",     # 正值颜色 (红)
    midpoint = 0,
    limits = c(-max_abs_beta, max_abs_beta),
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
  ) +
  
  # 4. 按任务分面
  facet_wrap(~Task, ncol = 1, scales = "free_y") +
  
  # 5. 美化主题和标签
  labs(
    title = "Beta heatmap (per-task Bonferroni; significant regions boxed)",
    x = "Age (years)",
    y = "Contrast"
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(), # 移除分面标签的背景
    strip.text = element_text(face = "bold", size=14, hjust = 0),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.key.width = unit(0.5, "cm"),
    axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=14)
  )

# 打印并保存图形
print(beta_heatmap)

ggsave(
  filename = paste0(FigureFolder, "/PHQ_y09_Beta_Heatmap_Plot.pdf"),
  plot = beta_heatmap,
  width = 8,
  height = 9,
  dpi = 300
)

cat("\n--- Beta heatmap plot saved successfully. ---\n")