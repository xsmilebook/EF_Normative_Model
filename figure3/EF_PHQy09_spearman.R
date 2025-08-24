# 清空环境
rm(list = ls())

# --- 0. 加载必要的 R 包 ---
library(tidyverse)
library(readxl)
library(openxlsx) # 用于写入结果

# --- 1. 文件路径设置 ---
# (保持您原来的设置不变)
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  # ... 
} else {
  interfileFolder <- "D:/datasets/yunfu/interfile_folder/temp"
  resultFolder <- "D:/datasets/yunfu/results/phq_y09_interaction" # 结果可以存到同一个文件夹
}
dir.create(resultFolder, showWarnings = FALSE)

# --- 2. 数据读取与准备 ---
# (这部分代码与您提供的一致)

cat("--- Loading and preparing data... ---\n")
GNGd_data <- readRDS(paste0(interfileFolder, '/GNGd_prime.deviations.rds'))
back1_data <- readRDS(paste0(interfileFolder, '/back1Acc.deviations.rds'))
back2_data <- readRDS(paste0(interfileFolder, '/back2Acc.deviations.rds'))

PHQ_data <- read_xlsx(paste0(interfileFolder, '/PHQ.xlsx')) %>%
  dplyr::select(用户ID, PHQ_y09) %>%  
  dplyr::mutate(用户ID = as.character(用户ID))

GNGd_data <- GNGd_data %>% dplyr::mutate(ID = as.character(x__ID))
back1_data <- back1_data %>% dplyr::mutate(ID = as.character(x__ID))
back2_data <- back2_data %>% dplyr::mutate(ID = as.character(x__ID))

# 合并 PHQ 数据
GNGd_data <- GNGd_data %>% inner_join(PHQ_data, by = c("ID" = "用户ID"))
back1_data <- back1_data %>% inner_join(PHQ_data, by = c("ID" = "用户ID"))
back2_data <- back2_data %>% inner_join(PHQ_data, by = c("ID" = "用户ID"))
cat("--- Data loading complete. ---\n\n")

# --- 3. 定义变量 ---
EFvars.set <- data.frame(
  varname = c("d_prime_deviationZ", "Oneback_acc_deviationZ", "Twoback_acc_deviationZ"),
  dataname = c("GNGd_data", "back1_data", "back2_data"),
  task_label = c("Go/No-Go", "1-back", "2-back")
)

# --- 4. 初始化结果存储 ---
spearman_results_df <- data.frame()

# --- 5. 执行 Spearman 相关分析循环 ---

cat("--- Starting Spearman correlation analysis... ---\n")
for (i in 1:nrow(EFvars.set)) {
  
  # 获取当前任务的变量名和数据框名
  ef_var <- EFvars.set$varname[i]
  data_name <- EFvars.set$dataname[i]
  task_name <- EFvars.set$task_label[i]
  current_data <- get(data_name)
  
  cat(paste("  Analyzing Task:", task_name, "\n"))
  
  # 过滤掉包含 NA 的行，以确保相关性分析的准确性
  clean_data <- current_data %>%
    filter(!is.na(.data[[ef_var]]) & !is.na(PHQ_y09))
  
  # 检查是否有足够的数据进行分析
  if (nrow(clean_data) < 5) {
    cat(paste("    WARNING: Insufficient data for", task_name, ". Skipping.\n"))
    next # 跳过本次循环
  }
  
  # 执行 Spearman 相关性检验
  # 使用 "everything" 来处理数据中有并列排名的情况
  cor_test_result <- cor.test(
    x = clean_data[[ef_var]], 
    y = clean_data$PHQ_y09,
    method = "spearman",
    exact = FALSE # 使用近似值计算p值，对于大数据集更快且准确
  )
  
  # 提取结果
  rho_value <- cor_test_result$estimate
  p_value <- cor_test_result$p.value
  n_obs <- nrow(clean_data)
  
  # 将结果存入一个临时数据框
  temp_df <- data.frame(
    Task = task_name,
    EF_Variable = ef_var,
    PHQ_Variable = "PHQ_y09",
    Spearman_Rho = rho_value,
    P_Value = p_value,
    N = n_obs
  )
  
  # 合并到最终结果中
  spearman_results_df <- rbind(spearman_results_df, temp_df)
}
cat("--- Analysis complete. ---\n\n")

# --- 6. 报告并保存结果 ---

if (nrow(spearman_results_df) > 0) {
  
  # 在控制台打印美化后的结果
  cat("==================================================\n")
  cat("           Spearman Correlation Results           \n")
  cat("==================================================\n")
  print(spearman_results_df, row.names = FALSE)
  cat("==================================================\n\n")
  
  # 保存结果到 Excel 文件
  output_file <- paste0(resultFolder, "/PHQ_y09_spearman_correlation_results.xlsx")
  write.xlsx(spearman_results_df, file = output_file, row.names = FALSE)
  
  cat(paste("Results have been saved to:", output_file, "\n"))
  
} else {
  cat("--- No results were generated. ---\n")
}