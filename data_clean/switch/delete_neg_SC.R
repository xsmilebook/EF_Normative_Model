library(ggplot2)
library(dplyr)
library(readxl)
library(gridExtra)
library(writexl)  # 添加这个包用于写入Excel

# --- 1. 数据读取 ---
switch <- read_excel("D:/datasets/yunfu/raw_data/switch_results/switch_dependent_var_new.xlsx")

# --- 2. 删除 SC_SR_acc 为负数的行 ---
switch_cleaned <- switch %>%
  filter(SC_SP_acc <= 0)

# 查看删除后的数据维度
cat("删除负数行后剩余数据行数:", nrow(switch_cleaned), "\n")

# --- 3. 保存为新的Excel文件 ---
write_xlsx(switch_cleaned, "D:/datasets/yunfu/raw_data/switch_results/switch_dependent_var_new_del_neg_SP.xlsx")

cat("已保存清理后的数据到新Excel文件\n")