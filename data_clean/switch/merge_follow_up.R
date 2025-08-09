# 加载必要的库
library(readxl)
library(writexl)
library(dplyr)

# 读取数据
switch <- read_excel("D:/datasets/yunfu/raw_data/switch_results/switch_dependent_followup.xlsx")
table2 <- read_excel('D:/datasets/yunfu/raw_data/yunfu_agebyID.xlsx')

# 选择table2中的指定列
table2_sub <- table2[, c('ID', 'Name', 'Age_day', 'Age_month', 'Age_year')]

# 合并数据（相当于MATLAB的outerjoin）
combinetable <- left_join(switch, table2_sub, by = c('ID'), suffix = c('', '.y'))

# 处理重复列名（如果有.x和.y后缀的列，保留原始列）
# 移除可能产生的重复列
combinetable <- combinetable[, !grepl("\\.y$", names(combinetable))]

# 如果有.x后缀的列，去掉后缀
names(combinetable) <- gsub("\\.x$", "", names(combinetable))

# 保存结果
write_xlsx(combinetable, "D:/datasets/yunfu/raw_data/switch_results/switch_dependent_followup_new.xlsx")

# 打印结果信息
cat("合并前数据行数:", nrow(switch), "\n")
cat("合并后数据行数:", nrow(combinetable), "\n")