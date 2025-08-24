library(dplyr)
library(ggplot2)

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


df_all <- df_all %>% filter(!is.na(PHQ_y09))

df_all <- df_all %>%
  mutate(Age_Group = ntile(Age_year, 3)) %>%
  mutate(Age_Group = factor(Age_Group, labels = c("stage: 1", "stage: 2", "stage: 3")))

# 按年龄组和性别统计人数
count_result <- df_all %>%
  group_by(Age_Group, PHQ_y09) %>%
  summarise(Count = n(), .groups = "drop")

# 查看结果
print(count_result)


figure <- ggplot(count_result, aes(x = Age_Group, y = Count, fill = PHQ_y09)) +
  geom_col() +
  labs(title = "Number of Participants by Age Group and PHQ_y09",
       x = "Age Group", y = "Count", fill = "PHQ_y09") +
  theme_minimal()


print(figure)


