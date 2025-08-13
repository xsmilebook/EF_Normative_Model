library(readxl)
library(dplyr)

# ==== 读取数据 ====
df_gng    <- read_excel("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/dataclean/cleandata/gonogo/gonogo_cleaned.xlsx")
df_2back  <- read_excel("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/dataclean/cleandata/2back/2back_cleaned.xlsx")
df_1back  <- read_excel("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/dataclean/cleandata/1back/1back_cleaned.xlsx")

# ==== 过滤清华大学西区 ====
df_gng    <- df_gng    %>% filter(School != "清华大学西区")
df_2back  <- df_2back  %>% filter(School != "清华大学西区")
df_1back  <- df_1back  %>% filter(School != "清华大学西区")

# ==== 提取 ID ====
id_gng    <- df_gng$ID
id_2back  <- df_2back$ID
id_1back  <- df_1back$ID

# ==== 检查重复 ID ====
dup_gng    <- id_gng[duplicated(id_gng)]
dup_2back  <- id_2back[duplicated(id_2back)]
dup_1back  <- id_1back[duplicated(id_1back)]

cat("gng 中重复 ID 个数:", length(unique(dup_gng)), "\n")
cat("2back 中重复 ID 个数:", length(unique(dup_2back)), "\n")
cat("1back 中重复 ID 个数:", length(unique(dup_1back)), "\n")

# ==== 检查交集 ====
common_ids3  <- Reduce(intersect, list(id_gng, id_2back, id_1back))
common_ids12 <- intersect(id_2back, id_1back)
common_idsG2 <- intersect(id_2back, id_gng)
common_idsG1 <- intersect(id_1back, id_gng)

cat("三个任务都参与的人数:", length(common_ids3), "\n")
cat("1back 和 2back 都参与的人数:", length(common_ids12), "\n")
cat("Go/No-Go 和 2back 都参与的人数:", length(common_idsG2), "\n")
cat("Go/No-Go 和 1back 都参与的人数:", length(common_idsG1), "\n")

# ==== 所有参与者总数 ====
all_ids <- unique(c(id_gng, id_2back, id_1back))
cat("参与总人数（去重后）:", length(all_ids), "\n")

# ==== 各任务男性数量 ====
n_M_gng    <- df_gng   %>% filter(Gender =="男") %>% nrow()
n_M_2back  <- df_2back %>% filter(Gender =="男") %>% nrow()
n_M_1back  <- df_1back %>% filter(Gender =="男") %>% nrow()

cat("Go/No-Go 中男性数量:", n_M_gng, "\n")
cat("2-back 中男性数量:", n_M_2back, "\n")
cat("1-back 中男性数量:", n_M_1back, "\n")

# ==== 各交集中男性数量 ====

n_M_common3 <- df_gng %>% filter(ID %in% common_ids3, Gender =="男") %>% nrow()
n_M_G1 <- df_gng %>% filter(ID %in% common_idsG1, Gender =="男") %>% nrow()
n_M_G2 <- df_gng %>% filter(ID %in% common_idsG2, Gender =="男") %>% nrow()
n_M_12 <- df_2back %>% filter(ID %in% common_ids12, Gender =="男") %>% nrow()

cat("三个任务都参与的男性数量:", n_M_common3, "\n")
cat("Go/No-Go 和 1-back 都参与的男性数量:", n_M_G1, "\n")
cat("Go/No-Go 和 2-back 都参与的男性数量:", n_M_G2, "\n")
cat("1-back 和 2-back 都参与的男性数量:", n_M_12, "\n")

# ==== 辅助函数：计算男性人数 ====
get_M_count <- function(all_ids, df_list){
  all_M_ids <- unlist(
    lapply(df_list, function(df){
      df %>% filter(Gender =="男") %>% pull(ID)
    })
  )
  length(intersect(unique(all_M_ids), all_ids))
}

# ==== 联合 ID 组合 ====
id_1or2back     <- union(id_1back, id_2back)
id_gng_or_1back <- union(id_gng, id_1back)
id_gng_or_2back <- union(id_gng, id_2back)
id_all          <- unique(c(id_gng, id_1back, id_2back))

# ==== 输出结果 ====
cat("\n====== 综合统计（去重）======\n")
cat("1back + 2back 总人数:", length(id_1or2back), "\n")
cat("男性人数:", get_M_count(id_1or2back, list(df_1back, df_2back)), "\n\n")

cat("gng + 1back 总人数:", length(id_gng_or_1back), "\n")
cat("男性人数:", get_M_count(id_gng_or_1back, list(df_gng, df_1back)), "\n\n")

cat("gng + 2back 总人数:", length(id_gng_or_2back), "\n")
cat("男性人数:", get_M_count(id_gng_or_2back, list(df_gng, df_2back)), "\n\n")

cat("所有任务至少参与一项 总人数:", length(id_all), "\n")
cat("男性人数:", get_M_count(id_all, list(df_gng, df_1back, df_2back)), "\n")





################################################Final################################################################################

library(readxl)
library(dplyr)

# ==== 读取数据 ====
df_gng    <- read_excel("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/combinedata/combined_questionnaire/Q_GNG.xlsx")
df_2back  <- read_excel("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/combinedata/combined_questionnaire/Q_2back.xlsx")
df_1back  <- read_excel("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/combinedata/combined_questionnaire/Q_1back.xlsx")


# ==== 过滤清华大学西区 ====
df_gng    <- df_gng    %>% filter(School != "清华大学西区")
df_2back  <- df_2back  %>% filter(School != "清华大学西区")
df_1back  <- df_1back  %>% filter(School != "清华大学西区")

# ==== 提取 ID ====
id_gng    <- df_gng$ID
id_2back  <- df_2back$ID
id_1back  <- df_1back$ID

# ==== 检查重复 ID ====
dup_gng    <- id_gng[duplicated(id_gng)]
dup_2back  <- id_2back[duplicated(id_2back)]
dup_1back  <- id_1back[duplicated(id_1back)]

cat("gng 中重复 ID 个数:", length(unique(dup_gng)), "\n")
cat("2back 中重复 ID 个数:", length(unique(dup_2back)), "\n")
cat("1back 中重复 ID 个数:", length(unique(dup_1back)), "\n")

# ==== 检查交集 ====
common_ids3  <- Reduce(intersect, list(id_gng, id_2back, id_1back))
common_ids12 <- intersect(id_2back, id_1back)
common_idsG2 <- intersect(id_2back, id_gng)
common_idsG1 <- intersect(id_1back, id_gng)

cat("三个任务都参与的人数:", length(common_ids3), "\n")
cat("1back 和 2back 都参与的人数:", length(common_ids12), "\n")
cat("Go/No-Go 和 2back 都参与的人数:", length(common_idsG2), "\n")
cat("Go/No-Go 和 1back 都参与的人数:", length(common_idsG1), "\n")

# ==== 所有参与者总数 ====
all_ids <- unique(c(id_gng, id_2back, id_1back))
cat("参与总人数（去重后）:", length(all_ids), "\n")

# ==== 各任务男性数量 ====
n_M_gng    <- df_gng   %>% filter(Sex == "M") %>% nrow()
n_M_2back  <- df_2back %>% filter(Sex == "M") %>% nrow()
n_M_1back  <- df_1back %>% filter(Sex == "M") %>% nrow()

cat("Go/No-Go 中男性数量:", n_M_gng, "\n")
cat("2-back 中男性数量:", n_M_2back, "\n")
cat("1-back 中男性数量:", n_M_1back, "\n")

# ==== 各交集中男性数量 ====

n_M_common3 <- df_gng %>% filter(ID %in% common_ids3, Sex == "M") %>% nrow()
n_M_G1 <- df_gng %>% filter(ID %in% common_idsG1, Sex == "M") %>% nrow()
n_M_G2 <- df_gng %>% filter(ID %in% common_idsG2, Sex == "M") %>% nrow()
n_M_12 <- df_2back %>% filter(ID %in% common_ids12, Sex == "M") %>% nrow()

cat("三个任务都参与的男性数量:", n_M_common3, "\n")
cat("Go/No-Go 和 1-back 都参与的男性数量:", n_M_G1, "\n")
cat("Go/No-Go 和 2-back 都参与的男性数量:", n_M_G2, "\n")
cat("1-back 和 2-back 都参与的男性数量:", n_M_12, "\n")

# ==== 辅助函数：计算男性人数 ====
get_M_count <- function(all_ids, df_list){
  all_M_ids <- unlist(
    lapply(df_list, function(df){
      df %>% filter(Sex == "M") %>% pull(ID)
    })
  )
  length(intersect(unique(all_M_ids), all_ids))
}

# ==== 联合 ID 组合 ====
id_1or2back     <- union(id_1back, id_2back)
id_gng_or_1back <- union(id_gng, id_1back)
id_gng_or_2back <- union(id_gng, id_2back)
id_all          <- unique(c(id_gng, id_1back, id_2back))

# ==== 输出结果 ====
cat("\n====== 综合统计（去重）======\n")
cat("1back + 2back 总人数:", length(id_1or2back), "\n")
cat("男性人数:", get_M_count(id_1or2back, list(df_1back, df_2back)), "\n\n")

cat("gng + 1back 总人数:", length(id_gng_or_1back), "\n")
cat("男性人数:", get_M_count(id_gng_or_1back, list(df_gng, df_1back)), "\n\n")

cat("gng + 2back 总人数:", length(id_gng_or_2back), "\n")
cat("男性人数:", get_M_count(id_gng_or_2back, list(df_gng, df_2back)), "\n\n")

cat("所有任务至少参与一项 总人数:", length(id_all), "\n")
cat("男性人数:", get_M_count(id_all, list(df_gng, df_1back, df_2back)), "\n")
################################################data of analysis################################################################################

library(readxl)

# 读取两个文件
df_gng <- read_excel("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/combinedata/combined_questionnaire/Q_GNG.xlsx")
df_2back <- read_excel("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/combinedata/combined_questionnaire/Q_2back.xlsx")
df_1back <- read_excel("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/combinedata/combined_questionnaire/Q_1back.xlsx")


# 提取 ID
id_gng    <- df_gng$ID
id_2back  <- df_2back$ID
id_1back  <- df_1back$ID

# 1. 检查 gng 中是否有重复 ID
dup_gng <- id_gng[duplicated(id_gng)]
cat("gng 中重复 ID 个数:", length(unique(dup_gng)), "\n")
print(unique(dup_gng))
# 2. 检查 2back 中是否有重复 ID
dup_2back <- id_2back[duplicated(id_2back)]
cat("2back 中重复 ID 个数:", length(unique(dup_2back)), "\n")
print(unique(dup_2back))

# 3. 检查 1back 中是否有重复 ID
dup_1back <- id_1back[duplicated(id_1back)]
cat("1back 中重复 ID 个数:", length(unique(dup_1back)), "\n")
#print(unique(dup_1back))

# 4. 检查三个 list 中重复 ID（交集）
common_ids3 <- Reduce(intersect, list(id_gng, id_2back, id_1back))
cat("三个文件中重复 ID 个数:", length(common_ids3), "\n")
#print(common_ids3)

common_ids12 <- intersect(id_2back, id_1back)
cat("gng和1back重复 ID 个数:", length(common_ids12), "\n")
#print(common_ids12)

common_idsG2 <- intersect(id_2back, id_gng)
cat("gng和2back ID 个数:", length(common_idsG2), "\n")
#print(common_idsG2)

common_idsG1 <- intersect(id_1back, id_gng)
cat("gng和1back ID 个数:", length(common_idsG1), "\n")
#print(common_idsG1)

# 统计所有唯一的 ID
all_ids <- unique(c(id_gng, id_2back, id_1back))

# 统计每个 ID 出现的次数
all_id_table <- table(c(id_gng, id_2back, id_1back))
cat("参与总人数:", length(all_id_table), "\n")






################################################Plot number of rawdata################################################################################
library(readxl)
library(ggplot2)
library(dplyr)

# 1. 读取数据
df_1back <- read_excel("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/dataclean/cleandata/1back/1back_deleteageconflit.xlsx")
df_2back <- read_excel("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/dataclean/cleandata/2back/2back_deleteageconflit.xlsx")
df_gng <- read_excel("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/dataclean/cleandata/gonogo/gonogo_deleteageconflit.xlsx")


# 2. 预处理：年龄离散化
if (!"Age_year" %in% names(df)) {
  stop("列名 'Age_year' 不存在，请检查拼写")
}
df$Age_year <- floor(df$Age_year)

# 3. 年龄分布柱状图
p_age <- ggplot(df, aes(x = Age_year)) +
  geom_bar(fill = "skyblue", color = "black", width = 0.8) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.3, size = 4, family = "Heiti SC") +
  scale_x_continuous(breaks = seq(min(df$Age_year), max(df$Age_year), by = 1)) +
  labs(title = "1-back年龄分布", x = "Age", y = "Number") +
  theme_minimal(base_size = 14) +
  theme(text = element_text(family = "Heiti SC"))

print(p_age)

library(readxl)
library(ggplot2)
library(dplyr)

# 1. 读取数据
df_1back <- read_excel("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/dataclean/cleandata/1back/1back_deleteageconflit.xlsx")
df_2back <- read_excel("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/dataclean/cleandata/2back/2back_deleteageconflit.xlsx")
df_gng   <- read_excel("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/dataclean/cleandata/gonogo/gonogo_deleteageconflit.xlsx")

# 读取基本信息表
info_df <- read_excel("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/combinedata/combined_questionnaire/YF_information.xlsx") %>%
  select(ID, Sex) %>%
  mutate(Sex = factor(Sex, levels = c("F", "M")))

# 合并 Sex 信息到每个任务数据
df_1back <- left_join(df_1back, info_df, by = "ID")
df_2back <- left_join(df_2back, info_df, by = "ID")
df_gng   <- left_join(df_gng,   info_df, by = "ID")

# 2. 预处理函数
plot_age_distribution_stacked <- function(df, task_name) {
  if (!all(c("Age_year", "Sex") %in% names(df))) {
    stop("列名 'Age_year' 或 'Sex' 不存在，请检查拼写")
  }
  
  df <- df %>%
    mutate(Age_year = floor(Age_year)) %>%
    filter(Age_year > 0 & Sex %in% c("F", "M"))
  
  # 设置Sex因子顺序，F在上，M在下（确保堆叠顺序）
  df$Sex <- factor(df$Sex, levels = c("F", "M"))
  
  # 计算每个年龄的总人数（用于标注）
  count_df <- df %>%
    group_by(Age_year) %>%
    summarise(n = n(), .groups = "drop")
  
  ggplot(df, aes(x = Age_year, fill = Sex)) +
    geom_bar(position = "stack", color = "black", width = 0.8) +
    
    # 添加总人数标签（只一个）
    geom_text(data = count_df, aes(x = Age_year, y = n, label = n),
              vjust = -0.3, size = 2, inherit.aes = FALSE) +
    
    scale_fill_manual(values = c("M" = "#A4C5DF", "F" = "#E5E9F2"), name = "Sex") +
    scale_x_continuous(breaks = seq(min(df$Age_year, na.rm = TRUE), max(df$Age_year, na.rm = TRUE), by = 1), expand = c(0, 0)) +
    labs(title = paste0(task_name, " distribution of age"), x = "Age (years)", y = "Number") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 9, hjust = 0.5, color = "black"),
      axis.text = element_text(colour = "black",size=6), 
      axis.title = element_text(size=6,face = "plain"),
      axis.line = element_line(colour = "black", size = 0.25),
      axis.ticks = element_line(colour = "black", size = 0.25), 
      axis.ticks.length = unit(0.05, "cm"),
      axis.text.y = element_text(colour = "black",size = 6),
      plot.background=element_rect(fill="white",color = NA),
      legend.title=element_text(size=6),
      legend.text=element_text(size=6),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey80"),
      legend.position = "none"
    )
}
p1 <- plot_age_distribution_stacked(df_1back, "1-back")
p2 <- plot_age_distribution_stacked(df_2back, "2-back")
p3 <- plot_age_distribution_stacked(df_gng,   "Go/No-Go")

print(p1)
print(p2)
print(p3)

ggsave("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/article_202505/figure_2507/1back_age_distribution.pdf", plot = p1, width = 18, height = 6, units="cm", dpi = 300)
ggsave("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/article_202505/figure_2507/2back_age_distribution.pdf", plot = p2, width = 18, height = 6, units="cm", dpi = 300)
ggsave("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/article_202505/figure_2507/gonogo_age_distribution.pdf", plot = p3, width = 18, height = 6, units="cm", dpi = 300)



# 4. 筛选18–19岁子集
df_18_19 <- df %>% filter(Age_year %in% c(18, 19))

# 5. 年级分布（18–19岁）
p_grade <- ggplot(df_18_19, aes(x = as.factor(Grade))) +
  geom_bar(fill = "lightgreen", color = "black", width = 0.7) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.3, size = 4, family = "Heiti SC") +
  labs(title = "18–19岁参与者的年级分布", x = "Grade", y = "Number") +
  theme_minimal(base_size = 14) +
  theme(text = element_text(family = "Heiti SC"))

print(p_grade)

# 6. 学校分布（18–19岁）
p_school <- ggplot(df_18_19, aes(x = School)) +
  geom_bar(fill = "salmon", color = "black", width = 0.7) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.3, size = 4, family = "Heiti SC") +
  labs(title = "18–19岁参与者的学校分布", x = "School", y = "Number") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Heiti SC"),
    text = element_text(family = "Heiti SC")
  )

print(p_school)





#GNG
# 1. 读取数据
df <- read_excel("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/dataclean/cleandata/gonogo/gonogo_deleteageconflit.xlsx")

# 2. 年龄离散化
if (!"Age_year" %in% names(df)) {
  stop("列名 'Age_year' 不存在，请检查拼写")
}
df$Age_year <- floor(df$Age_year)
df <- df %>% filter(Age_year >= 0) 

# 3. 年龄分布柱状图
p_age <- ggplot(df, aes(x = Age_year)) +
  geom_bar(fill = "skyblue", color = "black", width = 0.8) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.3, size = 4, family = "Heiti SC") +
  scale_x_continuous(breaks = seq(min(df$Age_year), max(df$Age_year), by = 1)) +
  labs(title = "Go/No-Go 年龄分布", x = "Age", y = "Number") +
  theme_minimal(base_size = 14) +
  theme(text = element_text(family = "Heiti SC"))

print(p_age)

# 4. 筛选 18–19 岁子集
df_18_19 <- df %>% filter(Age_year %in% c(18, 19))

# 5. 年级分布图（18–19 岁）
p_grade <- ggplot(df_18_19, aes(x = as.factor(Grade))) +
  geom_bar(fill = "lightgreen", color = "black", width = 0.7) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.3, size = 4, family = "Heiti SC") +
  labs(title = "18–19岁参与者的年级分布", x = "Grade", y = "Number") +
  theme_minimal(base_size = 14) +
  theme(text = element_text(family = "Heiti SC"))

print(p_grade)

# 6. 学校分布图（18–19 岁）
p_school <- ggplot(df_18_19, aes(x = School)) +
  geom_bar(fill = "salmon", color = "black", width = 0.7) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.3, size = 4, family = "Heiti SC") +
  labs(title = "18–19岁参与者的学校分布", x = "School", y = "Number") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Heiti SC"),
    text = element_text(family = "Heiti SC")
  )

print(p_school)
