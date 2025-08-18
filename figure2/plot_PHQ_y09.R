# 清理环境
rm(list=ls())

# 加载所需包
library(readxl)
library(tidyverse)
library(mgcv)
library(psych)
library(openxlsx)
# 新函数所需的包
library(ecostats)
library(pbkrtest)

# ===================================================================
# 函数定义：此处省略您提供的两个大型函数 anovaPB_ext 和 gam.fit.Independent.var
# 在本次分析中，这两个函数不会被用到，但我们保留它们以防万一
# ===================================================================

# modified anovaPB: return distribution of simulated statistics
anovaPB_ext <- function(objectNull, object, n.sim = 999, 
                        colRef = switch(class(object)[1], "lm" = 5, "lmerMod" = 6, "glmmTMB" = 6, 4),
                        rowRef = 2, ncpus = NULL, ...) {
  # ... 函数体和您提供的一样 ...
}

gam.fit.Independent.var <- function(dependentvar, dataname, smoothvar, interest.indep.var, covariates, stats_only = FALSE){  
  # ... 函数体和您提供的一样 ...
}


# ===================================================================
# 主要分析流程
# ===================================================================

# 设置路径
# ... (您的路径设置代码保持不变) ...
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- '/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psy/interfileFolder'
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psy/interfileFolder_back12before"
  functionFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psycode/Rcode/functions"
  resultFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psy/results_corr"
  FigureFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psy/results_corr"
} else {
  datapath <- 'D:/datasets/yunfu/raw_data'
  FigureFolder <- 'D:/datasets/yunfu/figures/fig1'
  interfileFolder <- "D:/datasets/yunfu/interfile_folder/temp"
  functionFolder <- "D:/code/EF_Normative_Model/functions"
  resultFolder <- "D:/datasets/yunfu/results/psy_corr"
}

# 确保您已经加载了更新后的 ordinalcorr 函数
# 注意：请确保 ordinalcorr_new.R 文件存在于您的 functionFolder 路径下
source(paste0(functionFolder, "/ordinalcorr_new.R")) 

# 读取数据
PHQ_data <- read_xlsx(paste0(interfileFolder, '/PHQ.xlsx'))
GNGd_data <- read_rds(paste0(interfileFolder, '/GNGd_prime.deviations.rds'))
back1_data <- read_rds(paste0(interfileFolder, '/back1Acc.deviations.rds'))
back2_data <- read_rds(paste0(interfileFolder, '/back2Acc.deviations.rds'))

# 数据预处理 (只保留 PHQ_y09)
PHQ_data <- PHQ_data %>%
  select(用户ID, PHQ_y09) %>%  
  mutate(用户ID = as.character(用户ID))

GNGd_data <- GNGd_data %>% mutate(x__ID = as.character(x__ID))
back1_data <- back1_data %>% mutate(x__ID = as.character(x__ID))
back2_data <- back2_data %>% mutate(x__ID = as.character(x__ID))

gngd_dev_col <- names(GNGd_data)[str_ends(names(GNGd_data), "deviationZ")]
back1_dev_col <- names(back1_data)[str_ends(names(back1_data), "deviationZ")]
back2_dev_col <- names(back2_data)[str_ends(names(back2_data), "deviationZ")]

# 创建用于分析的基础数据框
gngd_analysis_data <- GNGd_data %>% select(x__ID, Age_year, Gender, all_of(gngd_dev_col)) %>% rename(用户ID = x__ID) %>% inner_join(PHQ_data, by = "用户ID")
back1_analysis_data <- back1_data %>% select(x__ID, Age_year, Gender, all_of(back1_dev_col)) %>% rename(用户ID = x__ID) %>% inner_join(PHQ_data, by = "用户ID")
back2_analysis_data <- back2_data %>% select(x__ID, Age_year, Gender, all_of(back2_dev_col)) %>% rename(用户ID = x__ID) %>% inner_join(PHQ_data, by = "用户ID")

# === 按年龄段对 PHQ_y09 进行序数变量分析 ===

# 1. 定义任务列表和年龄段
tasks <- list(
  GNGd = list(data = gngd_analysis_data, col = gngd_dev_col),
  back1 = list(data = back1_analysis_data, col = back1_dev_col),
  back2 = list(data = back2_analysis_data, col = back2_dev_col)
)

age_groups <- list(
  "Early_Adolescence_11-13" = c(11, 13),
  "Mid_Adolescence_14-15" = c(14, 15),
  "Late_Adolescence_16-18" = c(16, 18)
)

# 2. 初始化一个列表来存储所有结果
all_results_list <- list()

# 3. 循环遍历每个任务和每个年龄段
cat("\n=== 正在按年龄段运行 PHQ_y09 (序数变量) 的相关分析... ===\n")

for (task_name in names(tasks)) {
  for (age_group_name in names(age_groups)) {
    
    cat(sprintf("--- 正在处理任务: %s, 年龄段: %s ---\n", task_name, age_group_name))
    
    # 获取当前任务的数据和年龄范围
    base_data <- tasks[[task_name]]$data
    dev_col <- tasks[[task_name]]$col
    age_range <- age_groups[[age_group_name]]
    
    # 按年龄段筛选数据，并为 ordinalcorr 准备
    current_analysis_data <- base_data %>%
      filter(Age_year >= age_range[1] & Age_year <= age_range[2]) %>%
      mutate(PHQ_y09 = as.integer(PHQ_y09) + 1) %>%
      na.omit()
    
    # 检查筛选后的数据是否为空
    if (nrow(current_analysis_data) < 10) { # 如果样本太少，则跳过
      cat("样本量过小，跳过分析。\n")
      next
    }
    
    # 执行 ordinalcorr 分析
    # 注意：我们将数据框直接传递给函数，而不是它的名字
    result <- ordinalcorr(
      dependentvar = "PHQ_y09", 
      dataname = "current_analysis_data", # 函数内部会使用 get() 获取这个数据框
      interest.indep.var = dev_col, 
      covariates = "Gender", 
      smoothvar = "Age_year", 
      knots = 3, 
      stats_only = TRUE
    )
    
    # 将结果添加到列表中，并附上任务和年龄段信息
    if (!is.null(result)) {
      result_df <- as.data.frame(result)
      result_df$dataset <- task_name
      result_df$age_group <- age_group_name
      all_results_list[[length(all_results_list) + 1]] <- result_df
    }
  }
}

# 4. 合并所有结果到一个数据框
final_results_df <- do.call(rbind, all_results_list)

# 5. 整理数据框
if (!is.null(final_results_df)) {
  # 重命名列以保持一致性
  final_results_df <- final_results_df %>%
    rename(t_value = gam.smooth.t, p_value_parametric = gam.smooth.pvalue, p_value_anova = anova.cov.pvalue)
  
  # 查看结果
  print(final_results_df)
  
  # 保存结果
  write.xlsx(final_results_df, paste0(resultFolder, "/PHQ_y09_by_age_group_correlations.xlsx"))
  
  cat("\n=== 所有分析已完成并保存到:", paste0(resultFolder, "/PHQ_y09_by_age_group_correlations.xlsx"), "===\n")
} else {
  cat("\n=== 分析未产生任何结果 ===\n")
}





# 假设 final_results_df 已经存在于您的 R 环境中

# --- 1. 准备用于绘制热图的数据 ---

# 【BUG修复】在数据准备的第一步，将必要的列转换为数值型
heatmap_data_raw <- final_results_df %>%
  mutate(
    correstimate = as.numeric(as.character(correstimate)),
    p_value_anova = as.numeric(as.character(p_value_anova))
  )

# 现在使用修正后的数据框进行后续操作
heatmap_data <- heatmap_data_raw %>%
  mutate(
    # 将任务变量(dataset)转换为有序因子，确保X轴顺序正确
    task_factor = factor(dataset, 
                         levels = c("GNGd", "back1", "back2"),
                         labels = c("Go/No-Go", "1-back", "2-back")),
    
    # 将年龄段变量(age_group)也转换为有序因子，确保Y轴顺序正确
    age_group_factor = factor(age_group,
                              levels = c("Late_Adolescence_16-18", 
                                         "Mid_Adolescence_14-15", 
                                         "Early_Adolescence_11-13"),
                              labels = c("16-18 years", "14-15 years", "11-13 years")),
    
    # 创建一个逻辑变量来标记结果是否显著
    is_significant = p_value_anova < 0.05
  )

# --- 2. 确定颜色填充的范围 ---
#min_corr <- min(heatmap_data$correstimate, na.rm = TRUE)
# max_corr <- max(heatmap_data$correstimate, na.rm = TRUE)
limit <- max(abs(c(-0.1, 0.1)), na.rm = TRUE)


# --- 3. 绘制热图 ---
Fig_heatmap_corr <- ggplot() +
  # 【警告修复】将 size 修改为 linewidth
  geom_tile(data = heatmap_data, aes(x = task_factor, y = age_group_factor, fill = correstimate), color = "white", linewidth=1.5) +  
  
  # 为显著的结果添加星号标记
  geom_text(data = filter(heatmap_data, is_significant == TRUE), 
            aes(x = task_factor, y = age_group_factor, label = "*"), 
            vjust = 0.8, hjust = 0.5, size = 12, color="black") +
  
  # 设置一个蓝-白-红的发散型颜色条
  scale_fill_gradient2(low = "#3B9AB2", mid = "white", high = "#E15759", 
                       midpoint = 0, limits = c(-limit, limit), name = "Correlation") +
  
  # 设置图表标题和坐标轴标签
  labs(title = "Correlations of EF and PHQ-y09 Across Age Groups", 
       x = "Executive Function Task", y = "Age Group") +
  
  # 使用简洁主题并进行微调
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, face="bold"),
    axis.text.y = element_text(size = 14, face="bold"),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, hjust = 0.5, face="bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.grid.major = element_blank()
  )

# 打印热图
print(Fig_heatmap_corr)

# (可选) 保存热图
ggsave(paste0(FigureFolder, "/PHQ_y09_corr_heatmap.pdf"), plot = Fig_heatmap_corr, width = 10, height = 6)





