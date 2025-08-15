# 清理环境
rm(list=ls())

# 加载必要的包
library(ggplot2)
library(tidyverse)
library(mgcv)
library(readxl)
library(parallel)
library(gamlss)
library(scales)
library(tableone)
library(openxlsx)
library(patchwork)
library(ggsignif) # 用于添加显著性标记

# --- 路径设置 (与原代码相同) ---
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- '/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/EF_Yunfu/interfileFolder_back12before'
  datapath_GNG <- '/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_all/interfileFolder'
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/EF_Yunfu/interfileFolder_back12before"
  functionFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/EF_Yunfu/code/Rcode/functions"
  resultFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/EF_Yunfu/results/corr"
}else{
  datapath <- 'D:/datasets/yunfu/raw_data'
  FigureFolder <- 'D:/datasets/yunfu/figures/fig1'
  interfileFolder <- "D:/datasets/yunfu/interfile_folder/temp"
  functionFolder <- "D:/code/EF_Normative_Model/functions"
  resultFolder <- "D:/datasets/yunfu/results/psy_corr"
}

# --- 数据导入与合并 (与原代码相同) ---
# source functions
#source(paste0(functionFolder, "/gam_varyingcoefficients.R"))
source(paste0(functionFolder, "/gamcog_withsmoothvar_deviation.R"))
#source(paste0(functionFolder, "/ordinalcorr_new.R"))

# import dataset
# PHQ_data <- read_xlsx(paste0(interfileFolder, '/PHQ.xlsx'))
GNGd_data <- read_rds(paste0(interfileFolder, '/Q_GNG_SBQR.rds'))
back1_data <- read_rds(paste0(interfileFolder, '/Q_1back_SBQR.rds'))
back2_data <- read_rds(paste0(interfileFolder, '/Q_2back_SBQR.rds'))

GNGd_data <- GNGd_data %>%
  filter(GNGd_data$SBQ_y01R != -1)
  
back1_data <- back1_data %>%
  filter(back1_data$SBQ_y01R != -1)

back2_data <- back2_data %>%
  filter(back2_data$SBQ_y01R != -1)


GNGd_data$Sex_factor <- as.factor(GNGd_data$Sex)
back1_data$Sex_factor <- as.factor(back1_data$Sex)
back2_data$Sex_factor <- as.factor(back2_data$Sex)

plot_data <- bind_rows(
  GNGd_data %>% 
    select(ID, SBQ_y01R, Age_year, Sex_factor, d_prime_deviationZ) %>%
    rename(acc_deviationZ = d_prime_deviationZ) %>%
    mutate(task = "GNGd"),
  
  back1_data %>% 
    select(ID, SBQ_y01R, Age_year, Sex_factor, Oneback_acc_deviationZ) %>%
    rename(acc_deviationZ = Oneback_acc_deviationZ) %>%
    mutate(task = "Back1"),
  
  back2_data %>% 
    select(ID, SBQ_y01R, Age_year, Sex_factor, Twoback_acc_deviationZ) %>%
    rename(acc_deviationZ = Twoback_acc_deviationZ) %>%
    mutate(task = "Back2")
) %>%
  filter(!is.na(SBQ_y01R))

plot_data$SBQ_y01R <- as.factor(plot_data$SBQ_y01R)

cat("--- 开始执行 ANOVA + TukeyHSD 统计分析 ---\n")

# 创建用于存储结果的空列表
ancova_posthoc_results <- list()
tasks_to_analyze <- unique(plot_data$task)

for (current_task in tasks_to_analyze) {
  cat(paste0("\n--- 正在分析任务: ", current_task, " ---\n"))
  data_clean <- plot_data %>% filter(task == current_task)
  data_clean$SBQ_y01R <- factor(data_clean$SBQ_y01R)
  
  # 拟合 ANOVA 模型
  ancova_model <- aov(acc_deviationZ ~ SBQ_y01R, data = data_clean)
  p_value <- anova(ancova_model)$`Pr(>F)`[1]
  
  cat(paste0("ANOVA 主效应 P 值: ", round(p_value, 5), "\n"))
  
  # 只有当 ANOVA 主效应显著时才进行事后分析
  if (!is.na(p_value) && p_value < 0.05) {
    cat("主效应显著，开始执行 Tukey's HSD 事后检验...\n")
    tukey_result <- TukeyHSD(ancova_model, "SBQ_y01R")
    tukey_df <- as.data.frame(tukey_result$SBQ_y01R)
    tukey_df$comparison <- rownames(tukey_df)
    tukey_df$Task <- current_task
    ancova_posthoc_results[[current_task]] <- tukey_df
  } else {
    cat("ANOVA 主效应不显著 (p >= 0.05)，跳过事后检验。\n")
  }
  
  # --- 新增功能：计算 Spearman 相关性 ---
  cat("---\n计算 Spearman 相关性:\n")
  # cor.test 需要数值输入，因此我们将因子转为数值
  cor_result <- cor.test(~ acc_deviationZ + as.numeric(as.character(SBQ_y01R)), 
                         data = data_clean, 
                         method = "spearman",
                         exact = FALSE) # exact=FALSE 避免因数据有重复值而产生警告
  
  # 提取并打印结果
  rho <- cor_result$estimate
  cor_p_value <- cor_result$p.value
  cat(paste0("Spearman's rho: ", round(rho, 3), "\n"))
  cat(paste0("P-value: ", round(cor_p_value, 5), "\n"))
}

# 组合所有事后分析结果
if (length(ancova_posthoc_results) > 0) {
  all_posthoc_results_df <- do.call(rbind, ancova_posthoc_results)
  rownames(all_posthoc_results_df) <- NULL
} else {
  all_posthoc_results_df <- data.frame() 
}


# --- 4. 转换统计结果以适配绘图函数 ---
cat("\n--- 正在格式化统计结果用于绘图 ---\n")

if (nrow(all_posthoc_results_df) > 0) {
  plot_signif_data <- all_posthoc_results_df %>%
    # 筛选出显著的结果
    filter(`p adj` < 0.05) %>%
    # 将 comparison 列 (如 "1-0") 拆分为 Group1 和 Group2
    separate(comparison, into = c("Group2", "Group1"), sep = "-") %>%
    # 重命名 p.adj 为 P_value，以匹配绘图函数
    rename(P_value = `p adj`) %>%
    # 添加一个与旧代码兼容的 Significant 列
    mutate(Significant = "Yes")
  
  print("以下显著性差异将被绘制在图上:")
  print(plot_signif_data)
} else {
  plot_signif_data <- data.frame() # 如果没有显著结果，创建一个空数据框
}


# --- 5. 最终版绘图函数 (手动控制，无需修改) ---
plot_task_with_signif <- function(data, task_name, significance_df, jitter_width = 0.15, scatter_offset = 0.5) {
  
  task_signif <- significance_df %>% 
    filter(Task == task_name, Significant == "Yes")
  
  task_data <- data %>% filter(task == task_name)
  p <- ggplot(task_data, aes(x = SBQ_y01R, y = acc_deviationZ)) +
    geom_boxplot(color = "black", outlier.shape = NA, alpha = 0.6, width = 0.6) +
    geom_jitter(aes(x = as.numeric(SBQ_y01R) + scatter_offset),
                color = "lightblue", alpha = 0.6, size = 1, 
                width = jitter_width, height = 0) +
    theme_minimal() +
    labs(title = task_name, x = "SBQ_y01R", y = "Deviation Z-score") +
    scale_x_discrete(expand = expansion(add = c(0.8, 0.8))) +
    coord_cartesian(ylim = c(-3, 3.6), clip = "off") 
  
  if (nrow(task_signif) > 0) {
    signif_manual_df <- task_signif %>%
      mutate(
        xmin = as.numeric(factor(Group1, levels = levels(task_data$SBQ_y01R))),
        xmax = as.numeric(factor(Group2, levels = levels(task_data$SBQ_y01R))),
        annotations = sapply(P_value, function(p) {
          if (p < 0.001) " "
          else if (p < 0.01) " "
          else if (p < 0.05) " "
          else " "
        }),
        y_position = seq(3.0, 3.3 + (nrow(.)-1)*0.3, length.out = nrow(.))
      )
    
    p <- p + geom_signif(
      data = signif_manual_df,
      aes(xmin = xmin, xmax = xmax, annotations = annotations, y_position = y_position),
      manual = TRUE,
      tip_length = 0.01,
      textsize = 5,
      vjust = 0.1
    )
  }
  
  return(p)
}


# --- 6. 生成并显示最终图形 ---
cat("\n--- 正在生成最终图形 ---\n")

p1 <- plot_task_with_signif(plot_data, "GNGd", plot_signif_data)
p2 <- plot_task_with_signif(plot_data, "Back1", plot_signif_data)
p3 <- plot_task_with_signif(plot_data, "Back2", plot_signif_data)

combined_plot <- p1 + p2 + p3 +
  plot_layout(ncol = 3) +
  plot_annotation(title = "Deviation Z-scores across SBQ scores for different tasks (ANOVA & TukeyHSD)")

print(combined_plot)