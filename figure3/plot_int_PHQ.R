# 清空环境（可选）
# rm(list = ls())

# --- 0. 加载所有必要的 R 包 ---
library(ggplot2)
library(tidyverse)
library(mgcv)
library(readxl)
library(ecostats) # <<<--- 重要：添加此包以使用 anovaPB
library(openxlsx)
# library(parallel) # 如果不用并行计算可以注释掉
# library(gamlss)
# library(scales)
# library(tableone)
# library(showtext)

# --- 1. 文件路径设置 ---
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  # ... (您的路径)
} else {
  datapath <- 'D:/datasets/yunfu/raw_data'
  interfileFolder <- "D:/datasets/yunfu/interfile_folder/temp"
  functionFolder <- "D:/code/EF_Normative_Model/functions"
  resultFolder <- "D:/datasets/yunfu/results/phq_sum_corr"
  FigureFolder <- 'D:/datasets/yunfu/figures/phq_sum_fig'
}
dir.create(resultFolder, showWarnings = FALSE)
dir.create(FigureFolder, showWarnings = FALSE)

# --- 2. 数据读取与预处理 ---

# 读取行为任务数据
GNGd_data <- readRDS(paste0(interfileFolder, '/GNGd_prime.deviations.rds'))
back1_data <- readRDS(paste0(interfileFolder, '/back1Acc.deviations.rds'))
back2_data <- readRDS(paste0(interfileFolder, '/back2Acc.deviations.rds'))

# 读取 PHQ 数据
PHQ_data <- read_xlsx(paste0(interfileFolder, '/PHQ.xlsx')) %>%
  select(用户ID, PHQ_sum) %>%  
  mutate(用户ID = as.character(用户ID))

# 数据类型统一
GNGd_data <- GNGd_data %>% mutate(x__ID = as.character(x__ID), Sex = as.factor(Sex))
back1_data <- back1_data %>% mutate(x__ID = as.character(x__ID), Sex = as.factor(Sex))
back2_data <- back2_data %>% mutate(x__ID = as.character(x__ID), Sex = as.factor(Sex))

# 合并 PHQ 数据
GNGd_data <- GNGd_data %>% inner_join(PHQ_data, by = c("x__ID" = "用户ID"))
back1_data <- back1_data %>% inner_join(PHQ_data, by = c("x__ID" = "用户ID"))
back2_data <- back2_data %>% inner_join(PHQ_data, by = c("x__ID" = "用户ID"))

# (推荐) 检查合并并重命名ID列
if(nrow(GNGd_data) == 0) stop("ERROR: GNGd_data is empty after merging. Check ID columns.")



# 加载自定义函数
# 确保你提供的函数保存在这个路径
source(paste0(functionFolder, '/gam_varyingcoefficients.R')) 

# --- 3. 定义变量与参数 ---
psyc_variables_continous <- c("PHQ_sum") 

EFvars.set <- as.data.frame(matrix(c(
  "d_prime_deviationZ", "GNGd",
  "Oneback_acc_deviationZ", "back1",
  "Twoback_acc_deviationZ", "back2"
), byrow = TRUE, ncol = 2, dimnames = list(NULL, c("varname", "dataname"))))

# 对因变量进行 Z-score 标准化
GNGd_data[, "PHQ_sum_z"] <- scale(GNGd_data[, "PHQ_sum"])
back1_data[, "PHQ_sum_z"] <- scale(back1_data[, "PHQ_sum"])
back2_data[, "PHQ_sum_z"] <- scale(back2_data[, "PHQ_sum"])

psyc_variables_continous_z <- "PHQ_sum_z"

# 设置模型参数
knots <- 3
set_fx <- FALSE 
increments <- 1000
draws <- 1000
return_posterior_coefficients <- TRUE

# --- 4. 初始化结果存储 ---
int.results.df <- data.frame()
significant_interactions <- data.frame()
nosignificant_interactions <- data.frame()

# --- 5. 主循环：执行 GAM 分析 ---
for (i in 1:nrow(EFvars.set)) {
  int_var <- EFvars.set$varname[i]
  dataname_str <- EFvars.set$dataname[i]
  dataname <- paste0(dataname_str, "_data")
  
  for (j in 1:length(psyc_variables_continous_z)) {
    dependentvar <- psyc_variables_continous_z[j]
    smooth_var <- "Age_year"
    covariates <- "Sex"
    
    # 调用 GAM 变系数模型函数
    result <- gam.varyingcoefficients(
      dependentvar = dependentvar, 
      dataname = dataname, 
      smooth_var = smooth_var, 
      int_var = int_var, 
      covariates = covariates,
      knots = knots, 
      set_fx = set_fx, 
      increments = increments, 
      draws = draws, 
      return_posterior_coefficients = return_posterior_coefficients
    )
    
    # 检查函数是否返回了有效结果
    if (!is.null(result)) {
      # 提取结果 (与函数内部的 full.results 对应)
      stats_results <- as.data.frame(result[[1]])
      slope_data <- result[[2]]
      
      # 整理并存储结果
      temp_df <- data.frame(
        dependentvar = dependentvar,
        smooth_var = smooth_var,
        int_var = int_var,
        dataname = dataname_str,
        # *** 关键修正：确保列名与函数输出完全匹配 ***
        IntpartialRsq = as.numeric(stats_results$IntpartialRsq),
        anova.pvalue = as.numeric(stats_results$anova.int.pvalue),
        boots.pvalue = as.numeric(stats_results$anova.pvalues) # anova.pvalues 是自助法的p值
      )
      
      int.results.df <- rbind(int.results.df, temp_df)
      
      # 根据自助法p值(boots.pvalue)判断显著性
      if (temp_df$boots.pvalue < 0.05) {
        significant_interactions <- rbind(significant_interactions, 
                                          data.frame(dependentvar = dependentvar, 
                                                     smooth_var = smooth_var,
                                                     int_var = int_var, 
                                                     dataname = dataname_str, 
                                                     Interaction = temp_df$IntpartialRsq,
                                                     slope_data = I(list(slope_data))))
      } else {
        nosignificant_interactions <- rbind(nosignificant_interactions, 
                                            data.frame(dependentvar = dependentvar, 
                                                       smooth_var = smooth_var,
                                                       int_var = int_var, 
                                                       dataname = dataname_str, 
                                                       Interaction = temp_df$IntpartialRsq,
                                                       slope_data = I(list(slope_data))))
      }
    }
  }
}

# --- 6. 保存与可视化 (这部分代码无需修改，逻辑正确) ---

# 保存完整的统计结果
if(nrow(int.results.df) > 0) {
  write.xlsx(int.results.df, file = paste0(resultFolder, "/PHQ_sum_interaction_results.xlsx"), row.names = FALSE)
}
if(nrow(significant_interactions) > 0) {
  write_rds(significant_interactions, file = paste0(resultFolder, "/PHQ_sum_int_significant.rds"))
}
if(nrow(nosignificant_interactions) > 0) {
  write_rds(nosignificant_interactions, file = paste0(resultFolder, "/PHQ_sum_int_nosignificant.rds"))
}

# --- 6.1 热图可视化 ---
if(nrow(int.results.df) > 0) {
  lwth <- min(int.results.df$IntpartialRsq, na.rm = TRUE)
  upth <- max(int.results.df$IntpartialRsq, na.rm = TRUE)
  
  plot_data <- int.results.df
  plot_data$sig <- plot_data$boots.pvalue < 0.05
  plot_data$int_var <- factor(plot_data$int_var, 
                              levels = c("d_prime_deviationZ", "Oneback_acc_deviationZ", "Twoback_acc_deviationZ"), 
                              labels = c("Go/No-go", "1-back", "2-back"))
  
  custom_colors <- c("#f5dfdb", "#edb8b0", "#e69191", "#c25759")
  
  Fig_heatmap <- ggplot() +
    geom_tile(data = plot_data, aes(x = int_var, y = dependentvar, fill = IntpartialRsq), color = "white") +  
    geom_text(data = subset(plot_data, sig == TRUE), 
              aes(x = int_var, y = dependentvar, label = "*"), 
              vjust = 0.8, hjust = 0.5, size = 9) +
    scale_fill_gradientn(colors = custom_colors, limits = c(lwth, upth), name = "Partial R²") +
    scale_y_discrete(labels = c("PHQ_sum_z" = "PHQ Total Score")) +
    labs(title = "Interaction of Executive Functions and Age on Depression Symptoms", 
         x = "Executive Function Tasks", y = "Mental Health") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title = element_text(size = 16),
      plot.title = element_text(size = 18, hjust = 0.5),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      panel.grid.major = element_blank()
    )
  
  print(Fig_heatmap)
  ggsave(paste0(FigureFolder, "/PHQ_sum_heatmap.pdf"), plot = Fig_heatmap, width = 18, height = 8, units = "cm")
}



# --- 6.2 斜率变化图 (Post-hoc anlysis) ---
# 定义任务和变量的映射，用于美化标题
task_mapping <- c("GNGd" = "Go/No-Go", "back1" = "1-back", "back2" = "2-back")
variable_mapping <- c("PHQ_sum_z" = "PHQ Total Score")

# 绘制显著交互的斜率图
if (nrow(significant_interactions) > 0) {
  for (k in 1:nrow(significant_interactions)) {
    # 提取当前结果的数据
    current_slope_data <- significant_interactions$slope_data[[k]]
    current_task <- significant_interactions$dataname[k]
    current_variable <- significant_interactions$dependentvar[k]
    
    modified_task_name <- task_mapping[current_task]
    modified_variable_name <- variable_mapping[current_variable]
    
    # 计算中位数和95%置信区间
    slope_summary <- current_slope_data %>%
      group_by(Age_year) %>%
      summarise(
        median_slope = median(slope, na.rm = TRUE),
        lower_95CI = quantile(slope, probs = 0.025, na.rm = TRUE),
        upper_95CI = quantile(slope, probs = 0.975, na.rm = TRUE)
      )
    
    # 绘图
    slope_plot <- ggplot(slope_summary, aes(x = Age_year)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_line(aes(y = median_slope), color = "black", size = 0.8) + 
      geom_ribbon(aes(ymin = lower_95CI, ymax = upper_95CI), alpha = 0.2, fill = "blue") + 
      labs(title = paste0(modified_variable_name, " ~ ", modified_task_name), 
           x = "Age (years)", 
           y = "Slope (Effect of EF)") +
      scale_x_continuous(limits = c(11, 18), breaks = seq(11, 18, by = 1)) +
      theme_classic() +
      theme(
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        plot.title = element_text(size = 12, hjust = 0.5)
      )
    
    print(slope_plot)
    ggsave(paste0(FigureFolder, "/significant_slope_", current_task, "_", current_variable, ".pdf"), plot = slope_plot, width = 12, height = 10, units = "cm")
  }
}

# --- 7. 复现图片风格的可视化 ---

# 首先，确保 ggplot2 和一个新的包 patchwork 已经加载
# patchwork 用于轻松地将两张图合并
if (!require("patchwork")) install.packages("patchwork")
library(patchwork)

# --- 准备绘图所需的数据 ---
# (假设 int.results.df 和 significant_interactions 已经存在于你的环境中)

# --- PART 1: 绘制条形图 (复现图 A) ---

# 准备条形图的数据
plot_data_A <- int.results.df %>%
  mutate(
    # 创建显著性标签，p<0.05则为'*'，否则为空
    sig_label = ifelse(boots.pvalue < 0.05, "*", ""),
    # 将任务变量转换为因子，以确保绘图顺序正确
    int_var_factor = factor(int_var, 
                            levels = c("d_prime_deviationZ", "Oneback_acc_deviationZ", "Twoback_acc_deviationZ"),
                            labels = c("Go/No-Go", "1-back", "2-back"))
  )

# 创建图 A
plot_A <- ggplot(plot_data_A, aes(x = int_var_factor, y = IntpartialRsq, fill = int_var_factor)) +
  geom_col() + # 使用 geom_col() 直接绘制数值
  # 在条形图上方添加显著性标记 '*'
  geom_text(aes(label = sig_label), vjust = -0.5, size = 8, color = "black") +
  # 手动设置填充颜色，模仿示例图的由浅到深
  scale_fill_manual(values = c("Go/No-Go" = "#dbe6eb", "1-back" = "#a1c1d5", "2-back" = "#3a7ca5")) +
  # 设置Y轴标签为 β (beta)，并调整范围
  scale_y_continuous(expression(beta), limits = c(0, 0.055), expand = c(0, 0)) +
  # 设置X轴标签
  labs(x = "EF task", y = expression(beta)) +
  # 使用经典主题，并进行精细调整
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none", # 隐藏图例
    axis.line = element_line(size = 0.5),
    axis.ticks = element_line(size = 0.5)
  )

# 打印图 A 查看效果
print(plot_A)


# --- PART 2: 绘制斜率图 (复现图 B) ---

# 从显著结果中选择一个来绘制，例如第一个显著的结果
# 如果没有显著结果，你可以从 nosignificant_interactions 中选择一个
if (nrow(significant_interactions) > 0) {
  plot_data_B_raw <- significant_interactions[1, ]
} else {
  # 如果没有显著的，就用第一个不显著的来演示
  print("没有发现显著的交互作用，将使用一个不显著的结果来绘制示例图B。")
  plot_data_B_raw <- nosignificant_interactions[1, ]
}

# 提取并汇总斜率数据
slope_summary_B <- plot_data_B_raw$slope_data[[1]] %>%
  group_by(Age_year) %>%
  summarise(
    median_slope = median(slope, na.rm = TRUE),
    lower_95CI = quantile(slope, probs = 0.025, na.rm = TRUE),
    upper_95CI = quantile(slope, probs = 0.975, na.rm = TRUE)
  )

# 创建标题 (例如 "PHQ_sum_z ~ 1-back")
plot_title_B <- paste("PHQ Total Score ~", task_mapping[plot_data_B_raw$dataname])


# --- PART 2: 绘制斜率图 (修正 Y 轴格式) ---

# (假设 slope_summary_B, plot_title_B 等变量已经存在)

# *** 关键修正：创建一个自定义的标签格式化函数 ***
format_slope_labels <- function(y) {
  # 使用 sprintf 进行精确格式化
  # ifelse(条件, 如果为真, 如果为假)
  # 如果 y < 0，则正常显示两位小数 (例如 -0.04)
  # 否则 (y >= 0)，在前面加上 '+' 号 (例如 +0.04, +0.00)
  ifelse(y < 0, sprintf("%.2f", y), sprintf("+%.2f", y))
}


# --- 创建图 B (应用新的标签函数) ---
plot_B <- ggplot(slope_summary_B, aes(x = Age_year)) +
  # 添加 y=0 的红色虚线参考线
  geom_hline(yintercept = 0, linetype = "dashed", color = "#d62828", size = 1) +
  # 绘制95%置信区间的灰色带
  geom_ribbon(aes(ymin = lower_95CI, ymax = upper_95CI), alpha = 0.2, fill = "gray50") +
  # 绘制斜率中位数变化的黑色实线
  geom_line(aes(y = median_slope), color = "black", size = 1.2) +
  # 设置坐标轴标签和标题
  labs(
    title = plot_title_B,
    x = "Age (years)",
    y = "Slope"
  ) +
  # 调整X轴
  scale_x_continuous(limits = c(11, 18), breaks = seq(11, 18, by = 2), expand = c(0.01, 0.01)) +
  
  # *** 应用修正后的自定义函数到 Y 轴 ***
  scale_y_continuous(labels = format_slope_labels) +
  
  # 使用经典主题，并进行精细调整
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.line = element_line(size = 0.5),
    axis.ticks = element_line(size = 0.5)
  )

# 打印图 B 查看效果
print(plot_B)


# --- PART 3: 将两张图合并并保存 ---

# 使用 patchwork 包的 '+' 操作符来合并图形
final_plot <- plot_A + plot_B + 
  plot_annotation(tag_levels = 'A') # 自动添加 'A', 'B' 标签

# 打印最终合并的图
print(final_plot)

# 保存最终的图片
ggsave(
  filename = paste0(FigureFolder, "/Combined_Interaction_Plot.pdf"),
  plot = final_plot,
  width = 10, # 宽度 (英寸)
  height = 4, # 高度 (英寸)
  dpi = 300
)



# 加载 ggplot2 和 dplyr 包
library(ggplot2)
library(dplyr)

# --- 1. 根据您的图片手动创建数据框 ---
# 我们从图片中提取需要绘图的核心数据
results_df <- tibble(
  parcel = c("PHQ_y09", "PHQ_y09", "PHQ_y09", "PHQ_sum", "PHQ_sum", "PHQ_sum"),
  interest.indep.var = c("d_prime_deviationZ", "Oneback_acc_deviationZ", "Twoback_acc_deviationZ", 
                         "d_prime_deviationZ", "Oneback_acc_deviationZ", "Twoback_acc_deviationZ"),
  anova.cov.pvalue = c(0.3606, 0.6478, 0.8831, 0.3802, 0.1318, 0.3620),
  beta = c(-0.01654, 0.00835, -0.00340, -0.03446, 0.06020, 0.04566)
)


# --- 2. 准备用于绘图的数据 ---
plot_data <- results_df %>%
  mutate(
    # 将任务变量转换为有序的因子，确保绘图顺序和标签正确
    task_label = factor(interest.indep.var, 
                        levels = c("d_prime_deviationZ", "Oneback_acc_deviationZ", "Twoback_acc_deviationZ"),
                        labels = c("Go/No-Go", "1-back", "2-back")),
    
    # 创建显著性标签：如果p值小于0.05，则为'*'，否则为空字符串
    sig_label = ifelse(anova.cov.pvalue < 0.05, "*", ""),
    
    # 为了让标签 '*' 显示在正确的位置（条形图顶端），我们需要计算其y轴坐标
    # 对于负数的beta值，标签应该在0的稍下方
    label_y_position = ifelse(beta >= 0, beta + 0.005, -0.005)
  )

# --- 3. 使用 ggplot2 绘制柱状图 ---
beta_barplot <- ggplot(plot_data, aes(x = task_label, y = beta, fill = task_label)) +
  # 绘制柱状图
  geom_col(position = "dodge") +
  
  # 添加 y=0 的水平参考线，以清晰地区分正负效应
  geom_hline(yintercept = 0, color = "gray40") +
  
  # 在柱状图的顶端添加显著性标记 '*'
  geom_text(aes(y = label_y_position, label = sig_label), 
            vjust = 0.5, size = 8, color = "black") +
  
  # 使用 facet_wrap 为 PHQ_y09 和 PHQ_sum 创建两个独立的子图
  facet_wrap(~ parcel, scales = "free_y") +
  
  # 手动设置填充颜色
  scale_fill_manual(values = c("Go/No-Go" = "#dbe6eb", "1-back" = "#a1c1d5", "2-back" = "#3a7ca5")) +
  
  # 设置图表标签
  labs(
    title = "Beta Coefficients of EF Tasks on PHQ Scores",
    x = "EF Task",
    y = expression(beta) # 将y轴标签显示为希腊字母 β
  ) +
  
  # 使用一个简洁的主题并进行微调
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"), # 子图标题的样式
    strip.background = element_blank(), # 移除子图标题的背景
    legend.position = "none" # 隐藏图例
  )

# 打印图表
print(beta_barplot)

# (可选) 保存图表到文件2
ggsave(file.path(FigureFolder,"beta_coefficients_barplot.pdf"), plot = beta_barplot, width = 8, height = 5)




# --- PART 2: (最终修改版) 循环生成所有任务的斜率图，并使用指定的 Beta 值作为参考 ---
c




