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
library(polycor) # ordinalcorr 函数需要这个包

# --- 这是我们最终修正完善后的 ordinalcorr 函数 ---
ordinalcorr <- function(dependentvar, dataname, interest.indep.var, covariates, smoothvar, knots, set_fx = FALSE, stats_only = FALSE){
  
  # 从环境中获取数据框
  gam.data <- get(dataname)
  
  # --- [诊断信息 1] 函数入口 ---
  cat("\n--- [DIAGNOSTIC INFO] Entering ordinalcorr ---\n")
  cat("Data source (dataname):", dataname, "\n")
  cat("Dimensions of incoming data:", nrow(gam.data), "rows,", ncol(gam.data), "cols\n")
  cat("Frequency table of dependent var '", dependentvar, "' BEFORE internal filtering:\n")
  print(table(gam.data[[dependentvar]], useNA = "ifany"))
  cat("---------------------------------------------\n")
  
  # --- 数据预处理 ---
  
  # 过滤掉关键变量为NA的行
  gam.data <- gam.data %>%
    filter(!is.na(.data[[dependentvar]]), !is.na(.data[[interest.indep.var]]))
  
  # 【关键修复】将因变量转换为从1开始的整数编码，以满足 ocat() 的严格要求。
  # 1. as.factor() 会根据值的顺序创建水平 (levels)。
  # 2. as.integer() 会返回这些水平的整数索引 (从1开始)。
  # 这是将任何一组有序数值/字符转换为 1:k 整数编码的最稳健方法。
  gam.data[[dependentvar]] <- as.integer(as.factor(gam.data[[dependentvar]]))
  
  # --- [诊断信息 2] 内部过滤和转换后 ---
  cat("\n--- [DIAGNOSTIC INFO] After internal filtering and transformation ---\n")
  cat("Dimensions of data after processing:", nrow(gam.data), "rows,", ncol(gam.data), "cols\n")
  cat("Frequency table of dependent var '", dependentvar, "' AFTER transformation (now should be 1, 2, ...):\n")
  print(table(gam.data[[dependentvar]], useNA = "ifany"))
  cat("----------------------------------------------------------------------\n")
  
  # 动态计算类别数量 R
  R <- length(unique(gam.data[[dependentvar]]))
  
  # 安全检查，如果因变量没有变异性，则返回NULL
  if (R < 2) {
    warning(paste("Dependent variable", dependentvar, "has no variance in the provided data subset. Skipping this analysis."))
    return(NULL)
  }
  
  # --- GAM 模型分析 ---
  
  parcel <- as.character(dependentvar)
  
  # 构建模型公式
  modelformula <- as.formula(sprintf("%s ~ %s + %s + s(%s, k=%d, fx=F)", dependentvar, interest.indep.var, covariates, smoothvar, knots))
  modelformula.null <- as.formula(sprintf("%s ~ %s + s(%s, k=%d, fx=F)", dependentvar, covariates, smoothvar, knots))
  
  # 定义正确的 family 参数
  gam.family <- ocat(R = R)
  
  # 拟合模型
  gam.model <- gam(modelformula, method="REML", data = gam.data, family = gam.family)
  gam.model.null <- gam(modelformula.null, method="REML", data = gam.data, family = gam.family)
  
  # 计算模型摘要
  gam.results <- summary(gam.model)
  
  # --- 提取统计量 ---
  
  gam.smooth.t <- gam.results$p.table[2, 3]
  gam.smooth.pvalue <- gam.results$p.table[2, 4]
  beta <- gam.results$p.table[2, 1]
  
  anova.cov.pvalue <- anova.gam(gam.model.null, gam.model, test='Chisq')$`Pr(>Chi)`[2]
  
  # --- 相关性分析 ---
  
  res1 <- gam.data[[dependentvar]]
  
  residual_formula <- as.formula(sprintf("%s ~ %s + s(%s, k=%d, fx=F)", interest.indep.var, covariates, smoothvar, knots))
  
  residual_model <- gam(residual_formula, data = gam.data)
  res2 <- residuals(residual_model)
  
  pcorr_test_result <- tryCatch({
    polycor::polyserial(res2, res1, ML=TRUE, std.err=TRUE)
  }, warning = function(w) { NULL }, error = function(e) { NULL })
  
  if (!is.null(pcorr_test_result)) {
    correstimate <- pcorr_test_result$rho
  } else {
    correstimate <- NA
  }
  
  corrmethod <- "polyserial"
  samplesize <- nrow(gam.data)
  
  stats.results <- cbind(parcel, interest.indep.var, gam.smooth.t, gam.smooth.pvalue, anova.cov.pvalue, corrmethod, correstimate, samplesize, beta)
  
  # --- 返回结果 ---
  
  if (stats_only) {
    return(stats.results)
  } else {
    data.results <- list()
    data.results[[1]] <- as.data.frame(stats.results)
    data.results[[2]] <- data.frame(Dependent.var = res1, Independent.res = res2)
    return(data.results)
  }
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

# 读取数据
SBQ_data <- read_xlsx(paste0(interfileFolder, '/SBQ_R_with_xID.xlsx'))
GNGd_data <- read_rds(paste0(interfileFolder, '/GNGd_prime.deviations.rds'))
back1_data <- read_rds(paste0(interfileFolder, '/back1Acc.deviations.rds'))
back2_data <- read_rds(paste0(interfileFolder, '/back2Acc.deviations.rds'))


# --- 数据预处理 ---
# 从SBQ数据中选择我们需要的列，并重命名ID列以方便合并
SBQ_data_selected <- SBQ_data %>%
  rename(x__ID = "用户ID") %>%
  select(x__ID, SBQ_y02, SBQ_y03a) %>%  
  mutate(x__ID = as.character(x__ID))

# 统一行为数据ID格式
GNGd_data <- GNGd_data %>% mutate(x__ID = as.character(x__ID))
back1_data <- back1_data %>% mutate(x__ID = as.character(x__ID))
back2_data <- back2_data %>% mutate(x__ID = as.character(x__ID))

# 提取每个任务的自变量列名
gngd_dev_col <- names(GNGd_data)[str_ends(names(GNGd_data), "deviationZ")]
back1_dev_col <- names(back1_data)[str_ends(names(back1_data), "deviationZ")]
back2_dev_col <- names(back2_data)[str_ends(names(back2_data), "deviationZ")]

# 创建用于分析的基础数据框
gngd_analysis_data <- GNGd_data %>% 
  select(x__ID, Age_year, Gender, all_of(gngd_dev_col)) %>% 
  inner_join(SBQ_data_selected, by = "x__ID")

back1_analysis_data <- back1_data %>% 
  select(x__ID, Age_year, Gender, all_of(back1_dev_col)) %>% 
  inner_join(SBQ_data_selected, by = "x__ID")

back2_analysis_data <- back2_data %>% 
  select(x__ID, Age_year, Gender, all_of(back2_dev_col)) %>% 
  inner_join(SBQ_data_selected, by = "x__ID")


# --- 对整个样本进行序数变量分析 (不划分年龄段) ---

# 1. 定义要分析的因变量和任务列表
dependent_vars <- c("SBQ_y02", "SBQ_y03a")

tasks <- list(
  GNGd = list(data = gngd_analysis_data, col = gngd_dev_col),
  back1 = list(data = back1_analysis_data, col = back1_dev_col),
  back2 = list(data = back2_analysis_data, col = back2_dev_col)
)

# 2. 初始化一个列表来存储所有结果
all_results_list <- list()

# 3. 循环遍历每个因变量和每个任务
cat("\n=== 正在对整个样本运行 SBQ (序数变量) 的相关分析... ===\n")

for (current_dep_var in dependent_vars) {
  cat(sprintf("\n*** 正在处理因变量: %s ***\n", current_dep_var))
  
  for (task_name in names(tasks)) {
    
    cat(sprintf("--- 任务: %s ---\n", task_name))
    
    dev_col <- tasks[[task_name]]$col
    
    # 定义当前分析所必需的列
    required_cols <- c(current_dep_var, dev_col, "Gender", "Age_year")
    
    # 准备当前分析所需的数据
    current_analysis_data <- tasks[[task_name]]$data %>%
      # 使用 drop_na() 只移除必需列中包含NA的行
      drop_na(all_of(required_cols)) %>%
      mutate(!!sym(current_dep_var) := as.integer(!!sym(current_dep_var)) + 1)
    
    if (nrow(current_analysis_data) < 10) {
      cat("样本量过小，跳过分析。\n")
      next
    }
    
    # 调用函数进行分析
    result <- ordinalcorr(
      dependentvar = current_dep_var, 
      dataname = "current_analysis_data",
      interest.indep.var = dev_col, 
      covariates = "Gender", 
      smoothvar = "Age_year", 
      knots = 3, 
      stats_only = TRUE
    )
    
    if (!is.null(result)) {
      result_df <- as.data.frame(result)
      result_df$dataset <- task_name
      all_results_list[[length(all_results_list) + 1]] <- result_df
    }
  }
}

# 4. 合并所有结果到一个数据框
final_results_df <- do.call(rbind, all_results_list)

# 5. 整理并保存
if (!is.null(final_results_df) && nrow(final_results_df) > 0) {
  final_results_df <- final_results_df %>%
    rename(t_value = gam.smooth.t, p_value_parametric = gam.smooth.pvalue, p_value_anova = anova.cov.pvalue)
  
  print(final_results_df)
  
  output_filename <- paste0(resultFolder, "/SBQ_total_sample_correlations.xlsx")
  write.xlsx(final_results_df, output_filename)
  
  cat("\n=== 分析结果已保存到:", output_filename, "===\n")
} else {
  cat("\n=== 分析未产生任何结果 ===\n")
}


# ===================================================================
# 柱状图可视化
# ===================================================================
if (!is.null(final_results_df) && nrow(final_results_df) > 0) {
  
  # --- 1. 准备用于绘制图表的数据 ---
  plot_data_raw <- final_results_df %>%
    mutate(
      correstimate = as.numeric(as.character(correstimate)),
      p_value_anova = as.numeric(as.character(p_value_anova))
    )
  
  plot_data <- plot_data_raw %>%
    mutate(
      task_factor = factor(dataset, 
                           levels = c("GNGd", "back1", "back2"),
                           labels = c("Go/No-Go", "1-back", "2-back")),
      
      is_significant = p_value_anova < 0.05,
      
      # 创建显著性标签 '*'
      sig_label = ifelse(is_significant, "*", ""),
      
      # 计算标签的Y轴位置，使其总是在柱子顶端外侧
      label_y_position = ifelse(correstimate >= 0, correstimate + 0.005, correstimate - 0.005)
    )
  
  # --- 2. 绘制分面柱状图 ---
  Fig_barplot_corr <- ggplot(plot_data, aes(x = task_factor, y = correstimate, fill = task_factor)) +
    geom_col() +
    
    # 添加 y=0 的水平参考线
    geom_hline(yintercept = 0, color = "gray40") +
    
    # 添加显著性标记 '*'
    geom_text(aes(y = label_y_position, label = sig_label), 
              vjust = ifelse(plot_data$correstimate >= 0, 0, 1), # 根据正负调整垂直对齐
              size = 8, color="black") +
    
    # 使用分面来为 SBQ_y02 和 SBQ_y03a 创建子图
    facet_wrap(~ parcel, ncol = 2) +
    
    # 手动设置填充颜色
    scale_fill_manual(values = c("Go/No-Go" = "#dbe6eb", "1-back" = "#a1c1d5", "2-back" = "#3a7ca5")) +
    
    # 设置图表标题和坐标轴标签
    labs(title = "Correlations of EF and SBQ Items", 
         x = "Executive Function Task", y = "Correlation (ρ)") +
    
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 12, face="bold"),
      axis.text.y = element_text(size = 12, face="bold"),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, hjust = 0.5, face="bold"),
      legend.position = "none", # 隐藏图例
      strip.text = element_text(size = 14, face = "bold"), 
      strip.background = element_rect(fill="grey90", color = NA)
    )
  
  # 打印图表
  print(Fig_barplot_corr)
  
  # (可选) 保存图表
  ggsave(paste0(FigureFolder, "/SBQ_corr_barplot.pdf"), plot = Fig_barplot_corr, width = 10, height = 6)
  
  cat("\n=== 柱状图已生成并保存 ===\n")
}
