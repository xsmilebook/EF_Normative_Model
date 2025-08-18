# correlation between ordinal variables & continuous variables.
# age effects need to be removed from the continous variables.
library(polycor)
library(tidyverse)
library(mgcv)
library(psych)

ordinalcorr <- function(dependentvar, dataname, interest.indep.var, covariates, smoothvar, knots, set_fx = FALSE, stats_only = FALSE){
  
  # 从环境中获取数据框
  gam.data <- get(dataname)
  
  # --- 数据预处理 ---
  # 1. 确保因变量是整数
  gam.data[[dependentvar]] <- as.integer(gam.data[[dependentvar]])
  
  # 过滤掉关键变量为NA的行
  gam.data <- gam.data %>%
    filter(!is.na(.data[[dependentvar]]), !is.na(.data[[interest.indep.var]]))
  
  # 2. 动态计算类别数量 R
  unique_responses <- sort(unique(gam.data[[dependentvar]]))
  R <- length(unique_responses)
  
  # 3. 检查并调整因变量，确保其从1开始
  if (min(unique_responses, na.rm = TRUE) < 1) {
    gam.data[[dependentvar]] <- gam.data[[dependentvar]] - min(unique_responses, na.rm = TRUE) + 1
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
  
  # 在正确的模型拟合后，再计算模型摘要
  gam.results <- summary(gam.model)
  
  # --- 提取统计量 ---
  
  gam.smooth.t <- gam.results$p.table[2, 3]
  gam.smooth.pvalue <- gam.results$p.table[2, 4]
  beta <- gam.results$p.table[2, 1]
  
  anova.cov.pvalue <- anova.gam(gam.model.null, gam.model, test='Chisq')$`Pr(>Chi)`[2]
  
  # --- 相关性分析 ---
  
  res1 <- gam.data[[dependentvar]]
  
  residual_formula <- as.formula(sprintf("%s ~ %s + s(%s, k=%d, fx=F)", interest.indep.var, covariates, smoothvar, knots))
  
  # 【BUG修复】使用 gam() 替代 lm() 来正确处理 s() 平滑项
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