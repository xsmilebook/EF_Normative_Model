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
# 函数定义：您提供的用于连续变量分析的新函数
# ===================================================================

# modified anovaPB: return distribution of simulated statistics
anovaPB_ext <- function(objectNull, object, n.sim = 999, 
                        colRef = switch(class(object)[1], "lm" = 5, "lmerMod" = 6, "glmmTMB" = 6, 4),
                        rowRef = 2, ncpus = NULL, ...) {
  
  # check model
  if (length(unlist(coef(objectNull))) > length(unlist(coef(object))))
    stop("The first object (null model) should be smaller than the alternative model.")
  
  respDimnames <- dimnames(model.response(model.frame(object)))
  
  if (is.null(ncpus)) {
    ncpus <- min(70, parallel::detectCores() - 2)
  }
  cat(paste("\nNote: This run will attempt to allocate", ncpus, "cores for permutation testing.\n"))
  
  targs <- match.call(expand.dots = FALSE)
  anovaFn <- anova
  statObs <- try(anova(objectNull, object, ...))
  
  if (inherits(statObs, "try-error")) {
    anovaFn <- function(objectNull, object, ...) {
      llAlt  <- logLik(object)
      llNull <- logLik(objectNull)
      table <- data.frame(
        df = c(attr(llNull, "df"), attr(llAlt, "df")),
        deviance = -2 * c(llNull, llAlt),
        LRT = c(NA, -2 * llNull + 2 * llAlt)
      )
      return(table)
    }
    statObs <- anovaFn(objectNull, object)
    statObs$P <- c(NA, NA)
    names(statObs)[4] <- 'Pr(>LRT)'
    colRef <- 3
    modelnamelist <- c(deparse(substitute(objectNull)), deparse(substitute(object)))
    Xnames <- c(paste(deparse(formula(objectNull), width.cutoff = 500), collapse = "\n"),
                paste(deparse(formula(object), width.cutoff = 500), collapse = "\n"))
    topnote <- paste(modelnamelist, ": ", Xnames, sep = "", collapse = "\n")
    title <- "Analysis of Deviance Table\n"
    rownames(statObs) <- modelnamelist
    attr(statObs, "heading") <- c(title, topnote)
  }
  
  stats <- rep(NA, n.sim + 1)
  stats[1] <- statObs[rowRef, colRef]
  
  if (inherits(object, c("lmerMod", "glmerMod"))) {
    cll <- object@call
    mf <- match.call(call = cll)
    dat <- if (.hasSlot(object, "data")) object@data else NULL
  } else {
    cll <- object$call
    mf <- match.call(call = cll)
    dat <- object$data
  }
  m <- match(c("formula", "data", "subset", "weights", "na.action", "etastart", "mustart", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)
  
  modelF <- try(eval(mf, parent.frame()), silent = TRUE)
  if (inherits(modelF, "try-error") | inherits(object, c("lmerMod", "glmerMod", "glmmTMB")))
    modelF <- model.frame(object)
  
  respName <- names(model.frame(object))[1]
  whichResp <- 1
  if (is.null(dat) == FALSE) {
    if (is.list(dat)) {
      whichAdd <- which(names(dat) %in% names(modelF) == FALSE)
      if(length(whichAdd) > 0) 
        for (iAdd in whichAdd) {
          if (!is.list(dat[[iAdd]])) modelF[[names(dat)[iAdd]]] <- dat[[iAdd]]
        }
    }
  }
  
  offs = NULL
  modelF$offs = try(model.offset(modelF))
  
  if (regexpr("(", respName, fixed = TRUE) > 0) {
    newResp <- sprintf("`%s`", respName)
    fm.update <- reformulate(".", response = newResp)
  } else {
    fm.update <- reformulate(".")
  }
  
  is.mva <- ncol(as.matrix(modelF[[whichResp]])) > 1
  yNew <- simulate(objectNull, n.sim)
  
  getStat <- function(iSim, yNew, objectNull, object, modelF, anovaFn, is.mva, fm.update, whichResp, respDimnames, rowRef, colRef) {
    modelF[[whichResp]] <- if (is.mva) yNew[,,iSim] else as.matrix(yNew[, iSim], dimnames = respDimnames)
    if (inherits(modelF$offs, "try-error") | is.null(modelF$offs)) {
      objectiNull = update(objectNull, formula = fm.update, data = modelF)
      objecti = update(object, formula = fm.update, data = modelF)
    } else {
      objectiNull = update(object, formula = fm.update, data = modelF, offset = offs)
      objecti = update(object, formula = fm.update, data = modelF, offset = offs)
    }
    return(anovaFn(objectiNull, objecti, ...)[rowRef, colRef])
  }
  
  if (ncpus > 1) {
    cl <- NULL
    tryCatch({
      cl <- parallel::makeCluster(ncpus)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(cl, c("getStat", "anovaFn", "yNew", "objectNull", "object", "modelF", "is.mva", "fm.update", "whichResp", "respDimnames", "rowRef", "colRef", as.character(cll[[1]])), envir = environment())
      statList <- parallel::clusterApplyLB(cl, 1:n.sim, getStat, yNew = yNew, objectNull = objectNull, object = object, modelF = modelF, anovaFn = anovaFn, is.mva = is.mva, fm.update = fm.update, whichResp = whichResp, respDimnames = respDimnames, rowRef = rowRef, colRef = colRef)
      stats[-1] <- unlist(statList)
    }, error = function(e) {
      warning("Parallel processing failed, falling back to serial execution: ", e$message)
      for (iBoot in 2:(n.sim + 1)) {
        stats[iBoot] <- getStat(iBoot - 1, yNew = yNew, objectNull = objectNull, object = object, modelF = modelF, anovaFn = anovaFn, is.mva = is.mva, fm.update = fm.update, whichResp = whichResp, respDimnames = respDimnames, rowRef = rowRef, colRef = colRef)
      }
    })
  } else {
    for (iBoot in 2:(n.sim + 1)) {
      stats[iBoot] <- getStat(iBoot - 1, yNew = yNew, objectNull = objectNull, object = object, modelF = modelF, anovaFn = anovaFn, is.mva = is.mva, fm.update = fm.update, whichResp = whichResp, respDimnames = respDimnames, rowRef = rowRef, colRef = colRef)
    }
  }
  
  statReturn <- statObs[, 1:colRef, drop = FALSE]
  statReturn$`P-value` <- NA
  observed_stat <- stats[1]
  simulated_stats <- stats[-1]
  p_value <- mean(simulated_stats >= observed_stat - 1e-8, na.rm = TRUE)
  statReturn$`P-value`[rowRef] <- p_value
  
  hasP <- grep("Pr", colnames(statObs))
  if (length(hasP) > 0) colnames(statReturn)[ncol(statReturn)] <- colnames(statObs)[hasP[1]]
  
  attr(statReturn, "heading") <- attr(statObs, "heading")
  class(statReturn) <- c("anovaPB", class(statObs))
  attr(statReturn, "simulated_stats") <- simulated_stats
  attr(statReturn, "observed_stat")  <- observed_stat
  attr(statReturn, "n.sim")          <- n.sim
  attr(statReturn, "full_stats_vector") <- stats  
  
  return(statReturn)
}

gam.fit.Independent.var <- function(dependentvar, dataname, smoothvar, interest.indep.var, covariates, stats_only = FALSE){  
  gam.data <- get(dataname)
  gam.data <- gam.data %>%
    filter(!is.na(.data[[interest.indep.var]]), !is.na(.data[[dependentvar]]))
  
  Independent.var<-gam.data[[interest.indep.var]]
  indep_outlier_index <- which(Independent.var < mean(Independent.var) - 3 * sd(Independent.var) | Independent.var > mean(Independent.var) + 3 * sd(Independent.var))
  if (length(indep_outlier_index) > 0) { gam.data <- gam.data[-indep_outlier_index, ] }
  
  dependent.var.vec <- gam.data[[dependentvar]]
  dep_outlier_index <- which(dependent.var.vec < mean(dependent.var.vec) - 3 * sd(dependent.var.vec) | dependent.var.vec > mean(dependent.var.vec) + 3 * sd(dependent.var.vec))
  if (length(dep_outlier_index) > 0) { gam.data <- gam.data[-dep_outlier_index, ] }
  
  parcel <- as.character(dependentvar)
  modelformula <- as.formula(sprintf("%s ~ %s + %s + s(%s, k=3, fx=F)",dependentvar, interest.indep.var, covariates, smoothvar))
  modelformula.null<-as.formula(sprintf("%s ~ %s + s(%s, k=3, fx=F)",dependentvar, covariates, smoothvar))
  gam.model <- gam(modelformula, method="REML", data = gam.data)
  gam.results <- summary(gam.model)
  gam.model.null <- gam(modelformula.null, method="REML", data = gam.data)
  
  gam.indep.t <- gam.results$p.table[2,3]
  gam.indep.pvalue <- gam.results$p.table[2,4]
  
  if (stats_only){
    anova_results <- anovaPB_ext(gam.model.null, gam.model, n.sim = 1000, test = 'Chisq')
    anova.pvalues <- anova_results$`Pr(>Chi)`[2]
    simulated_stats <- attr(anova_results, "simulated_stats")
    observed_stat <- attr(anova_results, "observed_stat")
  } else {
    anova.pvalues <- NA
    simulated_stats <- NULL
    observed_stat <- NULL
  }
  
  sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
  sse.nullmodel <- sum((gam.model.null$y - gam.model.null$fitted.values)^2)
  partialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  if(gam.indep.t < 0){ partialRsq <- partialRsq*-1 }
  
  varcorformula1 <- as.formula(sprintf("%s ~ %s+s(%s, k=3, fx=F)", dependentvar, covariates, smoothvar))
  res1<-residuals(gam(varcorformula1, method="REML", data = gam.data))
  res2<-gam.data[[interest.indep.var]]
  res1.normtest <- ks.test(res1, "pnorm", mean = mean(res1), sd = sd(res1))
  res2.normtest <- ks.test(res2, "pnorm", mean = mean(res2), sd = sd(res2))
  
  corrmethod <- ifelse(res1.normtest$p.value > 0.01 & res2.normtest$p.value > 0.01, "pearson", "spearman")
  
  PCorr_Test <- corr.test(res1, res2, method=corrmethod)
  correstimate<-as.numeric(PCorr_Test$r)
  corrp <- as.numeric(PCorr_Test$p)
  
  samplesize <- nrow(gam.data)
  beta = gam.results$p.table[2,1]
  stats.results <- cbind(parcel, interest.indep.var, gam.indep.t, gam.indep.pvalue, anova.pvalues, partialRsq, corrmethod, correstimate, corrp, samplesize, beta)
  
  sim_results <- list(simulated_stats = simulated_stats, observed_stat = observed_stat, n_sim = if (!is.null(simulated_stats)) length(simulated_stats) else NA)
  
  if(stats_only == TRUE)
    return(list(stats = as.data.frame(stats.results), simulation = sim_results))
  else
    return(list(stats = as.data.frame(stats.results), simulation = sim_results, data = data.frame(SCres=res1, cogres=res2)))
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
source(paste0(functionFolder, "/ordinalcorr_new.R")) 

# 读取数据
PHQ_data <- read_xlsx(paste0(interfileFolder, '/PHQ.xlsx'))
GNGd_data <- read_rds(paste0(interfileFolder, '/GNGd_prime.deviations.rds'))
back1_data <- read_rds(paste0(interfileFolder, '/back1Acc.deviations.rds'))
back2_data <- read_rds(paste0(interfileFolder, '/back2Acc.deviations.rds'))

# 数据预处理
PHQ_data <- PHQ_data %>%
  select(用户ID, PHQ_y09, PHQ_sum) %>%  
  mutate(用户ID = as.character(用户ID))

GNGd_data <- GNGd_data %>% mutate(x__ID = as.character(x__ID))
back1_data <- back1_data %>% mutate(x__ID = as.character(x__ID))
back2_data <- back2_data %>% mutate(x__ID = as.character(x__ID))

gngd_dev_col <- names(GNGd_data)[str_ends(names(GNGd_data), "deviationZ")]
back1_dev_col <- names(back1_data)[str_ends(names(back1_data), "deviationZ")]
back2_dev_col <- names(back2_data)[str_ends(names(back2_data), "deviationZ")]

# 创建用于分析的数据框
gngd_analysis_data <- GNGd_data %>% select(x__ID, Age_year, Gender, all_of(gngd_dev_col)) %>% rename(用户ID = x__ID) %>% inner_join(PHQ_data, by = "用户ID")
back1_analysis_data <- back1_data %>% select(x__ID, Age_year, Gender, all_of(back1_dev_col)) %>% rename(用户ID = x__ID) %>% inner_join(PHQ_data, by = "用户ID")
back2_analysis_data <- back2_data %>% select(x__ID, Age_year, Gender, all_of(back2_dev_col)) %>% rename(用户ID = x__ID) %>% inner_join(PHQ_data, by = "用户ID")

# 为 PHQ_y09 (序数) 准备数据
gngd_analysis_data_y09 <- gngd_analysis_data %>% mutate(PHQ_y09 = as.integer(PHQ_y09) + 1) %>% na.omit()
back1_analysis_data_y09 <- back1_analysis_data %>% mutate(PHQ_y09 = as.integer(PHQ_y09) + 1) %>% na.omit()
back2_analysis_data_y09 <- back2_analysis_data %>% mutate(PHQ_y09 = as.integer(PHQ_y09) + 1) %>% na.omit()

# 为 PHQ_sum (连续) 准备数据
gngd_analysis_data_sum <- gngd_analysis_data %>% mutate(PHQ_sum = as.numeric(PHQ_sum)) %>% na.omit()
back1_analysis_data_sum <- back1_analysis_data %>% mutate(PHQ_sum = as.numeric(PHQ_sum)) %>% na.omit()
back2_analysis_data_sum <- back2_analysis_data %>% mutate(PHQ_sum = as.numeric(PHQ_sum)) %>% na.omit()
gngd_analysis_data_sum <- gngd_analysis_data_sum %>%
  mutate(Gender = factor(Gender))

back1_analysis_data_sum <- back1_analysis_data_sum %>%
  mutate(Gender = factor(Gender))
back2_analysis_data_sum <- back2_analysis_data_sum %>%
  mutate(Gender = factor(Gender))
# --- 在这里增加下面的代码 ---

# 对 PHQ_sum 进行 Z-score 标准化
gngd_analysis_data_sum$PHQ_sum_z <- scale(gngd_analysis_data_sum$PHQ_sum)
back1_analysis_data_sum$PHQ_sum_z <- scale(back1_analysis_data_sum$PHQ_sum)
back2_analysis_data_sum$PHQ_sum_z <- scale(back2_analysis_data_sum$PHQ_sum)

# --- 增加代码结束 ---


# === 分析 1: PHQ_y09 (序数变量) ===
cat("\n=== 正在运行 PHQ_y09 (序数变量) 的相关分析... ===\n")
gngd_results_y09 <- ordinalcorr(dependentvar = "PHQ_y09", dataname = "gngd_analysis_data_y09", interest.indep.var = gngd_dev_col, covariates = "Gender", smoothvar = "Age_year", knots = 3, stats_only = TRUE)
back1_results_y09 <- ordinalcorr(dependentvar = "PHQ_y09", dataname = "back1_analysis_data_y09", interest.indep.var = back1_dev_col, covariates = "Gender", smoothvar = "Age_year", knots = 3, stats_only = TRUE)
back2_results_y09 <- ordinalcorr(dependentvar = "PHQ_y09", dataname = "back2_analysis_data_y09", interest.indep.var = back2_dev_col, covariates = "Gender", smoothvar = "Age_year", knots = 3, stats_only = TRUE)
all_results_y09 <- rbind(cbind(as.data.frame(gngd_results_y09), dataset = "GNGd"), cbind(as.data.frame(back1_results_y09), dataset = "back1"), cbind(as.data.frame(back2_results_y09), dataset = "back2"))

# === 分析 2: PHQ_sum (连续变量, 使用新函数) ===
cat("\n=== 正在运行 PHQ_sum (连续变量) 的相关分析... ===\n")
# --- 修改：调用新函数并提取 $stats 部分 ---
gngd_results_sum <- gam.fit.Independent.var(dependentvar = "PHQ_sum_z", dataname = "gngd_analysis_data_sum", interest.indep.var = gngd_dev_col, covariates = "Gender", smoothvar = "Age_year", stats_only = TRUE)$stats
back1_results_sum <- gam.fit.Independent.var(dependentvar = "PHQ_sum_z", dataname = "back1_analysis_data_sum", interest.indep.var = back1_dev_col, covariates = "Gender", smoothvar = "Age_year", stats_only = TRUE)$stats
back2_results_sum <- gam.fit.Independent.var(dependentvar = "PHQ_sum_z", dataname = "back2_analysis_data_sum", interest.indep.var = back2_dev_col, covariates = "Gender", smoothvar = "Age_year", stats_only = TRUE)$stats
all_results_sum <- rbind(cbind(gngd_results_sum, dataset = "GNGd"), cbind(back1_results_sum, dataset = "back1"), cbind(back2_results_sum, dataset = "back2"))

# === 合并并保存所有结果 ===
# --- 修改：在合并前统一列名 ---
all_results_y09 <- all_results_y09 %>%
  rename(t_value = gam.smooth.t, p_value_parametric = gam.smooth.pvalue, p_value_anova = anova.cov.pvalue) %>%
  mutate(corrp = NA) # 添加 corrp 列以匹配

all_results_sum <- all_results_sum %>%
  rename(t_value = gam.indep.t, p_value_parametric = gam.indep.pvalue, p_value_anova = anova.pvalues)

# 合并两个分析的结果
final_results_df <- rbind(all_results_y09, all_results_sum)

# 整理数据类型
final_results_df <- final_results_df %>%
  mutate(
    across(everything(), as.character), # 先全部转为字符避免类型问题
    t_value = as.numeric(t_value),
    p_value_parametric = as.numeric(p_value_parametric),
    p_value_anova = as.numeric(p_value_anova),
    partialRsq = as.numeric(partialRsq),
    correstimate = as.numeric(correstimate),
    corrp = as.numeric(corrp),
    samplesize = as.numeric(samplesize),
    beta = as.numeric(beta)
  )

# 查看结果
print(final_results_df)

# 保存结果
write.xlsx(final_results_df, paste0(resultFolder, "/PHQ_all_deviation_correlations_new.xlsx"))

# 打印摘要信息
cat("\n=== 所有相关分析结果摘要 ===\n")
for(i in 1:nrow(final_results_df)) {
  cat(sprintf("因变量: %s, 数据集: %s\n", final_results_df$parcel[i], final_results_df$dataset[i]))
  cat(sprintf("自变量: %s\n", final_results_df$interest.indep.var[i]))
  cat(sprintf("相关系数 (方法: %s): %.3f (p=%.4f)\n", final_results_df$corrmethod[i], final_results_df$correstimate[i], final_results_df$corrp[i]))
  cat(sprintf("模型 ANOVA p值 (置换检验): %.4f\n", final_results_df$p_value_anova[i]))
  cat(sprintf("偏R²: %.3f\n", final_results_df$partialRsq[i]))
  cat("---\n")
}