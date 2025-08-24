# 清理工作环境
rm(list=ls())

# 加载必要的R包
library(ggplot2)
library(tidyverse)
library(mgcv)
library(readxl)
library(openxlsx) # 目标代码原有，予以保留
library(parallel)
library(gamlss)
library(scales)

# --- 1. 设置文件路径 ---
# 与目标代码保持一致的路径设置逻辑
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/EFdevelopment/Crsen_survey_GD_YunFu/data'
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/EFdevelopment/Crsen_survey_GD_YunFu/interfileFolder"
  functionFolder <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/EFdevelopment/Crsen_survey_GD_YunFu/Rcode_EFnorms/functions"
  resultFolder <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/EFdevelopment/Crsen_survey_GD_YunFu/results"
}else{
  datapath <- 'D:/datasets/yunfu/raw_data'
  FigureFolder <- 'D:/datasets/yunfu/figures/fig1'
  interfileFolder <- "D:/datasets/yunfu/interfile_folder/Normative_Model"
  functionFolder <- "D:/code/EF_Normative_Model/functions"
  resultFolder <- "D:/datasets/yunfu/results/gamlss/Normative_Model"
}

construct_gamlss <- function(dataname, dependentvar, smoothterm, covariates, mu.df, sigma.df, degree, distribution.fam,IDvar, quantile.vec, stratify=NULL, randomvar=NA){
  
  # get data
  gam.data <- get(dataname)
  covariates <- gsub(" ", "", covariates)
  if (! is.na(randomvar)){
    gam.data1 <- gam.data %>% select(all_of(c(dependentvar, smoothterm, unlist(strsplit(covariates, "+", fixed=T)), randomvar, IDvar))) %>% drop_na()
    gam.data1$randomvar <- as.factor(gam.data1[[randomvar]])
    gam.data2 <- gam.data1 %>% group_by(randomvar) %>%
      filter(n() > 30) %>%
      ungroup()
  }else{
    gam.data2 <- gam.data %>% select(c(dependentvar, smoothterm, unlist(strsplit(covariates, "+", fixed=T)), IDvar)) %>% drop_na()
  }
  
  con<-gamlss.control(n.cyc=200)
  gam.data2 <- as.data.frame(gam.data2)
  # ******** 这一步是必须的，让eval能找到数据 ********
  assign("gam.data2", gam.data2, envir = .GlobalEnv)
  
  # construct model
  if (! is.na(randomvar)){
    mod.mu.formula <- as.formula(sprintf("%s ~ bs(%s, df = %s, degree=%s) + %s + random(%s)", dependentvar, smoothterm, mu.df, degree, covariates, randomvar))
    mod.mu.null.formula <- as.formula(sprintf("%s ~ %s + random(%s)", dependentvar, covariates, randomvar))
  }else{
    mod.mu.formula <- as.formula(sprintf("%s ~ bs(%s, df = %s, degree=%s) + %s", dependentvar, smoothterm, mu.df, degree, covariates))
    mod.mu.null.formula <- as.formula(sprintf("%s ~ %s", dependentvar, covariates))
  }
  
  command <- paste0("mod.tmp <- gamlss(mod.mu.formula, sigma.formula =~ bs(", smoothterm, ", df = ", sigma.df, ", degree = ", degree, ") + ", 
                    covariates, 
                    ", nu.formula = ~1,family=", distribution.fam,", data=gam.data2, control=con)")
  
  command.null <- paste0("mod.null.tmp <- gamlss(mod.mu.null.formula, sigma.formula =~  ", covariates, 
                         ", nu.formula = ~1,family=", distribution.fam,", data=gam.data2, control=con)")
  
  eval(parse(text = command))
  eval(parse(text = command.null))
  
  # performance
  performance.tb <- data.frame(BIC = mod.tmp$sbc, converged = mod.tmp$converged, partialRsq = NA, Rsq = Rsq(mod.tmp))
  
  sse.model <- sum((mod.tmp$y - fitted(mod.tmp, what = "mu"))^2)
  sse.nullmodel <- sum((mod.null.tmp$y - fitted(mod.null.tmp, what = "mu"))^2)
  partialRsq <- (sse.nullmodel - sse.model) / sse.nullmodel
  
  # first derivative
  if(requireNamespace("gamlss.add", quietly = TRUE)){
    PEF <- gamlss.add::getPEF(mod.tmp, term = smoothterm, n.points = 1000, parameter = "mu", type = "response", plot = FALSE)
    x_values <- seq(min(gam.data2[[smoothterm]]), max(gam.data2[[smoothterm]]), length.out = 1000)
    PEF_test <- PEF(x_values, deriv = 1)
    direction <- sum(PEF_test) / abs(sum(PEF_test))
    performance.tb$partialRsq <- partialRsq * direction
  }
  
  # quantiles
  centiles_strat <- list()
  if(requireNamespace("gamlss.add", quietly = TRUE)){
    n_quantiles <- length(quantile.vec)
    n_points <- 1000
    x.tmp <- seq(min(gam.data2[[smoothterm]]), max(gam.data2[[smoothterm]]), length.out=n_points)
    for (l in 1:nlevels(as.factor(gam.data2[[stratify]]))){
      centile.tmp <- array(NA, dim=c(n_quantiles, n_points))
      for (q in 1:n_quantiles){
        fixed_list <- list(levels(as.factor(gam.data2[[stratify]]))[l])
        names(fixed_list) <- stratify
        command <- paste0("Qua <- gamlss.add::getQuantile(mod.tmp, quantile=quantile.vec[q], term = smoothterm, fixed.at = list(", stratify, "='", levels(as.factor(gam.data2[[stratify]]))[l], "'), n.points = 1000)")
        eval(parse(text = command))
        centile.tmp[q, ] <- Qua(x.tmp)
      }
      centiles_strat[[l]] <- centile.tmp
    }
  }
  
  
  sumlist <- list(performance.tb=performance.tb, mod.tmp=mod.tmp, centiles_strat=centiles_strat)
  
  return(sumlist)
}



GNG_data <- read_xlsx(paste0(datapath, "/Q_GNG.xlsx")) # 修改为read_xlsx与参考代码一致

# --- 3. 数据预处理与划分 ---
# 2折交叉验证，按性别（Sex）分层抽样
set.seed(925)
stratify <- "Sex" # 使用 "Sex" 与参考代码保持一致
stratify.var <- GNG_data[ ,stratify]
SPLIT <- split(1:NROW(GNG_data), stratify.var)
LAPPLY <- lapply(SPLIT,function(X){sample(x=X,size=length(X)*1/2,replace=F)})
index_set1_b <- unlist(LAPPLY)
GNG_data.set1 <- GNG_data[index_set1_b,]
# 修改数据划分逻辑，与参考代码一致
GNG_data.set2 <- GNG_data %>% filter(!(row_number() %in% index_set1_b))

# =========================================================================
# === 最终修正版的主脚本（从第4部分开始） ===
# =========================================================================

# =========================================================================
# === 最终修正版的主脚本（采纳您的方案） ===
# =========================================================================

# --- 4. 构建模型 (Set1训练, Set2测试) ---
# 准备数据子集
# 移除了模型不需要的 "School" 列，以确保训练集和测试集结构一致
GNG_data.set1 <- GNG_data.set1 %>% dplyr::select(c("Age_year", "Sex", "ID", "d_prime")) %>% drop_na()
GNG_data.set2 <- GNG_data.set2 %>% dplyr::select(c("Age_year", "Sex", "ID", "d_prime")) %>% drop_na()
GNG_data.set1$Sex <- as.factor(GNG_data.set1$Sex)
GNG_data.set2$Sex <- as.factor(GNG_data.set2$Sex)

# 设置模型参数
dataname <- "GNG_data.set1"
GNG_data.set1 <- as.data.frame(GNG_data.set1)
smoothterm <- "Age_year"
dependentvar <- "d_prime"
randomvar <- NA
mu.df <- sigma.df <- 2
degree <- 2
distribution.fam <- "SEP3"
IDvar <- "ID"
covariates <- "Sex"
quantile.vec <- c(0.01,0.025, 0.05, 0.25, 0.5, 0.75, 0.95,0.975, 0.99)

# 创建新的保存路径
save_subdir <- paste0(interfileFolder, "/GoNoGo_dprime_results")
if (!dir.exists(save_subdir)) {
  dir.create(save_subdir, recursive = TRUE)
}

# 训练模型set1
model_path_set1 <- paste0(save_subdir, "/GNG_model_set1.rds")
if(!file.exists(model_path_set1)){
  # 为了让函数内部能找到数据，在调用前赋值
  gam.data2 <- GNG_data.set1
  mod.set1 <- construct_gamlss(dataname, dependentvar, smoothterm, covariates, mu.df, sigma.df, degree, distribution.fam,IDvar, quantile.vec, stratify, randomvar=NA)
  saveRDS(mod.set1, model_path_set1)
}else{
  mod.set1 <- readRDS(model_path_set1)
}

# 查看模型性能
modelperformance.set1 <- mod.set1$performance.tb
print(paste(sum(modelperformance.set1$converged), "models converged for set 1.")) 

# 在set2上计算离差
# ******** 关键补充：根据您的指示，在predict前定义gam.data2 ********
gam.data2 <- GNG_data.set1 
mod.tmp <- mod.set1$mod.tmp

mu_pred <- predict(mod.tmp, newdata = GNG_data.set2, what = "mu", type = "response")
sigma_pred <- predict(mod.tmp, newdata = GNG_data.set2, what = "sigma", type = "response")
nu_pred <- predict(mod.tmp, newdata = GNG_data.set2, what = "nu", type = "response")
tau_pred <- predict(mod.tmp, newdata = GNG_data.set2, what = "tau", type = "response")

deviation.set2.df <- data.frame(ID=GNG_data.set2$ID)
observation <- GNG_data.set2[[dependentvar]]
centile <- pSEP3(observation, mu = mu_pred, sigma = sigma_pred, nu = nu_pred, tau = tau_pred)
deviation.set2.df[[paste0(dependentvar, "_centile")]] <- centile
deviation.set2.df[[paste0(dependentvar, "_deviationZ")]] <- qnorm(centile)
deviation.set2.df[[paste0(dependentvar, "_sigma")]] <- sigma_pred
deviation.set2.df[[paste0(dependentvar, "_mu")]] <- mu_pred
deviation.set2.df[[paste0(dependentvar, "_nu")]] <- nu_pred
deviation.set2.df[[paste0(dependentvar, "_tau")]] <- tau_pred

saveRDS(deviation.set2.df, paste0(save_subdir, "/GNG_dprime_set2_deviation_scores.rds"))

# --- 5. 构建模型 (Set2训练, Set1测试) ---
dataname <- "GNG_data.set2"
GNG_data.set2 <- as.data.frame(GNG_data.set2)

# 训练模型set2
model_path_set2 <- paste0(save_subdir, "/GNG_model_set2.rds")
if(!file.exists(model_path_set2)){
  # 为了让函数内部能找到数据，在调用前赋值
  gam.data2 <- GNG_data.set2
  mod.set2 <- construct_gamlss(dataname, dependentvar, smoothterm, covariates, mu.df, sigma.df, degree, distribution.fam,IDvar, quantile.vec, stratify, randomvar=NA)
  saveRDS(mod.set2, model_path_set2)
}else{
  mod.set2 <- readRDS(model_path_set2)
}

# 查看模型性能
modelperformance.set2 <- mod.set2$performance.tb
print(paste(sum(modelperformance.set2$converged), "models converged for set 2.")) 

# 在set1上计算离差
# ******** 关键补充：根据您的指示，在predict前定义gam.data2 ********
gam.data2 <- GNG_data.set2
mod.tmp <- mod.set2$mod.tmp

mu_pred <- predict(mod.tmp, newdata = GNG_data.set1, what = "mu", type = "response")
sigma_pred <- predict(mod.tmp, newdata = GNG_data.set1, what = "sigma", type = "response")
nu_pred <- predict(mod.tmp, newdata = GNG_data.set1, what = "nu", type = "response")
tau_pred <- predict(mod.tmp, newdata = GNG_data.set1, what = "tau", type = "response")

deviation.set1.df <- data.frame(ID=GNG_data.set1$ID)
observation <- GNG_data.set1[[dependentvar]]
centile <- pSEP3(observation, mu = mu_pred, sigma = sigma_pred, nu = nu_pred, tau = tau_pred)
deviation.set1.df[[paste0(dependentvar, "_centile")]] <- centile
deviation.set1.df[[paste0(dependentvar, "_deviationZ")]] <- qnorm(centile)
deviation.set1.df[[paste0(dependentvar, "_sigma")]] <- sigma_pred
deviation.set1.df[[paste0(dependentvar, "_mu")]] <- mu_pred
deviation.set1.df[[paste0(dependentvar, "_nu")]] <- nu_pred
deviation.set1.df[[paste0(dependentvar, "_tau")]] <- tau_pred

saveRDS(deviation.set1.df, paste0(save_subdir, "/GNG_dprime_set1_deviation_scores.rds"))

# --- 6. 合并结果并保存 ---
GNGd_prime_deviation.set1 <- deviation.set1.df %>% select(c("ID", ends_with("_centile"), ends_with("_deviationZ"), ends_with("_sigma"), ends_with("_mu"), ends_with("_nu"), ends_with("_tau")))
GNGd_prime_deviation.set2 <- deviation.set2.df %>% select(c("ID", ends_with("_centile"), ends_with("_deviationZ"), ends_with("_sigma"), ends_with("_mu"), ends_with("_nu"), ends_with("_tau")))

GNG_deviation_combined <- rbind(GNGd_prime_deviation.set1, GNGd_prime_deviation.set2) 
GNG_final_deviation_data <- merge(GNG_deviation_combined, GNG_data, by="ID")

# 保存最终结果
write.csv(GNG_data.set1, paste0(save_subdir,"/GNG_dprime_datasubset1.csv"), row.names = FALSE)
write.csv(GNG_data.set2, paste0(save_subdir,"/GNG_dprime_datasubset2.csv"), row.names = FALSE)
saveRDS(GNG_final_deviation_data, paste0(save_subdir,"/GNG_dprime_allsubjects_deviations.rds"))
write.csv(GNG_final_deviation_data, paste0(resultFolder,"/GNG_dprime_allsubjects_deviations.csv"), row.names = FALSE)

print("GoNoGo task processing complete. All files saved.")
# 清理全局环境中可能残留的gam.data2
if(exists("gam.data2")) rm(gam.data2)