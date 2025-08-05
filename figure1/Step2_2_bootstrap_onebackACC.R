rm(list=ls())
#library(ggplot2)
library(tidyverse)
library(mgcv)
library(readxl)
library(parallel)
library(gamlss)
library(scales)
# input directory
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- '/ibmgpfs/cuizaixu_lab/xuhaoshu/datasets/yunfu/raw_data'
  FigureFolder <- '/ibmgpfs/cuizaixu_lab/xuhaoshu/datasets/yunfu/figures/fig1'
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/datasets/yunfu/interfile_folder"
  functionFolder <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/EF_Normative_Model/functions"
  resultFolder <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/datasets/yunfu/results/gamlss/select_parameter"
}else{
  datapath <- 'D:/datasets/yunfu/raw_data'
  FigureFolder <- 'D:/datasets/yunfu/figures/fig1'
  interfileFolder <- "D:/datasets/yunfu/interfile_folder"
  functionFolder <- "D:/code/EF_Normative_Model/functions"
  resultFolder <- "D:/datasets/yunfu/results/gamlss/select_parameter"
  
}

# load data
back1_data <- read_xlsx(paste0(datapath, "/Q_1back.xlsx"))
modelsum_1backAcc_bsdf <- readRDS(paste0(resultFolder, "/modelsum_1backAcc_bsdf.rds"))
performance_1backAcc_bsdf <- read.csv(paste0(resultFolder, "/performance_1backAcc_bsdf.csv"))
model.1backAcc <- modelsum_1backAcc_bsdf[[which.min(performance_1backAcc_bsdf$BIC)]]
# source functions
source(paste0(functionFolder, "/Compare_distributions_gamlss_new.R"))
# get input
n <- commandArgs(trailingOnly = TRUE)
n <- as.numeric(n)
# run bootstrap
back1_data$Sex <- as.factor(back1_data$Sex)
print(n)
Base.Seed <- 925
dataname <- "back1_data"
smoothvar <- "Age_year"
model_obj <- model.1backAcc
stratify <- "Sex"

bootstrap.out <- Boot.Function(n, Base.Seed, dataname,smoothvar, model_obj, stratify)
mod.tmp <- bootstrap.out$mod.tmp
gam.data.subset <- bootstrap.out$gam.data.subset

# estimate quantiles
quantiles<-c(0.01, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.99)
n_quantiles <- length(quantiles)
n_points <- 1000
Centiles_male <- array(NA, dim = c(n_quantiles, n_points))
rownames(Centiles_male) <- paste0("quantiles", quantiles)
Centiles_female <- array(NA, dim = c(n_quantiles, n_points))
rownames(Centiles_female) <- paste0("quantiles", quantiles)
x_female <- seq(min(back1_data[[smoothvar]]), max(back1_data[[smoothvar]]), length.out = n_points)
x_male <- seq(min(back1_data[[smoothvar]]), max(back1_data[[smoothvar]]), length.out = n_points)

gam.data.subset <- as.data.frame(gam.data.subset)
for (i in 1:n_quantiles){
  command <- "Qua <- getQuantile(mod.tmp, quantile=quantiles[i], term = 'age', fixed.at = list(sex='F'), n.points = n_points)"
  eval(parse(text=command))
  Centiles_female[i,] <- Qua(x_female)
  
  command <- "Qua <- getQuantile(mod.tmp, quantile=quantiles[i], term = 'age', fixed.at = list(sex='M'), n.points = n_points)"
  eval(parse(text=command))
  Centiles_male[i,] <- Qua(x_male)
}

Centiles <- (Centiles_female + Centiles_male) /2

eps <- 1e-07 
x.tmp0 <- seq(min(back1_data[[smoothvar]]), max(back1_data[[smoothvar]]), length.out = n_points)
x.tmp1 <- x.tmp0 + eps
P50_function_M <- getQuantile(mod.tmp, quantile = 0.5, term = "age", fixed.at = list(Sex = 'M'), n.points = n_points)
P50_function_F <- getQuantile(mod.tmp, quantile = 0.5, term = "age", fixed.at = list(Sex = 'F'), n.points = n_points)

P50_values_tmp0_M <- P50_function_M(x.tmp0)
P50_values_tmp1_M <- P50_function_M(x.tmp1)
P50_derivative_M <- (P50_values_tmp1_M - P50_values_tmp0_M) / eps

P50_values_tmp0_F <- P50_function_F(x.tmp0)
P50_values_tmp1_F <- P50_function_F(x.tmp1)
P50_derivative_F <- (P50_values_tmp1_F - P50_values_tmp0_F) / eps
P50_derivative <- (P50_derivative_M + P50_derivative_F) / 2

newdata_male0 <- data.frame("age" = x.tmp0, Sex = factor("M", levels = c("M", "F")))
newdata_male1 <- data.frame("age" = x.tmp1, Sex = factor("M", levels = c("M", "F")))
newdata_female0 <- data.frame("age" = x.tmp0, Sex = factor("F", levels = c("M", "F")))
newdata_female1 <- data.frame("age" = x.tmp1, Sex = factor("F", levels = c("M", "F")))

sigma_predictions_male0 <- predict(mod.tmp, newdata = newdata_male0, what = "sigma", type = "response")
sigma_predictions_male1 <- predict(mod.tmp, newdata = newdata_male1, what = "sigma", type = "response")
sigma_predictions_female0 <- predict(mod.tmp, newdata = newdata_female0, what = "sigma", type = "response")
sigma_predictions_female1 <- predict(mod.tmp, newdata = newdata_female1, what = "sigma", type = "response")

sigma_derivative_male <- (sigma_predictions_male1 - sigma_predictions_male0) / eps
sigma_derivative_female <- (sigma_predictions_female1 - sigma_predictions_female0) / eps
sigma_derivative <- (sigma_derivative_male + sigma_derivative_female) / 2

predict_sigma_F <- data.frame("age" = x.tmp0, Sex = 'F')
predict_sigma_M <- data.frame("age" = x.tmp0, Sex = 'M')

sigma_pred_F <- predict(mod.tmp, what = "sigma", newdata = predict_sigma_F, type = "response",data = gam.data.subset)
sigma_pred_M <- predict(mod.tmp, what = "sigma", newdata = predict_sigma_M, type = "response",data = gam.data.subset)
sigma_pred <- (sigma_pred_F + sigma_pred_M) / 2
# mu predict
predict_F <- data.frame("age" = x.tmp0, Sex = 'F')
predict_M <- data.frame("age" = x.tmp0, Sex = 'M')

mu_pred_F <- predict(mod.tmp, what = "mu", newdata = predict_F, type = "response",data = gam.data.subset)
mu_pred_M <- predict(mod.tmp, what = "mu", newdata = predict_M, type = "response",data = gam.data.subset)
mu_pred <- (mu_pred_F + mu_pred_M) / 2

# Extract results
bootstrap_centiles <- list(
  Centiles_female = Centiles_female,
  Centiles_male = Centiles_male,
  Centiles_overall = Centiles,
  Age_points = x.tmp0,
  EFvar = "OnebackAcc",
  bootstrap_time = n,
  converged = mod.tmp$converged,
  mu_pred = mu_pred,
  mu_derivative = bootstrap.out$mu_derivative,
  sigma_derivative = bootstrap.out$sigma_derivative,
  P50_derivative = P50_derivative,
  sigma_derivative_all = sigma_derivative,
  sigma_derivative_male = sigma_derivative_male,
  sigma_derivative_female = sigma_derivative_female,
  sigma_pred = sigma_pred_M,
  sigma_pred_F = sigma_pred_F,
  sigma_pred_M = sigma_pred_M
)
# save out
if (! dir.exists(paste0(interfileFolder, "/bootstrap/1backAcc"))){
  dir.create(paste0(interfileFolder, "/bootstrap/1backAcc"))
}

saveRDS(bootstrap_centiles, paste0(interfileFolder, "/bootstrap/1backAcc/centile_bootstrap_", n, ".rds"))

print(paste("Bootsrap", n, "finished!"))
