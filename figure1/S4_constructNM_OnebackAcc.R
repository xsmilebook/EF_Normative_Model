rm(list=ls())
library(ggplot2)
library(tidyverse)
library(mgcv)
library(readxl)
library(parallel)
library(gamlss)
library(scales)
# input directory
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- '/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psy/data'
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psy/interfileFolder"
  functionFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psy/Rcode/functions"
  resultFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psy/results"
  FigureFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psy/FigureFolder"
}else{
  datapath <- 'D:/datasets/yunfu/raw_data'
  FigureFolder <- 'D:/datasets/yunfu/figures/fig1'
  interfileFolder <- "D:/datasets/yunfu/interfile_folder/Normative_Model"
  functionFolder <- "D:/code/EF_Normative_Model/functions"
  resultFolder <- "D:/datasets/yunfu/results/gamlss/Normative_Model"
}

# load data
source(paste0(functionFolder, "/Construct_gamlss_set_new.R"))
back1_data <- read_xlsx(paste0(datapath, "/Q_1back.xlsx"))

## 1st. 2 fold CV, sampling, set1: 1/2; set2: 1/2, statified by Sex
set.seed(925)
stratify <- "Sex"
stratify.var <- back1_data[ ,stratify]
SPLIT <- split(1:NROW(back1_data), stratify.var)
LAPPLY <- lapply(SPLIT,function(X){sample(x=X,size=length(X)*1/2,replace=F)})
index_set1_b <- unlist(LAPPLY)
back1_data.set1 <- back1_data[index_set1_b,]
back1_data.set2 <- back1_data %>% filter(!(row_number() %in% index_set1_b))#Subtract set1 from back1_data
## 2nd, construct models
back1_data.set1 <- back1_data.set1 %>% select(c("Age_year", "Sex", "School", "ID", "Oneback_acc")) %>% drop_na()
back1_data.set2 <- back1_data.set2 %>% select(c("Age_year", "Sex", "School", "ID", "Oneback_acc")) %>% drop_na()
back1_data.set1$Sex <- as.factor(back1_data.set1$Sex)
back1_data.set2$Sex <- as.factor(back1_data.set2$Sex)
### set1
dataname <- "back1_data.set1"
back1_data.set1 <- as.data.frame(back1_data.set1)
smoothterm <- "Age_year"
dependentvar <- "Oneback_acc"
randomvar <- NA
mu.df <- sigma.df <- 2
degree <- 2
distribution.fam <- "SEP2"
IDvar <- "ID"
covariates <- "Sex"
stratify <- "Sex"
quantile.vec <- c(0.01,0.025, 0.05, 0.25, 0.5, 0.75, 0.95,0.975, 0.99)

if (!file.exists(paste0(interfileFolder, "/1-back/GAMLSS.back1Accset1.sum.rds"))){
  mod.set1 <- construct_gamlss(dataname, dependentvar, smoothterm, covariates, mu.df, sigma.df, degree, distribution.fam, IDvar, quantile.vec, stratify, randomvar=NA)
  saveRDS(mod.set1, paste0(interfileFolder, "/1-back/GAMLSS.back1Accset1.sum.rds"))
}else{
  mod.set1 <- readRDS(paste0(interfileFolder, "/1-back/GAMLSS.back1Accset1.sum.rds"))
}


# performance
modelperformance.set1 <- mod.set1$performance.tb
print(paste(sum(modelperformance.set1$converged), "models converged.")) 
# replicate previous results
#SCrankcorr(modelperformance.set1, "partialRsq", ds.resolution)
# ds.resolution Interest.var r.spearman   p.spearman
# 1            12   partialRsq -0.4724948 1.255304e-05

# compute deviations
gam.data2 <- back1_data.set1
mod.tmp <- mod.set1$mod.tmp

mu_pred <- predict(mod.tmp, newdata = back1_data.set2, what = "mu", type = "response")
sigma_pred <- predict(mod.tmp, newdata = back1_data.set2, what = "sigma", type = "response")


dependentvar <- "Oneback_acc"
deviation.set2.df <- data.frame(ID=back1_data.set2$ID)
observation <- back1_data.set2[[dependentvar]]
centile <- pSEP2(observation, mu = mu_pred, sigma = sigma_pred)
deviation.set2.df[[paste0(dependentvar, "_centile")]] <- centile
deviation.set2.df[[paste0(dependentvar, "_deviationZ")]] <- qnorm(centile)
deviation.set2.df[[paste0(dependentvar, "_sigma")]] <- sigma_pred
deviation.set2.df[[paste0(dependentvar, "_mu")]] <- mu_pred


saveRDS(deviation.set2.df, paste0(interfileFolder, "/1-back/EF_back1Acc.set2_deviation.rds"))


# set2
dataname <- "back1_data.set2"
if (!file.exists(paste0(interfileFolder, "/1-back/GAMLSS.back1Accset2.sum.rds"))){
  mod.set2 <- construct_gamlss(dataname, dependentvar, smoothterm, covariates, mu.df, sigma.df, degree, distribution.fam,IDvar, quantile.vec, stratify,randomvar=NA)
  saveRDS(mod.set2, paste0(interfileFolder, "/1-back/GAMLSS.back1Accset2.sum.rds"))
}else{
  mod.set2 <- readRDS(paste0(interfileFolder, "/1-back/GAMLSS.back1Accset2.sum.rds"))
}

# performance
modelperformance.set2 <- mod.set2$performance.tb
print(paste(sum(modelperformance.set2$converged), "models converged.")) 

# compute deviations
gam.data2 <- back1_data.set2
mod.tmp <- mod.set2$mod.tmp
mu_pred <- predict(mod.tmp, newdata = back1_data.set1, what = "mu", type = "response")
sigma_pred <- predict(mod.tmp, newdata = back1_data.set1, what = "sigma", type = "response")


dependentvar <- "Oneback_acc"
deviation.set1.df <- data.frame(ID=back1_data.set1$ID)
observation <- back1_data.set1[[dependentvar]]
centile <- pSEP2(observation, mu = mu_pred, sigma = sigma_pred)
deviation.set1.df[[paste0(dependentvar, "_centile")]] <- centile
deviation.set1.df[[paste0(dependentvar, "_deviationZ")]] <- qnorm(centile)
deviation.set1.df[[paste0(dependentvar, "_sigma")]] <- sigma_pred
deviation.set1.df[[paste0(dependentvar, "_mu")]] <- mu_pred


saveRDS(deviation.set1.df, paste0(interfileFolder, "/1-back/EF_back1Acc.set1_deviation.rds"))

## 3rd. merge datasets
# The deviations will be averaged across the two models.
back1Acc_deviation.set1 <- deviation.set1.df %>% select(c("ID", ends_with("_centile"), ends_with("_deviationZ"), ends_with("_sigma"), ends_with("_mu"), ends_with("_nu")))
back1Acc_deviation.set2 <- deviation.set2.df %>% select(c("ID", ends_with("_centile"), ends_with("_deviationZ"), ends_with("_sigma"), ends_with("_mu"), ends_with("_nu")))

# Create an empty dataframe for averaged deviations
back1Acc_deviation <- data.frame(ID = back1Acc_deviation.set1$ID)

# Calculate the average centile and deviationZ for the variable
centilename <- "Oneback_acc_centile"
deviationname <- "Oneback_acc_deviationZ"

back1Acc_deviation <- rbind(back1Acc_deviation.set1, back1Acc_deviation.set2) 

# Merge with original data if necessary
back1Acc_deviation <- merge(back1Acc_deviation, back1_data, by="ID")

# Save the result
write.csv(back1_data.set1, paste0(interfileFolder,"/1-back/back1Acc.data1.csv"))
write.csv(back1_data.set2, paste0(interfileFolder,"/1-back/back1Acc.data2.csv"))
saveRDS(back1Acc_deviation, paste0(interfileFolder,"/1-back/back1Acc.deviations.rds"))
write.csv(back1Acc_deviation, paste0(interfileFolder,"/1-back/back1Acc.deviations.csv"))
