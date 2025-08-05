rm(list=ls())
#library(ggplot2)
library(tidyverse)
library(mgcv)
library(readxl)
library(parallel)
library(gamlss)
library(scales)
library(readxl)
# input directory
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- '/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psy/data'
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psy/interfileFolder"
  functionFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psy/Rcode/functions"
  resultFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psy/results"
}else{
  datapath <- 'D:/datasets/yunfu/raw_data'
  FigureFolder <- 'D:/datasets/yunfu/figures/fig1'
  interfileFolder <- "D:/datasets/yunfu/interfile_folder"
  functionFolder <- "D:/code/EF_Normative_Model/functions"
  resultFolder <- "D:/datasets/yunfu/results/gamlss/select_parameter"

}

# source functions
source(paste0(functionFolder, "/Compare_distributions_gamlss_new.R"))

GNG_data <- read_xlsx(paste0(datapath, "/Q_GNG.xlsx"))
back1_data <- read_xlsx(paste0(datapath, "/Q_1back.xlsx"))
back2_data <- read_xlsx(paste0(datapath, "/Q_2back.xlsx"))



# ##### Step1. select the best distribution
# # 1. Go/no-go dprime
# dataname <- "GNG_data"
# GNG_data$Sex <- as.factor(GNG_data$Sex)
# GNG_data$School <- as.factor(GNG_data$School)
# dependentvar <- "d_prime"
# smoothvar <- "Age_year"
# IDvar <- "ID"
# bs.df <- 3
# covariates <- "Sex"
# randomvar=NA
# 
# out_GNGd_prime_distribution <- gamlss_comparedistribution(dataname, dependentvar,
#                                                           smoothvar, IDvar, bs.df,
#                                                           covariates, randomvar=NA)
# names(out_GNGd_prime_distribution) <- c("modelsum", "performance")
# performance_GNGd_prime_distribution <- out_GNGd_prime_distribution$performance
# modelsum_GNGd_prime_distribution <- out_GNGd_prime_distribution$modelsum
# write.csv(performance_GNGd_prime_distribution, paste0(resultFolder, "/performance_GNGd_prime_distribution.csv"), row.names = F)
# saveRDS(modelsum_GNGd_prime_distribution, paste0(resultFolder, "/modelsum_GNGd_prime_distribution.rds"))
# 


# ##### Step1. select the best distribution
# # 2. 1back
# dataname <- "back1_data"
# back1_data$Sex <- as.factor(back1_data$Sex)
# back1_data$School <- as.factor(back1_data$School)
# dependentvar <- "Oneback_acc"
# smoothvar <- "Age_year"
# IDvar <- "ID"
# bs.df <- 3
# covariates <- "Sex"
# randomvar=NA
# 
# out_1backAcc_distribution <- gamlss_comparedistribution(dataname, dependentvar,
#                                                            smoothvar, IDvar, bs.df,
#                                                            covariates, randomvar=NA)
# names(out_1backAcc_distribution) <- c("modelsum", "performance")
# performance_1backAcc_distribution <- out_1backAcc_distribution$performance
# modelsum_1backAcc_distribution <- out_1backAcc_distribution$modelsum
# write.csv(performance_1backAcc_distribution, paste0(resultFolder, "/performance_1backAcc_distribution.csv"), row.names = F)
# saveRDS(modelsum_1backAcc_distribution, paste0(resultFolder, "/modelsum_1backAcc_distribution.rds"))



##### Step1. select the best distribution
# 3. 2back
dataname <- "back2_data"
back2_data$Sex <- as.factor(back2_data$Sex)
back2_data$School <- as.factor(back2_data$School)
dependentvar <- "Twoback_acc"
smoothvar <- "Age_year"
IDvar <- "ID"
bs.df <- 3
covariates <- "Sex"
randomvar=NA

out_2backAcc_distribution <- gamlss_comparedistribution(dataname, dependentvar,
                                                           smoothvar, IDvar, bs.df,
                                                           covariates, randomvar=NA)
names(out_2backAcc_distribution) <- c("modelsum", "performance")
performance_2backAcc_distribution <- out_2backAcc_distribution$performance
modelsum_2backAcc_distribution <- out_2backAcc_distribution$modelsum
write.csv(performance_2backAcc_distribution, paste0(resultFolder, "/performance_2backAcc_distribution.csv"), row.names = F)
saveRDS(modelsum_2backAcc_distribution, paste0(resultFolder, "/modelsum_2backAcc_distribution.rds"))
