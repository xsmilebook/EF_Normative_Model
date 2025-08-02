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
  datapath <- 'V:/yunfu_data/EF_results'
  FigureFolder <- 'Z:/DMRI_network_development/EF_development/data0715/Figures_check'
  interfileFolder <- "Z:/DMRI_network_development/EF_development/interfileFolder_check"
  functionFolder <- "Z:/DMRI_network_development/EF_development/data0715/Rcode_EFnorms/functions"
  resultFolder <- "Z:/DMRI_network_development/EF_development/data0715/results_check"
  
}

# source functions
source(paste0(functionFolder, "/Compare_distributions_gamlss_new.R"))

GNG_data <- read_xlsx(paste0(datapath, "/Q_GNG.xlsx"))
back1_data <- read_xlsx(paste0(datapath, "/Q_1back.xlsx"))
back2_data <- read_xlsx(paste0(datapath, "/Q_2back.xlsx"))

# ##### 1.Go No go
# ##### Step2. select proper model details (degree of freedom -- df)
# con<-gamlss.control(n.cyc=200)
# df.set <- matrix(c(3,3,3,
#                    3,4,3,
#                    3,5,3,
#                    3,6,3,
#                    4,3,3,
#                    4,4,3,
#                    4,5,3,
#                    4,6,3,
#                    5,3,3,
#                    5,4,3,
#                    5,5,3,
#                    5,6,3,
#                    6,3,3,
#                    6,4,3,
#                    6,5,3,
#                    6,6,3,
#                    2,2,2,
#                    2,3,2,
#                    2,4,2,
#                    2,5,2,
#                    2,6,2,
#                    3,2,2,
#                    3,3,2,
#                    3,4,2,
#                    3,5,2,
#                    3,6,2,
#                    4,2,2,
#                    4,3,2,
#                    4,4,2,
#                    4,5,2,
#                    4,6,2,
#                    5,3,2,
#                    5,4,2,
#                    5,5,2,
#                    5,6,2,
#                    6,3,2,
#                    6,4,2,
#                    6,5,2,
#                    6,6,2),
#                  byrow=TRUE,ncol=3,dimnames=list(NULL,c("mu","sigma","degree")))
# # 1. Go/no-go
# dataname <- "GNG_data"
# GNG_data$Sex <- as.factor(GNG_data$Sex)
# GNG_data$School <- as.factor(GNG_data$School)
# dependentvar <- "d_prime"
# smoothvar <- "Age_year"
# IDvar <- "ID"
# bs.df <- 3
# saveout_dir <- interfileFolder
# bs.df.set <- df.set
# covariates <- "Sex"
# distribution.fam <- "SEP3"
# out_GNGd_prime_bsdf <- gamlss_compare.bs.df(dataname, dependentvar, smoothvar, IDvar,
#                                             bs.df.set, covariates,
#                                             distribution.fam, randomvar=NA)
# names(out_GNGd_prime_bsdf) <- c("modelsum", "performance")
# performance_GNGd_prime_bsdf <- out_GNGd_prime_bsdf$performance
# modelsum_GNGd_prime_bsdf <- out_GNGd_prime_bsdf$modelsum
# write.csv(performance_GNGd_prime_bsdf, paste0(resultFolder, "/performance_GNGd_prime_bsdf.csv"), row.names = F)
# saveRDS(modelsum_GNGd_prime_bsdf, paste0(resultFolder, "/modelsum_GNGd_prime_bsdf.rds"))
# 
# 
# ##### 2.Oneback
# ##### Step2. select proper model details (degree of freedom -- df)
# con<-gamlss.control(n.cyc=200)
# df.set <- matrix(c(3,3,3,
#                    3,4,3,
#                    3,5,3,
#                    3,6,3,
#                    4,3,3,
#                    4,4,3,
#                    4,5,3,
#                    4,6,3,
#                    5,3,3,
#                    5,4,3,
#                    5,5,3,
#                    5,6,3,
#                    6,3,3,
#                    6,4,3,
#                    6,5,3,
#                    6,6,3,
#                    2,2,2,
#                    2,3,2,
#                    2,4,2,
#                    2,5,2,
#                    2,6,2,
#                    3,2,2,
#                    3,3,2,
#                    3,4,2,
#                    3,5,2,
#                    3,6,2,
#                    4,2,2,
#                    4,3,2,
#                    4,4,2,
#                    4,5,2,
#                    4,6,2,
#                    5,3,2,
#                    5,4,2,
#                    5,5,2,
#                    5,6,2,
#                    6,3,2,
#                    6,4,2,
#                    6,5,2,
#                    6,6,2),
#                  byrow=TRUE,ncol=3,dimnames=list(NULL,c("mu","sigma","degree")))
# # 3. 1back
# dataname <- "back1_data"
# back1_data$Sex <- as.factor(back1_data$Sex)
# back1_data$School <- as.factor(back1_data$School)
# dependentvar <- "Oneback_acc"
# smoothvar <- "Age_year"
# IDvar <- "ID"
# bs.df <- 3
# covariates <- "Sex"
# saveout_dir <- interfileFolder
# distribution.fam <- "SEP2"
# bs.df.set <- df.set
# out_1backAcc_bsdf <- gamlss_compare.bs.df(dataname, dependentvar, smoothvar, IDvar,
#                                           bs.df.set, covariates,
#                                           distribution.fam, randomvar=NA)
# names(out_1backAcc_bsdf) <- c("modelsum", "performance")
# performance_1backAcc_bsdf <- out_1backAcc_bsdf$performance
# modelsum_1backAcc_bsdf <- out_1backAcc_bsdf$modelsum
# write.csv(performance_1backAcc_bsdf, paste0(resultFolder, "/performance_1backAcc_bsdf.csv"), row.names = F)
# saveRDS(modelsum_1backAcc_bsdf, paste0(resultFolder, "/modelsum_1backAcc_bsdf.rds"))



##### 3.Twoback
##### Step2. select proper model details (degree of freedom -- df)
con<-gamlss.control(n.cyc=200)
df.set <- matrix(c(3,3,3,
                   3,4,3,
                   3,5,3,
                   3,6,3,
                   4,3,3,
                   4,4,3,
                   4,5,3,
                   4,6,3,
                   5,3,3,
                   5,4,3,
                   5,5,3,
                   5,6,3,
                   6,3,3,
                   6,4,3,
                   6,5,3,
                   6,6,3,
                   2,2,2,
                   2,3,2,
                   2,4,2,
                   2,5,2,
                   2,6,2,
                   3,2,2,
                   3,3,2,
                   3,4,2,
                   3,5,2,
                   3,6,2,
                   4,2,2,
                   4,3,2,
                   4,4,2,
                   4,5,2,
                   4,6,2,
                   5,3,2,
                   5,4,2,
                   5,5,2,
                   5,6,2,
                   6,3,2,
                   6,4,2,
                   6,5,2,
                   6,6,2),
                 byrow=TRUE,ncol=3,dimnames=list(NULL,c("mu","sigma","degree")))
# 4. 2back
dataname <- "back2_data"
back2_data$Sex <- as.factor(back2_data$Sex)
back2_data$School <- as.factor(back2_data$School)
dependentvar <- "Twoback_acc"
smoothvar <- "Age_year"
IDvar <- "ID"
bs.df <- 3
covariates <- "Sex"
saveout_dir <- interfileFolder
distribution.fam <- "SEP3"
bs.df.set <- df.set
out_2backAcc_bsdf <- gamlss_compare.bs.df(dataname, dependentvar, smoothvar, IDvar,
                                          bs.df.set, covariates,
                                          distribution.fam, randomvar=NA)
names(out_2backAcc_bsdf) <- c("modelsum", "performance")
performance_2backAcc_bsdf <- out_2backAcc_bsdf$performance
modelsum_2backAcc_bsdf <- out_2backAcc_bsdf$modelsum
write.csv(performance_2backAcc_bsdf, paste0(resultFolder, "/performance_2backAcc_bsdf.csv"), row.names = F)
saveRDS(modelsum_2backAcc_bsdf, paste0(resultFolder, "/modelsum_2backAcc_bsdf.rds"))