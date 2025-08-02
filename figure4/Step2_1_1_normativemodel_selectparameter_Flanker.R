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
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psy/interfileFolder/ABCD"
  functionFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psy/Rcode/functions"
  resultFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psy/results/ABCD"
}else{
  datapath <- 'V:/yunfu_data/EF_results'
  FigureFolder <- 'Z:/DMRI_network_development/EF_development/data0715/Figures_check'
  interfileFolder <- "Z:/DMRI_network_development/EF_development/interfileFolder_check"
  functionFolder <- "Z:/DMRI_network_development/EF_development/data0715/Rcode_EFnorms/functions"
  resultFolder <- "Z:/DMRI_network_development/EF_development/data0715/results_check"
  
}

# source functions
source(paste0(functionFolder, "/Compare_distributions_gamlss_new.R"))

Flanker_data <- read_xlsx(paste0(datapath, "/ABCD_Flanker.xlsx"))
# WM_baseline <- read_xlsx(paste0(datapath, "/WM_baseline.xlsx"))
# WM_development <- read_xlsx(paste0(datapath, "/WM_development.xlsx"))

##### Step1. select the best distribution
# 1. Flanker
dataname <- "Flanker_data"
Flanker_data$Sex <- as.factor(Flanker_data$Sex)
dependentvar <- "nihtbx_flanker_uncorrected"
smoothvar <- "Age_year"
IDvar <- "ID"
bs.df <- 3
covariates <- "Sex"
randomvar=NA

out_Flanker_distribution <- gamlss_comparedistribution(dataname, dependentvar,
                                                     smoothvar, IDvar, bs.df,
                                                     covariates, randomvar=NA)
names(out_Flanker_distribution) <- c("modelsum", "performance")
performance_Flanker_distribution <- out_Flanker_distribution$performance
modelsum_Flanker_distribution <- out_Flanker_distribution$modelsum
write.csv(performance_Flanker_distribution, paste0(resultFolder, "/performance_Flanker_distribution.csv"), row.names = F)
saveRDS(modelsum_Flanker_distribution, paste0(resultFolder, "/modelsum_Flanker_distribution.rds"))

