rm(list=ls())
library(ggplot2)
library(tidyverse)
library(mgcv)
library(readxl)
library(parallel)
library(gamlss)
library(scales)
library(tableone)
library(openxlsx)
library(ggradar)
# input directory
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/EFdevelopment/Crsen_survey_GD_YunFu/data/EF_results'
  demopath <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/EFdevelopment/Crsen_survey_GD_YunFu/data/rawdata_results0616'
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/EFdevelopment/Crsen_survey_GD_YunFu/interfileFolder"
  functionFolder <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/EFdevelopment/Crsen_survey_GD_YunFu/Rcode_EFnorms/functions"
  
}else{
  datapath <- '/Users/tanlirou/Documents/EF_yunfu_check/correlation/data'
  FigureFolder <- '/Users/tanlirou/Documents/EF_yunfu_check/correlation/results'
  interfileFolder <- "/Users/tanlirou/Documents/EF_yunfu_check/correlation/data"
  functionFolder <- "/Users/tanlirou/Documents/EF_yunfu_check/correlation/code/functions"
  resultFolder <- "/Users/tanlirou/Documents/EF_yunfu_check/correlation/results"
}

# source functions
#source(paste0(functionFolder, "/gam_varyingcoefficients.R"))
source(paste0(functionFolder, "/gamcog_new.R"))
source(paste0(functionFolder, "/ordinalcorr_new.R"))
# import dataset
#switch_data <- read_xlsx(paste0(datapath, '/Q_switch.xlsx'))
#head(switch_data)
GNGd_data <- read_csv(paste0(datapath, '/GNGd/GNGd_prime.deviations.csv'))
back1_data <- read_csv(paste0(datapath, '/1backAcc/back1Acc.deviations.csv'))
back2_data <- read_csv(paste0(datapath, '/2backACC/back2Acc.deviations.csv'))


## 1) set up variables
psyc_variables_continous <- c("SDQ_PB_sum", "SDQ_H_sum", "SDQ_CP_sum", "SDQ_PP_sum", "SDQ_ES_sum", "SDQ_sum")
psyc_variables_discrete <- c( "SDQ_cutoff")
# EF vars
EFvars.set <- matrix(c("d_prime_deviationZ", "GNGd",
                       "Oneback_acc_deviationZ", "back1",
                       "Twoback_acc_deviationZ", "back2"), byrow=TRUE,ncol=2,dimnames=list(NULL,c("varname","dataname")))
EFvars.set <- as.data.frame(EFvars.set)
## 2) convert variables class & describe variables

GNGd_data[,psyc_variables_discrete] <- lapply(GNGd_data[,psyc_variables_discrete], as.factor)
GNGd_data[,psyc_variables_continous] <- lapply(GNGd_data[,psyc_variables_continous], as.numeric)
schooltab <- unique(GNGd_data$School)
GNGd_data$school_fac <- factor(GNGd_data$School, levels=schooltab, labels=paste0("school", 1:length(schooltab)))
# 
# GNGCV_data[,psyc_variables_discrete] <- lapply(GNGCV_data[,psyc_variables_discrete], as.factor)
# GNGCV_data[,psyc_variables_continous] <- lapply(GNGCV_data[,psyc_variables_continous], as.numeric)
# schooltab <- unique(GNGCV_data$School)
# GNGCV_data$school_fac <- factor(GNGCV_data$School, levels=schooltab, labels=paste0("school", 1:length(schooltab)))

back1_data[,psyc_variables_discrete] <- lapply(back1_data[,psyc_variables_discrete], as.factor)
back1_data[,psyc_variables_continous] <- lapply(back1_data[,psyc_variables_continous], as.numeric)
schooltab <- unique(back1_data$School)
back1_data$school_fac <- factor(back1_data$School, levels=schooltab, labels=paste0("school", 1:length(schooltab)))

back2_data[,psyc_variables_discrete] <- lapply(back2_data[,psyc_variables_discrete], as.factor)
back2_data[,psyc_variables_continous] <- lapply(back2_data[,psyc_variables_continous], as.numeric)
schooltab <- unique(back2_data$School)
back2_data$school_fac <- factor(back2_data$School, levels=schooltab, labels=paste0("school", 1:length(schooltab)))
# describe
# GNGd
describe_tab_GNGd <- CreateTableOne(c(EFvars.set$varname[EFvars.set$dataname=="GNGd"], psyc_variables_continous, psyc_variables_discrete, "Age_year", "Sex"), data=GNGd_data, factorVars = psyc_variables_discrete, testNonNormal=T)
describe_tab_GNGd.continous <- as.data.frame(describe_tab_GNGd$ContTable[["Overall"]])
describe_tab_GNGd.discrete <- do.call(rbind, lapply(describe_tab_GNGd$CatTable[["Overall"]], as.data.frame))
# # GNGCV
# describe_tab_GNGCV <- CreateTableOne(c(EFvars.set$varname[EFvars.set$dataname=="GNGCV"], psyc_variables_continous, psyc_variables_discrete, "Age_year", "Sex"), data=GNGCV_data, factorVars = psyc_variables_discrete, testNonNormal=T)
# describe_tab_GNGCV.continous <- as.data.frame(describe_tab_GNGCV$ContTable[["Overall"]])
# describe_tab_GNGCV.discrete <- do.call(rbind, lapply(describe_tab_GNGCV$CatTable[["Overall"]], as.data.frame))
# back1
describe_tab_back1 <- CreateTableOne(c(EFvars.set$varname[EFvars.set$dataname=="back1"], psyc_variables_continous, psyc_variables_discrete, "Age_year", "Sex"), data=back1_data, factorVars = psyc_variables_discrete, testNonNormal=T)
describe_tab_back1.continous <- as.data.frame(describe_tab_back1$ContTable[["Overall"]])
describe_tab_back1.discrete <- do.call(rbind, lapply(describe_tab_back1$CatTable[["Overall"]], as.data.frame))
# back2
describe_tab_back2 <- CreateTableOne(c(EFvars.set$varname[EFvars.set$dataname=="back2"], psyc_variables_continous, psyc_variables_discrete, "Age_year", "Sex"), data=back2_data, factorVars = psyc_variables_discrete, testNonNormal=T)
describe_tab_back2.continous <- as.data.frame(describe_tab_back2$ContTable[["Overall"]])
describe_tab_back2.discrete <- do.call(rbind, lapply(describe_tab_back2$CatTable[["Overall"]], as.data.frame))
# save out
write.xlsx(list(GNGd_con=describe_tab_GNGd.continous,GNGd_dis=describe_tab_GNGd.discrete,back1_con=describe_tab_back1.continous,back1_dis=describe_tab_back1.discrete,back2_con=describe_tab_back2.continous,back2_dis=describe_tab_back2.discrete), paste0(resultFolder, "/description_interest_vars.xlsx"), rowNames=T)

## 3) Correlations in separate age periods
# continuos variables
n=0
corr.result.period <- list()
for (i in 1:nrow(EFvars.set)){
  EFvar.tmp <- EFvars.set$varname[i]
  dataname0 <- paste0(EFvars.set$dataname[i], "_data")
  data.tmp <- get(dataname0)
  data.tmp.EF <- data.tmp %>% filter(Age_year > 11 & Age_year <= 18)
  
  for (period in c("EF")){
    dataname <- paste0("data.tmp.", period)
    corr.result <- list()
    
    # 获取数据段并检查非NA值数量
    data_segment <- get(dataname)
    if (sum(!is.na(data_segment)) < 30) { 
      corr.result.df <- data.frame(period = period, correlation = NA)
      corr.result.period[[n]] <- corr.result.df
      next
    }
    
    for (x in 1:length(psyc_variables_continous)){
      psyvar.tmp <- psyc_variables_continous[x]
      dependentvar <- psyvar.tmp
      interest.indep.var <- EFvar.tmp
      covariates <- "Sex"
      knots=3
      
      result.tmp <- gam.fit.Independent.var(dependentvar, dataname,interest.indep.var, covariates, stats_only = T)
      result.tmp <- as.data.frame(result.tmp)
      result.tmp$dataname <- EFvars.set$dataname[i]
      corr.result[[x]] <- result.tmp
    }
    corr.result.df <- do.call(rbind, corr.result)
    corr.result.df$period <- period
    n=n+1
    corr.result.period[[n]] <- corr.result.df
  }
}
corr.result.period.df.con <- do.call(rbind, corr.result.period)
write.csv(corr.result.period.df.con, paste0(resultFolder, "/corr_EF_psych_continuous.result.csv"), row.names = F)
