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













