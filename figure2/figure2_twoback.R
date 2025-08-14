# setting environment
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
rm(list = ls())

# 2. load packages
library(tidyverse)
library(mgcv)
library(psych)
library(gamlss)
library(scales)
library(openxlsx)
library(ggplot2)
library(parallel)

# 3. receive parameters from command
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("error: psychological variable is lost!", call. = FALSE)
}
psyvar_arg <- args[1]

# 4. setting path
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- '/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psy/interfileFolder'
  functionFolder <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/EF_Normative_Model/functions"
  # create subfolder
  resultFolder <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/datasets/yunfu/results/EF_psy"
  FigureFolder <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/datasets/yunfu/figures/fig2"
} else {
  datapath <- '/Users/tanlirou/Documents/YF_EF_psy/EF_psy/correlation/data'
  functionFolder <- "/Users/tanlirou/Documents/YF_EF_psy/EF_psy/code_pure_202507/functions"
  resultFolder <- "/Users/tanlirou/Documents/YF_EF_psy/EF_psy/code_pure_202507/results/Figure2_corr_delete_Z/individual_results"
  FigureFolder <- '/Users/tanlirou/Documents/YF_EF_psy/EF_psy/code_pure_202507/results/Figure2_corr_delete_Z/individual_figures'
}

dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureFolder, showWarnings = FALSE, recursive = TRUE)

source(file.path(functionFolder, "gamcog_withsmoothvar_deviation.R"))

# 5. load data
# GNGd_data <- read_rds(paste0(datapath, '/Gonogo/GNGd_prime.deviations.rds'))
# back1_data <- read_rds(paste0(datapath, '/1-back/back1Acc.deviations.rds'))
back2_data <- read_rds(paste0(datapath, '/2-back/back2Acc.deviations.rds'))


# 6. setting variables
original_vars <- c("SDQ_PB_sum", "SDQ_H_sum", "SDQ_CP_sum", "SDQ_PP_sum", "SDQ_ES_sum")
psyc_variables_continous_all <- paste0(original_vars, "_z")

if (!psyvar_arg %in% psyc_variables_continous_all) {
  stop(paste("The variable '", psyvar_arg, "' is not a valid variable"), call. = FALSE)
}

EFvar <- "Twoback_acc_deviationZ"
dataname <- "twoback"

# 7. data preprocessing
back2_data[, original_vars] <- lapply(back2_data[, original_vars], as.numeric)
schooltab <- unique(back2_data$School)
back2_data$school_fac <- factor(back2_data$School, levels = schooltab, labels = paste0("school", 1:length(schooltab)))

standardize_clean <- function(df, vars) {
  for (var in vars) {
    x <- df[[var]]
    x <- x[!is.na(x)]
    mu <- mean(x, na.rm = TRUE)
    sd_val <- sd(x, na.rm = TRUE)
    valid_index <- which(abs(df[[var]] - mu) <= 3 * sd_val)
    z_varname <- paste0(var, "_z")
    df[[z_varname]] <- NA
    df[[z_varname]][valid_index] <- scale(df[[var]][valid_index])
  }
  return(df)
}

back2_data <- standardize_clean(back2_data, original_vars)

# 8. start analysis
cat(paste("Start process:", psyvar_arg, "\n"))
period <- "EF"
dataname0 <- "back2_data"

dependentvar <- psyvar_arg
interest.indep.var <- EFvar
smoothvar <- "Age_year"
covariates <- "Sex"

result.full <- gam.fit.Independent.var(
  dependentvar = dependentvar,
  dataname = dataname0,
  smoothvar = smoothvar,
  interest.indep.var = interest.indep.var,
  covariates = covariates,
  stats_only = TRUE
)

# 9. save results
result.tmp <- result.full$stats
result.tmp$dataname <- dataname
result.tmp$period <- period
result.tmp$psyvar <- psyvar_arg
output_csv_path <- paste0(resultFolder, "/corr_twoback_", psyvar_arg, ".csv")
write.csv(result.tmp, output_csv_path, row.names = FALSE)
cat(paste("Statistics have been saved to:", output_csv_path, "\n"))


# 10. save simulated results
if (!is.null(result.full$simulation)) {
  sim_result_single <- list(
    EFvar = EFvar,
    dataname = dataname,
    simulation = list(
      psyvar = psyvar_arg,
      simulated_stats = result.full$simulation$simulated_stats,
      observed_stat = result.full$simulation$observed_stat,
      n_sim = result.full$simulation$n_sim
    )
  )
  output_rds_path <- paste0(resultFolder, "/anova_simulation_", dataname0, psyvar_arg, ".rds")
  saveRDS(sim_result_single, output_rds_path)
  cat(paste("simulated results have been saved to:", output_rds_path, "\n"))
  
  # 11. generate distribution figure
  bootstrap_folder <- paste0(FigureFolder, "/bootstrap_distributions")
  dir.create(bootstrap_folder, showWarnings = FALSE, recursive = TRUE)
  
  clean_psyvar <- gsub("[/:*?\"<>|]", "_", psyvar_arg)
  file_path <- paste0(bootstrap_folder, "/sim_dist_twoback_", clean_psyvar, ".png")
  
  png(filename = file_path, width = 800, height = 600, res = 150)
  hist_data <- sim_result_single$simulation$simulated_stats
  observed <- sim_result_single$simulation$observed_stat
  p_value <- mean(hist_data >= observed, na.rm = TRUE)
  
  hist(hist_data, main = paste("Bootstrap Distribution\n", "twoback vs", psyvar_arg),
       xlab = "Deviance Difference", col = "#56B1F7", border = "white", breaks = 50)
  abline(v = observed, col = "#D55E00", lwd = 2.5)
  legend("topright", legend = paste("Observed =", round(observed, 3), "\nP-value =", format.pval(p_value, digits = 3)), bty = "n")
  dev.off()
  
  cat(paste("bootstrap figure have been saved to:", file_path, "\n"))
}

cat(paste("Analysis of variavle ", psyvar_arg, " has completed!\n"))