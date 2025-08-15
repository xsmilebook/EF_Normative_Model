# 1. setting environment
# Sys.setenv(OPENBLAS_NUM_THREADS = 1) 
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

# 3. setting path
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- '/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psy/interfileFolder'
  functionFolder <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/EF_Normative_Model/functions"
  resultFolder <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/datasets/yunfu/results/EF_psy"
  FigureFolder <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/datasets/yunfu/figures/fig2"
} else {
  datapath <- 'D:/datasets/yunfu/raw_data'
  functionFolder <- "D:/code/EF_Normative_Model/functions"
  interfileFolder <- "D:/datasets/yunfu/interfile_folder/temp"
  resultFolder <- "D:/datasets/yunfu/results/EF_psy"
  FigureFolder <- 'D:/datasets/yunfu/results/figures/fig2'
}

dir.create(resultFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(FigureFolder, showWarnings = FALSE, recursive = TRUE)

# 4. source the custom function
source(file.path(functionFolder, "gamcog_withsmoothvar_deviation.R"))

# 5. load data 
# GNGd_data <- read_rds(paste0(datapath, '/Gonogo/GNGd_prime.deviations.rds'))
back2_data <- read_rds(paste0(interfileFolder, '/back2Acc.deviations.rds'))

# 6. setting variables for loops
Time_ids <- paste0("Time_", 0:9)
# original_vars <- c("SDQ_PB_sum", "SDQ_H_sum", "SDQ_CP_sum", "SDQ_PP_sum", "SDQ_ES_sum")
original_vars <- c("SDQ_PP_sum")
psyc_variables_continous_all <- paste0(original_vars, "_z")

EFvar <- "Twoback_acc_deviationZ"
dataname <- "twoback"
dataname0 <- "back2_data" 

# 7. data preprocessing 
for (var in original_vars) {
  back2_data[[var]] <- as.numeric(back2_data[[var]])
}
schooltab <- unique(back2_data$School)
back2_data$school_fac <- factor(back2_data$School, levels = schooltab, labels = paste0("school", 1:length(schooltab)))

standardize_clean <- function(df, vars) {
  for (var in vars) {
    x <- df[[var]]
    if (sum(!is.na(x)) > 1) {
      x_clean <- x[!is.na(x)]
      mu <- mean(x_clean, na.rm = TRUE)
      sd_val <- sd(x_clean, na.rm = TRUE)
      if (!is.na(sd_val) && sd_val > 0) {
        valid_index <- which(!is.na(df[[var]]) & abs(df[[var]] - mu) <= 3 * sd_val)
        z_varname <- paste0(var, "_z")
        df[[z_varname]] <- NA
        df[[z_varname]][valid_index] <- scale(df[[var]][valid_index])
      } else {
        df[[paste0(var, "_z")]] <- NA
      }
    } else {
      df[[paste0(var, "_z")]] <- NA
    }
  }
  return(df)
}

back2_data <- standardize_clean(back2_data, original_vars)

# 8. start analysis loops
# ==============================================================================
# outer loop
for (Time_id in Time_ids) {
  
  # ============================================================================
  # inner loop
  for (psyvar_arg in psyc_variables_continous_all) {
    
    cat(paste("\n\n===== Starting Analysis for Time_id:", Time_id, "| PsyVar:", psyvar_arg, "=====\n"))
    
    if (!psyvar_arg %in% names(back2_data) || all(is.na(back2_data[[psyvar_arg]]))) {
      cat(paste("Skipping:", psyvar_arg, "- Variable not found or all values are NA.\n"))
      next
    }
    
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
    
    # 9. save statistics results
    result.tmp <- result.full$stats
    result.tmp$dataname <- dataname
    result.tmp$period <- "EF" 
    result.tmp$psyvar <- psyvar_arg
    result.tmp$Time_id <- Time_id 
    
    output_csv_path <- file.path(resultFolder, paste0("corr_", dataname, "_", psyvar_arg, "_", Time_id, ".csv"))
    write.csv(result.tmp, output_csv_path, row.names = FALSE)
    cat(paste("Statistics have been saved to:", output_csv_path, "\n"))
    
    # 10. save simulated results
    if (!is.null(result.full$simulation) && length(result.full$simulation$simulated_stats) > 0) {
      sim_result_single <- list(
        EFvar = EFvar,
        dataname = dataname,
        simulation = list(
          psyvar = psyvar_arg,
          Time_id = Time_id,
          simulated_stats = result.full$simulation$simulated_stats,
          observed_stat = result.full$simulation$observed_stat,
          n_sim = result.full$simulation$n_sim
        )
      )
      
      output_rds_path <- file.path(resultFolder, paste0("anova_simulation_", dataname, "_", psyvar_arg, "_", Time_id, ".rds"))
      saveRDS(sim_result_single, output_rds_path)
      cat(paste("Simulated results have been saved to:", output_rds_path, "\n"))
      
      # 11. generate distribution figure
      bootstrap_folder <- file.path(FigureFolder, "bootstrap_distributions")
      dir.create(bootstrap_folder, showWarnings = FALSE, recursive = TRUE)
      
      clean_psyvar <- gsub("[/:*?\"<>|]", "_", psyvar_arg)
      file_path <- file.path(bootstrap_folder, paste0("sim_dist_", dataname, "_", clean_psyvar, "_", Time_id, ".png"))
      
      png(filename = file_path, width = 800, height = 600, res = 150)
      hist_data <- sim_result_single$simulation$simulated_stats
      observed <- sim_result_single$simulation$observed_stat
      p_value <- mean(hist_data >= observed, na.rm = TRUE)
      
      hist(hist_data, main = paste("Bootstrap Distribution\n", dataname, "vs", psyvar_arg, "(", Time_id, ")"),
           xlab = "Deviance Difference", col = "#56B1F7", border = "white", breaks = 50)
      abline(v = observed, col = "#D55E00", lwd = 2.5)
      legend("topright", legend = paste("Observed =", round(observed, 3), "\nP-value =", format.pval(p_value, digits = 3, eps = 0.001)), bty = "n")
      dev.off()
      
      cat(paste("Bootstrap figure has been saved to:", file_path, "\n"))
      
    } else {
      cat("No simulation results to save or plot.\n")
    }
    
    cat(paste("Analysis for", psyvar_arg, "at", Time_id, "has completed!\n"))

    if(exists("result.full")) rm(result.full)
    if(exists("result.tmp")) rm(result.tmp)
    if(exists("sim_result_single")) rm(sim_result_single)
    if(exists("hist_data")) rm(hist_data)

    gc() 
    
  } 
} 

cat("\n\nAll analyses have been completed!\n")