rm(list = ls())
library(ggplot2)
library(tidyverse)
library(mgcv)
library(readxl)
library(parallel)
library(gamlss)
library(scales)
library(tableone)
library(openxlsx)
# Set file paths
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/EFdevelopment/Crsen_survey_GD_YunFu/data/EF_results'
  demopath <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/EFdevelopment/Crsen_survey_GD_YunFu/data/rawdata_results0616'
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/EFdevelopment/Crsen_survey_GD_YunFu/interfileFolder"
  functionFolder <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/EFdevelopment/Crsen_survey_GD_YunFu/Rcode_EFnorms/functions"
} else {
  datapath <- '/Users/tanlirou/Documents/YF_EF_psy/EF_psy/ABCD/correlation/data'
  FigureFolder <- '/Users/tanlirou/Documents/YF_EF_psy/EF_psy/ABCD/correlation/results_test'
  interfileFolder <- "/Users/tanlirou/Documents/EF_yunfu_check/ABCD/interfileFolder"
  functionFolder <- "/Users/tanlirou/Documents/YF_EF_psy/EF_psy/ABCD/correlation/funtions"
  resultFolder <- "/Users/tanlirou/Documents/YF_EF_psy/EF_psy/ABCD/correlation/results_test"
}

# load data
# Flanker_data <- read_csv(paste0(datapath, '/Flanker.deviations.csv'))
# df1 <- read_xlsx("/Users/tanlirou/Documents/EF_yunfu_check/ABCD/data_combine/ABCD_combine_cbcl.xlsx")
# df1_selected <- df1 %>% select(ID,cbcl_scr_syn_social_r,cbcl_scr_syn_attention_r, cbcl_scr_syn_internal_r,cbcl_scr_syn_external_r,cbcl_scr_syn_totprob_r )
# Flanker_data_new <- Flanker_data %>% left_join(df1_selected, by = "ID")
# write_csv(Flanker_data_new,paste0(datapath, '/New_Flanker.deviations.csv'))
Flanker_data <- read_csv(paste0(datapath, '/Flanker.deviations_addr.csv'))
Flanker_data$Sex <- as.factor(Flanker_data$Sex)
Flanker_data$Sex <- factor(Flanker_data$Sex, levels = c(1, 2), labels = c("M", "F"))
corr_results_file <- paste0(resultFolder, "/corr_EF_psych_continuous.result_ABCD.csv")
corr_results <- read_csv(corr_results_file)
source(paste0(functionFolder,"/gamm_varyingcoefficients_new.R"))
# Define variables
psyc_variables_continous <- c("cbcl_scr_syn_social_r","cbcl_scr_syn_attention_r","cbcl_scr_syn_internal_r","cbcl_scr_syn_external_r")
EFvars.set <- matrix(c("nihtbx_flanker_uncorrected_deviationZ", "Flanker"), byrow=TRUE, ncol=2, dimnames=list(NULL, c("varname", "dataname")))
EFvars.set <- as.data.frame(EFvars.set)

# Set parameters
knots <- 3
set_fx <- TRUE
increments <- 1000
draws <- 1000
return_posterior_coefficients <- T

# Initialize lists to store interaction effects and p-values
int.results.df <- data.frame(
  dependentvar = character(),
  Age = numeric(),
  parcel = character(),
  Interaction = numeric(),
  anova.cov.pvalue = numeric(),
  boots.pvalues = numeric(),
  anovap.fdr = numeric(),
  sig = logical(),
  stringsAsFactors = FALSE
)

# Initialize lists for significant interaction effects
significant_interactions <- data.frame(dependentvar = character(),
                                       Age = numeric(),
                                       parcel = character(),
                                       Interaction = numeric(),
                                       slope_data = I(list()),
                                       stringsAsFactors = FALSE)
nosignificant_interactions <- data.frame(dependentvar = character(),
                                         Age = numeric(),
                                         parcel = character(),
                                         Interaction = numeric(),
                                         slope_data = I(list()),
                                         stringsAsFactors = FALSE)

# Main loop
for (i in 1:nrow(EFvars.set)) {
  int_var <- EFvars.set$varname[i]
  dataname <- paste0(EFvars.set$dataname[i], "_data")
  p_values_task <- c()
  Interactions_task <- c()
  boots.pvalues_task <- c()
  int_pvals_task <- c()
  slope_data_list <- list()
  
  for (j in 1:length(psyc_variables_continous)) {
    dependentvar <- psyc_variables_continous[j]
    smooth_var <- "Age_year"
    covariates <- "Sex"
    bestmodel_row <- corr_results %>% filter(parcel == dependentvar & interest.indep.var == int_var)
    bestmodel <- as.character(bestmodel_row$bettermodel[1])
    Flanker_data <- Flanker_data[!is.na(Flanker_data[[dependentvar]]), ]
    
    result <- gamm.varyingcoefficients(dependentvar, dataname, smooth_var, int_var, covariates, bestmodel, knots, set_fx, increments, draws, return_posterior_coefficients)
    # Extract interaction effect p-values and estimates
    p_value <- as.numeric(result[[1]][1, "anova.int.pvalue"])
    Interaction <- as.numeric(result[[1]][1, "IntpartialRsq"])
    boots.pvalues <- as.numeric(result[[1]][1, "boots.pvalues"])
    slope_data <- result[[2]] 
    slope_data_list[[j]] <- data.frame(Age_year = slope_data[[smooth_var]], slope = slope_data$slope)
    
    # Store values in the lists
    p_values_task <- c(p_values_task, p_value)
    Interactions_task <- c(Interactions_task, Interaction)
    boots.pvalues_task <- c(boots.pvalues_task, boots.pvalues)
    
  
    if ("model" %in% names(result)) {
      modelobj <- result$model
      if (bestmodel == "GAM") {
        int_pval <- pbootint(modelobj, int_var)  #For GAMM, use pbootint
      } else if (bestmodel == "LM") {
        int_pval <- NA  #  Or use a specific test for LM
      }
    }
    
    temp_df <- data.frame(
      dependentvar = dependentvar,
      Age = smooth_var,
      parcel = int_var,
      Interaction = Interaction,
      anova.cov.pvalue = p_value,
      boots.pvalues = boots.pvalues,
      anovap.fdr = NA,
      sig = FALSE,
      stringsAsFactors = FALSE
    )
    
    int.results.df <- rbind(int.results.df, temp_df)
  }
  
  #FDR correction and update significance flags
  fdr_corrected_p <- p.adjust(boots.pvalues_task, method = "fdr")
  int.results.df[int.results.df$parcel == int_var, "anovap.fdr"] <- fdr_corrected_p
  int.results.df[int.results.df$parcel == int_var, "sig"] <- (fdr_corrected_p < 0.05)
  
  sig_indices <- which(fdr_corrected_p < 0.05)
  if (length(sig_indices) > 0) {
    for (idx in sig_indices) {
      significant_interactions <- rbind(significant_interactions, 
                                        data.frame(dependentvar = psyc_variables_continous[idx], 
                                                   Age = smooth_var,
                                                   parcel = int_var, 
                                                   Interaction = Interactions_task[idx],
                                                   slope_data = I(list(slope_data_list[[idx]]))))
    }
  }
  nosig_indices <- which(fdr_corrected_p > 0.05)
      if (length(nosig_indices) > 0) {
        for (idx in nosig_indices) {
          nosignificant_interactions <- rbind(nosignificant_interactions, 
                                              data.frame(dependentvar = psyc_variables_continous[idx], 
                                                         Age = smooth_var,
                                                         parcel = int_var, 
                                                         Interaction = Interactions_task[idx],
                                                         slope_data = I(list(slope_data_list[[idx]]))))
    }
  }
  write.xlsx(int.results.df, file = paste0(resultFolder, "/interaction_results_full.xlsx"), row.names = FALSE)
}
write_rds(significant_interactions, file = paste0(resultFolder, "/int_significant_ABCD.rds"))
write_rds(nosignificant_interactions, file = paste0(resultFolder, "/int_Nosignificant_ABCD.rds"))

int.results.df <- readRDS(paste0(resultFolder, "/int_significant_ABCD.rds"))
significant_interactions <- readRDS(paste0(resultFolder, "/int_significant_ABCD.rds"))

# 绘制热图
lwth <- min(int.results.df$Interaction, na.rm = TRUE)
upth <- max(int.results.df$Interaction, na.rm = TRUE)
y_levels <- c("cbcl_scr_syn_social_r","cbcl_scr_syn_attention_r","cbcl_scr_syn_internal_r","cbcl_scr_syn_external_r")

for (i in 1:nrow(EFvars.set)) {
  EFvar.tmp <- EFvars.set$varname[i]
  dataname <- EFvars.set$dataname[i]
  corr.result.tmp <- int.results.df[which(int.results.df$parcel == EFvar.tmp & int.results.df$dataname == dataname), ]
  corr.result.tmp$sig <- as.logical(corr.result.tmp$sig)  # 转换 sig 列为逻辑值
  corr.result.tmp.sig <- corr.result.tmp[corr.result.tmp$sig == TRUE, ]

  # Plot heatmap
  Fig <- ggplot() +
    geom_tile(data = corr.result.tmp, aes(x = parcel, y = dependentvar, fill = Interaction), color = "white") +
    geom_text(data = corr.result.tmp.sig, aes(x = parcel, y = dependentvar, label = "*"), vjust = 0.75, hjust = 0.5, size = 9) +
    scale_fill_distiller(type = "seq", palette = "Reds", limits = c(lwth, upth), direction = 1) +
    scale_y_discrete(limits = y_levels,
                     labels = c("cbcl_scr_syn_social_r"= "Social","cbcl_scr_syn_attention_r" = "Attention", 
                                "cbcl_scr_syn_internal_r"= "Internal","cbcl_scr_syn_external_r" = "External")) +
    scale_x_discrete(labels = c("nihtbx_flanker_uncorrected_deviationZ" = "Flanker")) +
    labs(title = paste0("Interaction Effects in ", dataname),
         x = "Tasks", y = "Mental Health") +
    theme(axis.line = element_blank(),
          aspect.ratio = 1.2,
          axis.text.x = element_text(size = 12, hjust = 0.5),
          axis.text.y = element_text(size = 12, hjust = 0.5, vjust = 0.5),
          axis.title = element_text(size = 18),
          plot.title = element_text(size = 18, hjust = 0.5),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          panel.background = element_rect(fill = NA),
          panel.grid.major = element_line(linewidth = 0),
          panel.grid.minor = element_line(linewidth = 1))

  print(Fig)
  ggsave(paste0(FigureFolder, "/corr_int_", EFvar.tmp, "_tasks.pdf"), plot = Fig, width = 16, height = 12, units = "cm")
}


# Customize the order of scales
y_levels <- c("cbcl_scr_syn_social_r","cbcl_scr_syn_attention_r","cbcl_scr_syn_internal_r","cbcl_scr_syn_external_r")  #Scale order

# Create a combined task data frame
combined_data <- int.results.df
combined_data$sig <- as.logical(combined_data$sig)  #Convert sig column to logical values
combined_data$parcel <- factor(combined_data$parcel, levels = c("nihtbx_flanker_uncorrected_deviationZ"), 
                               labels = c( "Flanker"))  # Set task order and labels

####post plot
# Iterate through all significant interaction effects
variable_mapping <- c(
  "cbcl_scr_syn_social_r"= "Social","cbcl_scr_syn_attention_r" = "Attention", 
  "cbcl_scr_syn_internal_r"= "Internal","cbcl_scr_syn_external_r" = "External"
)
for (k in 1:nrow(significant_interactions)) {
  # Get the data for the current significant result
  current_slope_data <- significant_interactions$slope_data[[k]]
  current_task <- significant_interactions$dataname[k]
  current_variable <- significant_interactions$dependentvar[k]
  
  modified_variable_name <- ifelse(current_variable %in% names(variable_mapping), variable_mapping[current_variable], current_variable)
  # Calculate the median and confidence intervals
  slope_summary <- current_slope_data %>%
    group_by(Age_year) %>%
    summarise(
      median_slope = median(slope, na.rm = TRUE),
      lower_95CI = quantile(slope, probs = 0.025, na.rm = TRUE),
      upper_95CI = quantile(slope, probs = 0.975, na.rm = TRUE)
    )
  
  # Plot the slope changes over age
  slope_plot <- ggplot(slope_summary, aes(x = Age_year)) +
    geom_line(aes(y = median_slope), color = "black", size = 1) + 
    geom_ribbon(aes(ymin = lower_95CI, ymax = upper_95CI), alpha = 0.3, fill = "grey80") + 
    labs(title = paste0(modified_variable_name, " ~ Flanker"), 
         x = "", 
         y = "Slope") +
    scale_x_continuous(name="", limits = c(8.8, 16),breaks = seq(9, 16, by = 1))+
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.5, color = "black"),  
      axis.text = element_text(size = 8.5,color = "black"), 
      axis.title = element_text(size = 8.5), 
      plot.title = element_text(size = 8.5, hjust = 0.5),
      panel.grid.major = element_blank(),    
      panel.grid.minor = element_blank(),    
      panel.border = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.5),  
      axis.ticks.length = unit(0.2, "cm") 
    )
  
  print(slope_plot)
  ggsave(paste0(FigureFolder, "/significant_slope_change_", current_task, "_", current_variable, ".pdf"), plot = slope_plot, width = 9, height = 7, units = "cm")
}


