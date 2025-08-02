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

# set path
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/EFdevelopment/Crsen_survey_GD_YunFu/data/EF_results'
  demopath <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/EFdevelopment/Crsen_survey_GD_YunFu/data/rawdata_results0616'
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/EFdevelopment/Crsen_survey_GD_YunFu/interfileFolder"
  functionFolder <- "/ibmgpfs/cuizaixu_lab/xuxiaoyu/EFdevelopment/Crsen_survey_GD_YunFu/Rcode_EFnorms/functions"
  
} else {
  datapath <- '/Users/tanlirou/Documents/EF_yunfu_check/correlation/data'
  FigureFolder <- '/Users/tanlirou/Documents/EF_yunfu_check/correlation/results_int'
  interfileFolder <- "/Users/tanlirou/Documents/EF_yunfu_check/correlation/data"
  functionFolder <- "/Users/tanlirou/Documents/EF_yunfu_check/correlation/code/functions"
  resultFolder <- "/Users/tanlirou/Documents/EF_yunfu_check/correlation/results_int"
}

# load data
GNGd_data <- read_csv(paste0(datapath, '/GNGd/GNGd_prime.deviations.csv'))
back1_data <- read_csv(paste0(datapath, '/1backAcc/back1Acc.deviations.csv'))
back2_data <- read_csv(paste0(datapath, '/2backACC/back2Acc.deviations.csv'))
GNGd_data$Sex <- as.factor(GNGd_data$Sex)
back1_data$Sex <- as.factor(back1_data$Sex)
back2_data$Sex <- as.factor(back2_data$Sex)
source("~/Documents/EF_yunfu_check/code/functions/gam_varyingcoefficients.R")
# Define variables
psyc_variables_continous <- c("SDQ_PB_sum", "SDQ_H_sum", "SDQ_CP_sum", "SDQ_PP_sum", "SDQ_ES_sum", "SDQ_sum")
EFvars.set <- matrix(c("d_prime_deviationZ", "GNGd",
                       "Oneback_acc_deviationZ", "back1",
                       "Twoback_acc_deviationZ", "back2"), byrow=TRUE, ncol=2, dimnames=list(NULL, c("varname", "dataname")))
EFvars.set <- as.data.frame(EFvars.set)

# Set parameters
knots <- 3
set_fx <- FALSE 
increments <- 1000
draws <- 1000
return_posterior_coefficients <- T

# Initialize list to store interaction effects and p-values
int.results.df <- data.frame(
  dependentvar = character(),
  Age = numeric(),
  parcel = character(),
  dataname = character(),
  Interaction = numeric(),
  anova.cov.pvalue = numeric(),
  anova.pvalue = numeric(),
  anovap.fdr = numeric(),
  sig = logical(),
  stringsAsFactors = FALSE
)

# Initialize list to store significant interactions
significant_interactions <- data.frame(dependentvar = character(),
                                       Age = numeric(),
                                       parcel = character(),
                                       dataname = character(),
                                       Interaction = numeric(),
                                       slope_data = I(list()),
                                       stringsAsFactors = FALSE)
nosignificant_interactions <- data.frame(dependentvar = character(),
                                         Age = numeric(),
                                         parcel = character(),
                                         dataname = character(),
                                         Interaction = numeric(),
                                         slope_data = I(list()),
                                         stringsAsFactors = FALSE)

#Main loop
for (i in 1:nrow(EFvars.set)) {
  int_var <- EFvars.set$varname[i]
  dataname <- paste0(EFvars.set$dataname[i], "_data")
  p_values_task <- c()
  anova.pvalues_task <- c()
  Interactions_task <- c()
  slope_data_list <- list()
  
  for (j in 1:length(psyc_variables_continous)) {
    dependentvar <- psyc_variables_continous[j]
    smooth_var <- "Age_year"
    result <- gam.varyingcoefficients(dependentvar, dataname, smooth_var, int_var, covariates = "Sex", knots, set_fx, increments, draws, return_posterior_coefficients)
    # Extract p-values and estimates for interaction effects
    p_value <- as.numeric(result[[1]][1, "anova.int.pvalue"])
    Interaction <- as.numeric(result[[1]][1, "IntpartialRsq"])
    anova.pvalues <- as.numeric(result[[1]][1, "anova.pvalues"])
    slope_data <- result[[2]] 
    slope_data_list[[j]] <- data.frame(Age_year = slope_data[[smooth_var]], slope = slope_data$slope)
    #  Store values in the list
    p_values_task <- c(p_values_task, p_value)
    anova.pvalues_task <- c(anova.pvalues_task, anova.pvalues)
    Interactions_task <- c(Interactions_task, Interaction)

    temp_df <- data.frame(
      dependentvar = dependentvar,
      Age = smooth_var,
      parcel = int_var,
      dataname = EFvars.set$dataname[i],
      Interaction = Interaction,
      anova.cov.pvalue = p_value,
      anova.pvalues = anova.pvalues,
      anovap.fdr = NA,
      sig = FALSE,
      stringsAsFactors = FALSE
    )
    
    int.results.df <- rbind(int.results.df, temp_df)
  }
  
  # FDR correction and update significance markers
  fdr_corrected_p <- p.adjust(anova.pvalues_task, method = "fdr")
  int.results.df[int.results.df$parcel == int_var, "anovap.fdr"] <- fdr_corrected_p
  int.results.df[int.results.df$parcel == int_var, "sig"] <- (fdr_corrected_p < 0.05)
  
  sig_indices <- which(fdr_corrected_p < 0.05)
  if (length(sig_indices) > 0) {
    for (idx in sig_indices) {
      significant_interactions <- rbind(significant_interactions, 
                                        data.frame(dependentvar = psyc_variables_continous[idx], 
                                                   Age = smooth_var,
                                                   parcel = int_var, 
                                                   dataname = EFvars.set$dataname[i], 
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
                                                   dataname = EFvars.set$dataname[i], 
                                                   Interaction = Interactions_task[idx],
                                                   slope_data = I(list(slope_data_list[[idx]]))))
    }
  }
  write.xlsx(int.results.df, file = paste0(resultFolder, "/interaction_results_full.xlsx"), row.names = FALSE)
}


# Plot heatmap
lwth <- min(int.results.df$Interaction, na.rm = TRUE)
upth <- max(int.results.df$Interaction, na.rm = TRUE)
y_levels <- c("SDQ_PB_sum", "SDQ_H_sum", "SDQ_CP_sum", "SDQ_PP_sum", "SDQ_ES_sum", "SDQ_sum")

for (i in 1:nrow(EFvars.set)) {
  EFvar.tmp <- EFvars.set$varname[i]
  dataname <- EFvars.set$dataname[i]
  corr.result.tmp <- int.results.df[which(int.results.df$parcel == EFvar.tmp & int.results.df$dataname == dataname), ]
  corr.result.tmp$sig <- as.logical(corr.result.tmp$sig)  # Convert sig column to logical
  corr.result.tmp.sig <- corr.result.tmp[corr.result.tmp$sig == TRUE, ]
  # Plot heatmap
  Fig <- ggplot() +
    geom_tile(data = corr.result.tmp, aes(x = parcel, y = dependentvar, fill = Interaction), color = "white") +
    geom_text(data = corr.result.tmp.sig, aes(x = parcel, y = dependentvar, label = "*"), vjust = 0.75, hjust = 0.5, size = 9) +
    scale_fill_distiller(type = "seq", palette = "Reds", limits = c(lwth, upth), direction = 1) +
    scale_y_discrete(limits = y_levels,
                     labels = c("SDQ_PB_sum" = "Prosocial Behavior","SDQ_H_sum" = "Hyperactivity",
                                "SDQ_CP_sum" = "Conduct Problems","SDQ_PP_sum" = "Peer Problems",
                                "SDQ_ES_sum" = "Emotional Symptoms","SDQ_sum" = "Total SDQ Score")) +
    scale_x_discrete(labels = c("d_prime_deviationZ" = "Go/No-go",
                                "Oneback_acc_deviationZ" = "1-back",
                                "Twoback_acc_deviationZ" = "2-back")) +
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
  #ggsave(paste0(FigureFolder, "/corr_int_", EFvar.tmp, "_tasks.pdf"), plot = Fig, width = 16, height = 12, units = "cm")
}


# Determine the upper and lower limits of the color scale
lwth <- min(int.results.df$Interaction, na.rm = TRUE)
upth <- max(int.results.df$Interaction, na.rm = TRUE)

# Custom order for the scales
y_levels <- c("SDQ_sum", "SDQ_ES_sum", "SDQ_PP_sum", "SDQ_CP_sum", "SDQ_H_sum", "SDQ_PB_sum")  # Order of the scales

# Create the combined task data frame
combined_data <- int.results.df
combined_data$sig <- as.logical(combined_data$sig)  # Convert sig column to logical values
combined_data$parcel <- factor(combined_data$parcel, levels = c( "d_prime_deviationZ", "Oneback_acc_deviationZ", "Twoback_acc_deviationZ"), 
                               labels = c("Go/No-go ", "1-back", "2-back"))  # Set task order and labels

# Define custom color gradient
custom_colors <- c( "#f5dfdb", "#edb8b0", "#e69191", "#c25759")  

# Plot the combined heatmap
Fig <- ggplot() +
  geom_tile(data = combined_data, aes(x = parcel, y = dependentvar, fill = Interaction), color = "white") +  
  geom_text(data = combined_data[combined_data$sig == TRUE, ], 
            aes(x = parcel, y = dependentvar, label = "*"), 
            vjust = 0.8, hjust = 0.5, size = 9) +
  scale_fill_gradientn(colors = custom_colors, limits = c(lwth, upth)) +  # Use custom color gradient
  scale_y_discrete(limits = y_levels,  
                   labels = c("SDQ_sum" = "SDQ Total Score","SDQ_ES_sum" = "Emotional Symptoms", 
                              "SDQ_PP_sum" = "Peer Problems","SDQ_CP_sum" = "Conduct Problems",
                              "SDQ_H_sum" = "Hyperactivity","SDQ_PB_sum" = "Prosocial Behavior")) +
  labs(title = "Correlation of Executive Functions and Mental Health across Age", 
       x = "Tasks", y = "Mental Health") +
  theme(axis.line = element_blank(),
        aspect.ratio = 1.6,
        axis.text.x = element_text(size = 14, hjust = 0.5),
        axis.text.y = element_text(size = 14, hjust = 0.5, vjust = 0.5),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 22, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linewidth = 0),
        panel.grid.minor = element_line(linewidth = 1))
# Save plot
print(Fig)
ggsave(paste0(FigureFolder, "/corr_int.pdf"), plot = Fig, width = 18, height =24 , units = "cm")


# Custom order for the scales
y_levels <- c("APSS_sum", "SDQ_PB_sum", "SDQ_H_sum", "SDQ_CP_sum", "SDQ_PP_sum", "SDQ_ES_sum", "SDQ_sum")  # Order of the scales

# Create the combined task data frame
combined_data <- int.results.df
combined_data$sig <- as.logical(combined_data$sig)  # Convert sig column to logical values
combined_data$parcel <- factor(combined_data$parcel, levels = c("Twoback_acc_deviationZ", "Oneback_acc_deviationZ",  "d_prime_deviationZ"), 
                               labels = c("2back", "1back", "Go/No-go"))  # Set task order and labels

#### Post plot
# Loop through all significant interactions
for (k in 1:nrow(significant_interactions)) {
  # Get the data for the current significant result
  current_slope_data <- significant_interactions$slope_data[[k]]
  current_task <- significant_interactions$dataname[k]
  current_variable <- significant_interactions$dependentvar[k]
  
  # Calculate the median and confidence intervals
  slope_summary <- current_slope_data %>%
    group_by(Age_year) %>%
    summarise(
      median_slope = median(slope, na.rm = TRUE),
      lower_95CI = quantile(slope, probs = 0.025, na.rm = TRUE),
      upper_95CI = quantile(slope, probs = 0.975, na.rm = TRUE)
    )
  
  # Plot the trend of slope changes with age
  slope_plot <- ggplot(slope_summary, aes(x = Age_year)) +
    geom_line(aes(y = median_slope), color = "black", size = 1) + 
    geom_ribbon(aes(ymin = lower_95CI, ymax = upper_95CI), alpha = 0.3, fill = "grey80") + 
    labs(title = paste0(modified_variable_name, " ~ ", modified_task_name), 
         x = "", 
         y = "Slope") +
    scale_x_continuous(name="", limits = c(11, 18),breaks = seq(11, 18, by = 1))+
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.5, color = "black"),  
      axis.text = element_text(size = 8.5,color = "black"), 
      axis.title = element_text(size = 8.5), 
      plot.title = element_text(size = 8.5, hjust = 0.5),
      panel.grid.major = element_blank(),    # 删除主网格线
      panel.grid.minor = element_blank(),    # 删除次网格线
      panel.border = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.5),  # 设置小线段颜色和大小
      axis.ticks.length = unit(0.2, "cm") 
    )
  
  print(slope_plot)
  # Save plot
  ggsave(paste0(FigureFolder, "/significant_slope_change_", current_task, "_", current_variable, ".pdf"), plot = slope_plot, width = 16, height = 12, units = "cm")
}

# Loop through all non-significant interactions
for (k in 1:nrow(nosignificant_interactions)) {
  # Get the data for the current non-significant result
  current_slope_data <- nosignificant_interactions$slope_data[[k]]
  current_task <- nosignificant_interactions$dataname[k]
  current_variable <- nosignificant_interactions$dependentvar[k]
  
  # Calculate the median and confidence intervals
  slope_summary <- current_slope_data %>%
    group_by(Age_year) %>%
    summarise(
      median_slope = median(slope, na.rm = TRUE),
      lower_95CI = quantile(slope, probs = 0.025, na.rm = TRUE),
      upper_95CI = quantile(slope, probs = 0.975, na.rm = TRUE)
    )
  
  # Plot the trend of slope changes with age
  slope_plot <- ggplot(slope_summary, aes(x = Age_year)) +
    geom_line(aes(y = median_slope), color = "black", size = 1) + 
    geom_ribbon(aes(ymin = lower_95CI, ymax = upper_95CI), alpha = 0.3, fill = "grey80") + 
    labs(title =paste0(modified_variable_name, " ~ ", modified_task_name), 
         x = "", 
         y = "Slope") +
    scale_x_continuous(name="", limits = c(11, 18),breaks = seq(11, 18, by = 1))+
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.5, color = "black"),  
      axis.text = element_text(size = 8.5,color = "black"), 
      axis.title = element_text(size = 8.5), 
      plot.title = element_text(size = 8.5, hjust = 0.5),
      panel.grid.major = element_blank(),    # 删除主网格线
      panel.grid.minor = element_blank(),    # 删除次网格线
      panel.border = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.5),  # 设置小线段颜色和大小
      axis.ticks.length = unit(0.2, "cm") 
    )
  print(slope_plot)
  # Save plot
  ggsave(paste0(FigureFolder, "/Nosignificant_slope_change_", current_task, "_", current_variable, ".pdf"), plot = slope_plot, width = 16, height = 12, units = "cm")
}