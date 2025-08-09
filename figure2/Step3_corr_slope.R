
rm(list=ls())

library(tidyverse)
library(mgcv)
library(psych)
library(gamlss)
library(scales)
library(tableone)
library(openxlsx)
library(ggplot2)
library(patchwork) 

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
  interfileFolder <- "D:/datasets/yunfu/interfile_folder/Normative_Model"
  functionFolder <- "D:/code/EF_Normative_Model/functions"
  resultFolder <- "D:/datasets/yunfu/results/gamlss/Normative_Model"
}

source(paste0(functionFolder, "/gamcog_withsmoothvar_deviation.R"))

#head(switch_data)
GNGd_data <- read_rds(file.path(interfileFolder, "GNGd_prime", 'GNGd_prime.deviations.rds'))
back1_data <- read_rds(file.path(interfileFolder, "1-back", 'back1Acc.deviations.rds'))
back2_data <- read_rds(file.path(interfileFolder, "2-back", 'back2Acc.deviations.rds'))

print(paste0(nrow(GNGd_data)," ", nrow(back1_data)," ", nrow(back2_data)))

GNGd_data <- GNGd_data %>%
  filter(
    d_prime_deviationZ > (mean(d_prime_deviationZ) - 3 * sd(d_prime_deviationZ)),
    d_prime_deviationZ < (mean(d_prime_deviationZ) + 3 * sd(d_prime_deviationZ))
  )
back1_data <- back1_data %>%
  filter(
    Oneback_acc_deviationZ > (mean(Oneback_acc_deviationZ) - 3 * sd(Oneback_acc_deviationZ)),
    Oneback_acc_deviationZ < (mean(Oneback_acc_deviationZ) + 3 * sd(Oneback_acc_deviationZ))
  )
back2_data <- back2_data %>%
  filter(
    Twoback_acc_deviationZ > (mean(Twoback_acc_deviationZ) - 3 * sd(Twoback_acc_deviationZ)),
    Twoback_acc_deviationZ < (mean(Twoback_acc_deviationZ) + 3 * sd(Twoback_acc_deviationZ))
  )

print(paste0(nrow(GNGd_data)," ", nrow(back1_data)," ", nrow(back2_data)))

psyc_variables_continous <- c("SDQ_PB_sum", "SDQ_H_sum", "SDQ_CP_sum", "SDQ_PP_sum", "SDQ_ES_sum")
psyc_stats <- data.frame(
  variable = character(),
  original_n = numeric(),
  filtered_n = numeric(),
  removed_n = numeric(),
  stringsAsFactors = FALSE
)

for (psyc_item in psyc_variables_continous){
  if(psyc_item %in% names(GNGd_data)) {
    original_n <- nrow(GNGd_data)
    mean_val <- mean(GNGd_data[[psyc_item]], na.rm = TRUE)
    sd_val <- sd(GNGd_data[[psyc_item]], na.rm = TRUE)
    
    GNGd_data_filtered <- GNGd_data %>%
      filter(
        !!sym(psyc_item) > (mean_val - 3 * sd_val),
        !!sym(psyc_item) < (mean_val + 3 * sd_val)
      )
    
    filtered_n <- nrow(GNGd_data_filtered)
    removed_n <- original_n - filtered_n
    
    psyc_stats <- rbind(psyc_stats, data.frame(
      variable = paste0(psyc_item, "_GNGd"),
      original_n = original_n,
      filtered_n = filtered_n,
      removed_n = removed_n,
      stringsAsFactors = FALSE
    ))
  }
  
  if(psyc_item %in% names(back1_data)) {
    original_n <- nrow(back1_data)
    mean_val <- mean(back1_data[[psyc_item]], na.rm = TRUE)
    sd_val <- sd(back1_data[[psyc_item]], na.rm = TRUE)
    
    back1_data_filtered <- back1_data %>%
      filter(
        !!sym(psyc_item) > (mean_val - 3 * sd_val),
        !!sym(psyc_item) < (mean_val + 3 * sd_val)
      )
    
    filtered_n <- nrow(back1_data_filtered)
    removed_n <- original_n - filtered_n
    
    psyc_stats <- rbind(psyc_stats, data.frame(
      variable = paste0(psyc_item, "_back1"),
      original_n = original_n,
      filtered_n = filtered_n,
      removed_n = removed_n,
      stringsAsFactors = FALSE
    ))
  }
  
  if(psyc_item %in% names(back2_data)) {
    original_n <- nrow(back2_data)
    mean_val <- mean(back2_data[[psyc_item]], na.rm = TRUE)
    sd_val <- sd(back2_data[[psyc_item]], na.rm = TRUE)
    
    back2_data_filtered <- back2_data %>%
      filter(
        !!sym(psyc_item) > (mean_val - 3 * sd_val),
        !!sym(psyc_item) < (mean_val + 3 * sd_val)
      )
    
    filtered_n <- nrow(back2_data_filtered)
    removed_n <- original_n - filtered_n
    
    psyc_stats <- rbind(psyc_stats, data.frame(
      variable = paste0(psyc_item, "_back2"),
      original_n = original_n,
      filtered_n = filtered_n,
      removed_n = removed_n,
      stringsAsFactors = FALSE
    ))
  }
}



## 1) set up variables
psyc_variables_continous <- c("SDQ_PB_sum", "SDQ_H_sum", "SDQ_CP_sum", "SDQ_PP_sum", "SDQ_ES_sum")
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

back1_data[,psyc_variables_discrete] <- lapply(back1_data[,psyc_variables_discrete], as.factor)
back1_data[,psyc_variables_continous] <- lapply(back1_data[,psyc_variables_continous], as.numeric)
schooltab <- unique(back1_data$School)
back1_data$school_fac <- factor(back1_data$School, levels=schooltab, labels=paste0("school", 1:length(schooltab)))

back2_data[,psyc_variables_discrete] <- lapply(back2_data[,psyc_variables_discrete], as.factor)
back2_data[,psyc_variables_continous] <- lapply(back2_data[,psyc_variables_continous], as.numeric)
schooltab <- unique(back2_data$School)
back2_data$school_fac <- factor(back2_data$School, levels=schooltab, labels=paste0("school", 1:length(schooltab)))
###3# 对 SDQ 连续变量进行标准化：mean = 0, sd = 1
GNGd_data[, paste0(psyc_variables_continous, "_z")] <- scale(GNGd_data[, psyc_variables_continous])
back1_data[, paste0(psyc_variables_continous, "_z")] <- scale(back1_data[, psyc_variables_continous])
back2_data[, paste0(psyc_variables_continous, "_z")] <- scale(back2_data[, psyc_variables_continous])

psyc_variables_continous <- paste0(psyc_variables_continous, "_z")

# describe
# GNGd
describe_tab_GNGd <- CreateTableOne(c(EFvars.set$varname[EFvars.set$dataname=="GNGd"], psyc_variables_continous, psyc_variables_discrete, "Age_year", "Sex"), data=GNGd_data, factorVars = psyc_variables_discrete, testNonNormal=T)
describe_tab_GNGd.continous <- as.data.frame(describe_tab_GNGd$ContTable[["Overall"]])
describe_tab_GNGd.discrete <- do.call(rbind, lapply(describe_tab_GNGd$CatTable[["Overall"]], as.data.frame))
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
  data.tmp.EF <- data.tmp
  
  for (period in c("EF")){
    dataname <- paste0("data.tmp.", period)
    corr.result <- list()
    
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
      smoothvar <- "Age_year"
      covariates <- "Sex"
      knots=3
      
      result.tmp <- gam.fit.Independent.var(dependentvar, dataname,  smoothvar, interest.indep.var, covariates, stats_only = T)
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
write.csv(corr.result.period.df.con, paste0(resultFolder, "/corr_EF_psych_continuous_result_slope.csv"), row.names = F)

corr.result.period.df.con <- read_csv("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/correlation_change/results250515/corr_EF_psych_continuous_result_slope.csv")
#corr.result.period.df.con <- read_csv(paste0(resultFolder, "/corr_EF_psych_continuous.result.csv"))
corr.result.period.df <- rbind(corr.result.period.df.con)
## plot
corr.result.period.df$correstimate <- as.numeric(corr.result.period.df$correstimate)
lwth <- min(corr.result.period.df$correstimate, na.rm = TRUE)
upth <- max(corr.result.period.df$correstimate, na.rm = TRUE)
y_levels <- c("SDQ_PB_sum", "SDQ_H_sum", "SDQ_CP_sum", "SDQ_PP_sum",  "SDQ_ES_sum")  # 确保顺序
# Initialize the result list
updated_results <- list()

# Loop through each EFvar
for (i in 1:nrow(EFvars.set)) {
  EFvar.tmp <- EFvars.set$varname[i]
  dataname <- EFvars.set$dataname[i]
  
  # Extract data for the corresponding variable
  corr.result.tmp <- corr.result.period.df[which(corr.result.period.df$interest.indep.var == EFvar.tmp & corr.result.period.df$dataname == dataname), ]
  corr.result.tmp$anova.pvalues <- as.numeric(corr.result.tmp$anova.pvalues)
  
  # Calculate FDR
  corr.result.tmp$anovap.fdr <- p.adjust(corr.result.tmp$anova.pvalues, method = "fdr")
  corr.result.tmp$sig <- (corr.result.tmp$anovap.fdr < 0.05)
  
  # Update the main data frame
  updated_results[[i]] <- corr.result.tmp
  
  # Set y-axis order
  corr.result.tmp$parcel <- factor(corr.result.tmp$parcel, levels = y_levels)
  
  # Plot the figure
  Fig <- ggplot() +
    geom_tile(data = corr.result.tmp, aes(x = period, y = parcel, fill = correstimate), color = "white") +
    geom_text(data = corr.result.tmp[corr.result.tmp$sig == TRUE, ], aes(x = period, y = parcel, label = "*"), vjust = 0.7, hjust = 0.5, size = 6) +
    scale_fill_distiller(type = "seq", palette = "RdBu", limits = c(lwth, upth), direction = -1) +
    scale_y_discrete(limits = y_levels,
                     labels = c("SDQ_PB_sum" = "Prosocial Behavior","SDQ_H_sum" = "Hyperactivity",
                                "SDQ_CP_sum" = "Conduct Problems","SDQ_PP_sum" = "Peer Problems",
                                "SDQ_ES_sum" = "Emotional Symptoms")) +
    labs(title = paste0("Correlation between Executive functions and ", EFvar.tmp, " in ", dataname),
         x = "Periods of adolescence", y = "Psychiatric scores") +
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
  #ggsave(paste0(FigureFolder, "/", dataname, "/corr_continuouspsych_", EFvar.tmp, "_3periods.pdf"), plot = Fig, width = 14, height = 18, units = "cm")
}

# Combine all updated results
updated_results_df <- do.call(rbind, updated_results)

# Save as Excel or CSV file
write.csv(updated_results_df, file = paste0("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/correlation_change/results250515/corr_results_with_fdr.csv"), row.names = FALSE)


combined_results <- data.frame()
for (i in 1:nrow(EFvars.set)) {
  EFvar.tmp <- EFvars.set$varname[i]
  dataname <- EFvars.set$dataname[i]
  corr.result.tmp <- corr.result.period.df[which(corr.result.period.df$interest.indep.var == EFvar.tmp & corr.result.period.df$dataname == dataname), ]
  corr.result.tmp$anova.pvalues <- as.numeric(corr.result.tmp$anova.pvalues)
  corr.result.tmp$correstimate <- as.numeric(corr.result.tmp$correstimate)
  corr.result.tmp$anovap.fdr <- p.adjust(corr.result.tmp$anova.pvalues, method = "fdr")
  corr.result.tmp$sig <- (corr.result.tmp$anovap.fdr < 0.05)
  corr.result.tmp$period <- factor(corr.result.tmp$period, levels = c("EF"))
  combined_results <- rbind(combined_results, corr.result.tmp)
}



combined_results <- read.csv("/Users/tanlirou/Documents/YF_EF_psy/EF_psy/correlation_change/resutls250515/corr_results_with_fdr.csv")
## Data preparation
y_levels <- c("SDQ_ES_sum_z", "SDQ_PP_sum_z", "SDQ_CP_sum_z", "SDQ_H_sum_z", "SDQ_PB_sum_z")
combined_results$correstimate <- as.numeric(combined_results$correstimate)

# Custom labels for psychiatric variables
psy_labels <- c("SDQ_ES_sum_z" = "Emotional
Symptoms",
                "SDQ_PP_sum_z" = "Peer
Problems", 
                "SDQ_CP_sum_z" = "Conduct
Problems",
                "SDQ_H_sum_z" = "Hyperactivity", 
                "SDQ_PB_sum_z" = "Prosocial
Behavior")

# Apply new labels
combined_results$parcel <- factor(combined_results$parcel,
                                  levels = y_levels,
                                  labels = psy_labels)
combined_results$Task <- factor(combined_results$dataname, levels = c("GNGd", "back1", "back2"),
                                labels = c("Go/No-go", "1-back", "2-back"))

# Manually define task colors
task_colors <- c("Go/No-go" = "#E5E9F2",
                 "1-back" = "#A4C5DF",
                 "2-back" = "#4980B5")

# Calculate the position for significance markers
combined_results$significance <- ifelse(combined_results$sig, "*", "")  # Mark with a star
combined_results$label_y <- ifelse(combined_results$slope > 0, 
                                   combined_results$slope + 0.01,  # Place the star above the bar top
                                   combined_results$slope - 0.02)  # Place the star below the bar bottom

## Plotting the figure
y_limits <- c(-0.15, 0.12)  # Symmetric y-limits
vline_positions <- seq(1.5, length(unique(combined_results$parcel)) - 0.5, by = 1)

# Fig <- ggplot(data = combined_results, aes(x = parcel, y = correstimate, fill = Task, color = Task)) +
#   geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.5, size = 0) +  # Fill colors and borders
#   geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.4) +  # Add horizontal line at y = 0
#   geom_text(aes(label = significance, y = label_y), 
#             position = position_dodge(width = 0.7), size = 5, color = "black") +  # Add significance stars
#   scale_fill_manual(values = task_colors) +  # Manually set task fill colors
#   scale_color_manual(values = task_colors) +  # Manually set task border colors
#   scale_y_continuous(limits = y_limits, breaks = seq(-0.15, 0.15, by = 0.05), labels = scales::number_format()) +  # Set y-axis limits and ticks
#   labs(title = "Correlation between Executive Functions and Mental Health",
#        x = "",
#        y = "Correlation Coefficient",
#        color = "Tasks",
#        fill = "Tasks") +  # Legend for tasks
#   theme_minimal() +
#   theme(axis.line.y = element_blank(),  # Remove y-axis
#         axis.title = element_text(size = 8.5),
#         axis.text.x = element_text(size = 8.5, hjust = 0.5, color = "black"),  # Ensure x-axis labels are centered
#         axis.text.y = element_text(size = 8.5, color = "black"),  
#         plot.title = element_text(size = 8.5, hjust = 0.5),
#         legend.title = element_text(size = 8.5),
#         legend.text = element_text(size = 8.5),
#         legend.position = "bottom",
#         legend.key.size = unit(0.4, "cm"), 
#         panel.grid.major.y = element_line(color = "grey90", linetype = "solid", size = 0.2),  # Add light lines as horizontal grid
#         panel.grid.major.x = element_blank(),  # Remove x-axis grid lines
#         panel.grid.minor = element_blank()) + # Remove minor grid lines
#   annotate("segment", x = vline_positions, xend = vline_positions, y = -0.005, yend = 0, color = "black", size = 0.4)



Fig <- ggplot(data = combined_results, aes(x = parcel, y = slope, fill = Task, color = Task)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.75, size = 0) +  # Fill colors and borders
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.25) +  # Add horizontal line at y = 0
  geom_text(aes(label = significance, y = label_y), 
            position = position_dodge(width = 0.8), size = 5, color = "black") +  # Add significance stars
  scale_fill_manual(values = task_colors) +  # Manually set task fill colors
  scale_color_manual(values = task_colors) +  # Manually set task border colors
  scale_y_continuous(limits = y_limits, breaks = seq(-0.15, 0.1, by = 0.05), labels = scales::number_format()) +  # Set y-axis limits and ticks
  labs(title = "",
       x = "",
       y = "β",
       color = "Tasks",
       fill = "Tasks") +  # Legend for tasks
  theme_minimal() +
  theme(axis.line.x = element_line(color = "black", size = 0.25),  # Add x-axis line
        axis.line.y = element_line(color = "black", size = 0.25),  # Add y-axis line
        axis.title = element_text(size = 9),
        axis.text.x = element_text(size = 9, hjust = 0.5, color = "black"),  # Ensure x-axis labels are centered
        axis.text.y = element_text(size = 9, color = "black"),  
        axis.ticks.x = element_line(color = "black", size = 0.25),
        axis.ticks.y = element_line(color = "black", size = 0.25),
        axis.ticks.length = unit(0.05, "cm"),
        plot.title = element_text(size = 9, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        #legend.position = c(0.95, 0.15),
        legend.position =c(0.25, 0.9),
        legend.direction = "horizontal" ,
        legend.margin = margin(t = -10),
        legend.key.size = unit(0.3, "cm"), 
        panel.grid.major.y = element_blank(),  # Remove x-axis grid lines
        panel.grid.major.x = element_blank(),  # Remove x-axis grid lines
        panel.grid.minor = element_blank()) 

print(Fig)
ggsave(paste0(FigureFolder, "/correlation_beta_sdq5.pdf"), plot = Fig, width = 15, height = 8, units = "cm")
