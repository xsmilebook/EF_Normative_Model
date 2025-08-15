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
# input directory
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- '/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/EF_Yunfu/interfileFolder_back12before'
  datapath_GNG <- '/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_all/interfileFolder'
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/EF_Yunfu/interfileFolder_back12before"
  functionFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/EF_Yunfu/code/Rcode/functions"
  resultFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/EF_Yunfu/results/corr"
}else{
  datapath <- 'D:/datasets/yunfu/raw_data'
  FigureFolder <- 'D:/datasets/yunfu/figures/fig1'
  interfileFolder <- "D:/datasets/yunfu/interfile_folder/temp"
  functionFolder <- "D:/code/EF_Normative_Model/functions"
  resultFolder <- "D:/datasets/yunfu/results/psy_corr"
}

# source functions
#source(paste0(functionFolder, "/gam_varyingcoefficients.R"))
source(paste0(functionFolder, "/gamcog_withsmoothvar_deviation.R"))
#source(paste0(functionFolder, "/ordinalcorr_new.R"))
# import dataset
#switch_data <- read_xlsx(paste0(datapath, '/Q_switch.xlsx'))
#head(switch_data)
PHQ_data <- read_xlsx(paste0(interfileFolder, '/PHQ.xlsx'))
GNGd_data <- read_rds(paste0(interfileFolder, '/GNGd_prime.deviations.rds'))
back1_data <- read_rds(paste0(interfileFolder, '/back1Acc.deviations.rds'))
back2_data <- read_rds(paste0(interfileFolder, '/back2Acc.deviations.rds'))

PHQ_data <- PHQ_data %>%
  select(用户ID, PHQ_y09) %>%  
  mutate(用户ID = as.character(用户ID))

GNGd_data <- GNGd_data %>% mutate(x__ID = as.character(x__ID))
back1_data <- back1_data %>% mutate(x__ID = as.character(x__ID))
back2_data <- back2_data %>% mutate(x__ID = as.character(x__ID))

GNGd_data <- inner_join(
  GNGd_data,
  PHQ_data,
  by = c("x__ID" = "用户ID")
)

back1_data <- inner_join(
  back1_data,
  PHQ_data,
  by = c("x__ID" = "用户ID")
)

back2_data <- inner_join(
  back2_data,
  PHQ_data,
  by = c("x__ID" = "用户ID")
)


## 1) set up variables
psyc_variables_continous <- c("PHQ_y09")

# EF vars
EFvars.set <- matrix(c("d_prime_deviationZ", "GNGd",
                       "Oneback_acc_deviationZ", "back1",
                       "Twoback_acc_deviationZ", "back2"), byrow=TRUE,ncol=2,dimnames=list(NULL,c("varname","dataname")))
EFvars.set <- as.data.frame(EFvars.set)
## 2) convert variables class & describe variables
GNGd_data[,psyc_variables_discrete] <- lapply(GNGd_data[,psyc_variables_discrete], as.factor)
GNGd_data[psyc_variables_continous] <- lapply(GNGd_data[psyc_variables_continous], as.numeric)

schooltab <- unique(GNGd_data$School)
GNGd_data$school_fac <- factor(GNGd_data$School, levels=schooltab, labels=paste0("school", 1:length(schooltab)))

back1_data[,psyc_variables_discrete] <- lapply(back1_data[,psyc_variables_discrete], as.factor)
back1_data[psyc_variables_continous] <- lapply(back1_data[psyc_variables_continous], as.numeric)
schooltab <- unique(back1_data$School)
back1_data$school_fac <- factor(back1_data$School, levels=schooltab, labels=paste0("school", 1:length(schooltab)))

back2_data[,psyc_variables_discrete] <- lapply(back2_data[,psyc_variables_discrete], as.factor)
back2_data[psyc_variables_continous] <- lapply(back2_data[psyc_variables_continous], as.numeric)
schooltab <- unique(back2_data$School)
back2_data$school_fac <- factor(back2_data$School, levels=schooltab, labels=paste0("school", 1:length(schooltab)))
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
write.csv(corr.result.period.df.con, paste0(resultFolder, "/corr_EF_psych_continuous.result.csv"), row.names = F)



corr.result.period.df.con <- read_csv(paste0(resultFolder, "/corr_EF_psych_continuous.result.csv"))
corr.result.period.df <- rbind(corr.result.period.df.con)
## plot
corr.result.period.df$correstimate <- as.numeric(corr.result.period.df$correstimate)
lwth <- min(corr.result.period.df$correstimate, na.rm = TRUE)
upth <- max(corr.result.period.df$correstimate, na.rm = TRUE)
y_levels <- c("PHQ_y09")  # 确保顺序
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
                     labels = c("PHQ_y09" = "Depressive Symptoms (PHQ-9)")) +
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
  ggsave(paste0(FigureFolder, "/", dataname, "/corr_continuouspsych_", EFvar.tmp, "_3periods.pdf"), plot = Fig, width = 14, height = 18, units = "cm")
}

# Combine all updated results
updated_results_df <- do.call(rbind, updated_results)

# Save as Excel or CSV file
write.csv(updated_results_df, file = paste0(FigureFolder, "/corr_results_with_fdr.csv"), row.names = FALSE)


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

## Data preparation
y_levels <- c("PHQ_y09")
combined_results$correstimate <- as.numeric(combined_results$correstimate)

# Custom labels for psychiatric variables
psy_labels <- c("PHQ_y09" = "PHQ_y09")

# Apply new labels
combined_results$parcel <- factor(combined_results$parcel,
                                  levels = y_levels,
                                  labels = psy_labels)
combined_results$Task <- factor(combined_results$dataname, levels = c("GNGd", "back1", "back2"),
                                labels = c("Go/No-go", "1-back", "2-back"))

# Manually define task colors
task_colors <- c("Go/No-go" = "#E4E9F3",
                 "1-back" = "#A4C5DF",
                 "2-back" = "#4980B5")

# Calculate the position for significance markers
combined_results$significance <- ifelse(combined_results$sig, "*", "")  # Mark with a star
combined_results$label_y <- ifelse(combined_results$correstimate > 0, 
                                   combined_results$correstimate + 0.01,  # Place the star above the bar top
                                   combined_results$correstimate - 0.02)  # Place the star below the bar bottom

## Plotting the figure
y_limits <- c(-0.15, 0.15)  # Symmetric y-limits
vline_positions <- seq(1.5, length(unique(combined_results$parcel)) - 0.5, by = 1)

Fig <- ggplot(data = combined_results, aes(x = parcel, y = correstimate, fill = Task, color = Task)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.5, size = 0) +  # Fill colors and borders
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.4) +  # Add horizontal line at y = 0
  geom_text(aes(label = significance, y = label_y), 
            position = position_dodge(width = 0.7), size = 5, color = "black") +  # Add significance stars
  scale_fill_manual(values = task_colors) +  # Manually set task fill colors
  scale_color_manual(values = task_colors) +  # Manually set task border colors
  scale_y_continuous(limits = y_limits, breaks = seq(-0.15, 0.15, by = 0.05), labels = scales::number_format()) +  # Set y-axis limits and ticks
  labs(title = "Correlation between Executive Functions and Mental Health",
       x = "",
       y = "Correlation Coefficient",
       color = "Tasks",
       fill = "Tasks") +  # Legend for tasks
  theme_minimal() +
  theme(axis.line.y = element_blank(),  # Remove y-axis
        axis.title = element_text(size = 8.5),
        axis.text.x = element_text(size = 8.5, hjust = 0.5, color = "black"),  # Ensure x-axis labels are centered
        axis.text.y = element_text(size = 8.5, color = "black"),  
        plot.title = element_text(size = 8.5, hjust = 0.5),
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8.5),
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"), 
        panel.grid.major.y = element_line(color = "grey90", linetype = "solid", size = 0.2),  # Add light lines as horizontal grid
        panel.grid.major.x = element_blank(),  # Remove x-axis grid lines
        panel.grid.minor = element_blank()) + # Remove minor grid lines
  annotate("segment", x = vline_positions, xend = vline_positions, y = -0.005, yend = 0, color = "black", size = 0.4)


if (length(unique(combined_results$parcel)) > 1) {
  vline_positions <- seq(1.5, length(unique(combined_results$parcel)) - 0.5, by = 1)
  Fig <- Fig + annotate("segment", x = vline_positions, xend = vline_positions, y = -0.005, yend = 0, color = "black", size = 0.4)
}

print(Fig)
ggsave(paste0(FigureFolder, "/correlation_barplot_sdq5_withSDQsum.pdf"), plot = Fig, width = 18, height = 10, units = "cm")

print(Fig)
ggsave(paste0(FigureFolder, "/correlation_barplot_sdq5_withSDQsum.pdf"), plot = Fig, width = 18, height = 10, units = "cm")