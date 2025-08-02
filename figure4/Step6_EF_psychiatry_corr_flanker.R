library(ggplot2)
library(tidyverse)
library(mgcv)
library(readxl)
library(parallel)
library(gamlss)
library(scales)
library(tableone)
library(openxlsx)
library(writexl)
library(dplyr)
# input directory
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- '/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/ABCD_verify/interfileFolder'
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/ABCD_verify/interfileFolder"
  functionFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/ABCD_verify/Rcode/functions_corr"
  resultFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/ABCD_verify/interfileFolder/corr"
}else{
  datapath <- '/Users/tanlirou/Documents/YF_EF_psy/EF_psy/ABCD/correlation/data'
  FigureFolder <- '/Users/tanlirou/Documents/YF_EF_psy/EF_psy/ABCD/correlation/？'
  interfileFolder <- "/Users/tanlirou/Documents/EF_yunfu_check/ABCD/interfileFolder"
  functionFolder <- "/Users/tanlirou/Documents/YF_EF_psy/EF_psy/ABCD/correlation/funtions"
  resultFolder <- "/Users/tanlirou/Documents/YF_EF_psy/EF_psy/ABCD/correlation/results_test"
}

# source functions
source(paste0(functionFolder, "/gamm_factor_interaction_new.R"))
# import dataset
#switch_data <- read_xlsx(paste0(datapath, '/Q_switch.xlsx'))
#head(switch_data)
# Flanker_data <- read_csv(paste0(datapath, '/Flanker.deviations.csv'))
# cbcl_data <- read_xlsx(paste0(datapath,'/mh_p_cbcl_r.xlsx'))
# # load data
# cbcl_data <- cbcl_data %>%
#   select(ID, cbcl_scr_syn_totprob_r, cbcl_scr_syn_social_r, cbcl_scr_syn_attention_r, 
#          cbcl_scr_syn_internal_r, cbcl_scr_syn_external_r)
# 
# merged_data <- left_join(Flanker_data, cbcl_data, by = "ID")
# write.csv(merged_data, paste0(datapath,'/Flanker.deviations_addr.csv'), row.names = FALSE)

Flanker_data <- read_csv(paste0(datapath, '/Flanker.deviations_addr.csv'))

## 1) set up variables
psyc_variables_continous <- c("cbcl_scr_syn_internal_r","cbcl_scr_syn_social_r",
                              "cbcl_scr_syn_external_r","cbcl_scr_syn_attention_r")
# EF vars
EFvar <- "nihtbx_flanker_uncorrected_deviationZ"
## 2) convert variables class & describe variables
Flanker_data[, c(psyc_variables_continous, EFvar)] <- lapply(Flanker_data[, c(psyc_variables_continous, EFvar)], as.numeric)
site_id_ltab <- unique(Flanker_data$site_id_l)
Flanker_data$site_id_l_fac <- factor(Flanker_data$site_id_l, levels=site_id_ltab, labels=paste0("site_id_l", 1:length(site_id_ltab)))

###
describe_tab_Flanker <- CreateTableOne(c(EFvar, psyc_variables_continous, "Age_year", "Sex"), 
                                       data = Flanker_data, testNonNormal = TRUE)
describe_tab_Flanker.continous <- as.data.frame(describe_tab_Flanker$ContTable[["Overall"]])
# save out
write.xlsx(list(Flanker_con=describe_tab_Flanker.continous), paste0(resultFolder, "/description_interest_vars.xlsx"), rowNames=T)

## 3) Correlations in separate age periods
# continuos variables
corr.result <- list()
for (psyvar.tmp in psyc_variables_continous) {
  dependentvar <- psyvar.tmp
  interest.indep.var <- EFvar
  covariates <- "Sex"
  smoothvar <- "Age_year"
  
  # Perform analysis using gamm.smooth.predict.interaction
  result.tmp <- gamm.smooth.predict.interaction(
    dependentvar = dependentvar,
    dataname = "Flanker_data",
    smoothvar = smoothvar,
    interest.indep.var = interest.indep.var,
    covariates = covariates
  )
  
  result.tmp <- as.data.frame(result.tmp)
  corr.result[[psyvar.tmp]] <- result.tmp
}

# Combine results
corr.result.df <- do.call(rbind, corr.result)
write.csv(corr.result.df, paste0(resultFolder, "/corr_EF_psych_continuous.result_ABCD.csv"), row.names = FALSE)

corr.result.df <- read_csv(paste0(resultFolder, "/corr_EF_psych_continuous.result_ABCD.csv"))
corr.result.df$correstimate <- as.numeric(corr.result.df$correstimate)
corr.result.df$anovap.fdr <- p.adjust(corr.result.df$boots.pvalues, method = "fdr")
corr.result.df$sig <- (corr.result.df$anovap.fdr < 0.05)
corr.result.df$significance <- ifelse(corr.result.df$sig, "*", "")

write.csv(corr.result.df, paste0(resultFolder, "/corr_EF_psych_continuous.result_ABCD_withfdr.csv"), row.names = FALSE)

# Define labels and limits
psy_labels <- c("cbcl_scr_syn_internal_r" = "Internalizing Symptoms", 
                "cbcl_scr_syn_social_r" = "Social",
                "cbcl_scr_syn_external_r" = "Externalizing Symptoms",
                "cbcl_scr_syn_attention_r" = "Attention")

corr.result.df$parcel <- factor(corr.result.df$parcel,
                                levels = psyc_variables_continous,
                                labels = psy_labels)

y_limits <- c(-0.125, 0.05)
corr.result.df$label_y <- ifelse(
  corr.result.df$correstimate >= 0, 
  corr.result.df$correstimate + 0.01, 
  corr.result.df$correstimate - 0.02
)

vline_positions <- seq(1.5, length(unique(corr.result.df$parcel)) - 0.5, by = 1)

Fig <- ggplot(data = corr.result.df, aes(x = parcel, y = correstimate, color = Task, fill = Task)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.5, size = 1, 
           color = NA, show.legend = T) +  
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.4) +  # Add a horizontal line at y = 0 with a light gray color as a reference line
  geom_text(aes(label = significance, y = label_y), 
            position = position_dodge(width = 0.5), size = 5, color = "black") + 
  scale_fill_manual(values = task_colors) +  
  scale_color_manual(values = task_colors) +  
  scale_y_continuous(limits = y_limits, breaks = seq(-0.1, 0.05, by = 0.05), labels = scales::number_format()) + 
  labs(title = "Correlation between Flanker and Mental Health",x = "", y = "correlation coefficient",color = "Tasks") +  
  theme_minimal() +
  theme(axis.line.y = element_blank(),  
        axis.title = element_text(size = 8.5),
        axis.text.x = element_text(size = 8.5, hjust = 0.5,color = "black"),  
        axis.text.y = element_text(size = 8.5,color = "black"),  
        plot.title = element_text(size = 8.5, hjust = 0.5),
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8.5),
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"), 
        panel.grid.major.y = element_line(color = "gray90", linetype = "solid", size = 0.2),  
        panel.grid.major.x = element_blank(),  
        panel.grid.minor = element_blank()) +
  annotate("segment", x = vline_positions, xend = vline_positions, y = -0.005, yend = 0, color = "black", size = 0.4)

print(Fig)
ggsave(paste0(FigureFolder, "/correlation_barplot.pdf"), plot = Fig, width = 10, height = 7, units = "cm")
