# import packages
library(tidyverse)
library(ggplot2)
library(patchwork)

# 2. setting path
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  resultFolder_individual <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/datasets/yunfu/results/EF_psy"
  resultFolder_final <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/EF_Normative_Model/anova/results"
  FigureFolder_final <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/EF_Normative_Model/anova/figures"
} else {
  resultFolder_individual <- "/Users/tanlirou/Documents/YF_EF_psy/EF_psy/code_pure_202507/results/Figure2_corr_delete_Z/individual_results"
  resultFolder_final <- "/Users/tanlirou/Documents/YF_EF_psy/EF_psy/code_pure_202507/results/Figure2_corr_delete_Z"
  FigureFolder_final <- '/Users/tanlirou/Documents/YF_EF_psy/EF_psy/code_pure_202507/results/Figure2_corr_delete_Z'
}

# 3. merge results
cat("read .csv file from:", resultFolder_individual, "\n")
csv_files <- list.files(path = resultFolder_individual, pattern = "^corr_GNGd_.*\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) {
  stop("error: failed to find any file in specific folder.")
}

corr.result.df <- map_dfr(csv_files, read.csv)
cat(paste("merge results: ", nrow(corr.result.df), ".\n"))

# 4. fdr
corr.result.df <- corr.result.df %>%
  mutate(
    correstimate = as.numeric(correstimate),
    anova.pvalues = as.numeric(anova.pvalues),
    corrp = as.numeric(corrp),
    anovap.bonf = p.adjust(anova.pvalues, method = "bonferroni"),
    sig = anovap.bonf < 0.05
  )

# 5. heatmap
y_levels <- c("SDQ_PB_sum_z", "SDQ_H_sum_z", "SDQ_CP_sum_z", "SDQ_PP_sum_z", "SDQ_ES_sum_z")
y_labels <- c("SDQ_ES_sum_z" = "Emotional Symptoms", "SDQ_PP_sum_z" = "Peer Problems", 
              "SDQ_CP_sum_z" = "Conduct Problems", "SDQ_H_sum_z" = "Hyperactivity", 
              "SDQ_PB_sum_z" = "Prosocial Behavior")

corr.result.df$parcel <- factor(corr.result.df$parcel, levels = rev(y_levels))

Fig <- ggplot(corr.result.df, aes(x = period, y = parcel, fill = correstimate)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(data = filter(corr.result.df, sig == TRUE), 
            aes(label = "*"), vjust = 0.75, size = 8, color = "black") +
  scale_fill_distiller(palette = "RdBu", direction = -1, name = "Estimate") +
  scale_y_discrete(labels = y_labels) +
  labs(
    title = "Correlation between GNGd and Psychiatric Scores",
    x = NULL, y = "Psychiatric Scores"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    panel.grid = element_blank()
  )

# 6. save results
final_plot_path <- paste0(FigureFolder_final, "/GNGd_correlation_heatmap_final.pdf")
ggsave(final_plot_path, plot = Fig, width = 8, height = 7)
cat("heat map was saved to:", final_plot_path, "\n")

final_csv_path <- paste0(resultFolder_final, "/GNGd_corr_results_with_bonf_aggregated.csv")
write.csv(corr.result.df, file = final_csv_path, row.names = FALSE)
cat("Bonferroni correction results: ", final_csv_path, "\n")
