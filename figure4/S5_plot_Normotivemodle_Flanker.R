rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(mgcv)
library(openxlsx)
library(parallel)
library(gamlss)
library(scales)
library(psych)
library(reshape)
library(cowplot)
library(writexl)
library(extrafont)
library(patchwork)
# input directory
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- '/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/data0715'
  demopath <- '/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/data0715'
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/interfileFolder_Acc"
  functionFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/Rcode_EFnorms/functions"
  resultFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/results_Acc"
}else{
  datapath <- '/Users/tanlirou/Documents/YF_EF_psy/EF_psy/ABCD/interfileFolder'
  FigureFolder <- '/Users/tanlirou/Documents/YF_EF_psy/EF_psy/ABCD/figureFolder'
  interfileFolder <- '/Users/tanlirou/Documents/YF_EF_psy/EF_psy/ABCD/interfileFolder'
  functionFolder <- '/Users/tanlirou/Documents/YF_EF_psy/EF_psy/ABCD/code/functions'
  resultFolder <- "/Users/tanlirou/Documents/YF_EF_psy/EF_psy/ABCD/results"
}
sumFlanker_deviation <- readRDS(paste0(interfileFolder,"/Flanker.deviations.rds"))
modFlanker.set1.sum <- readRDS(paste0(interfileFolder, "/GAMLSS.Flankerset1.sum.rds"))
modFlanker.set2.sum <- readRDS(paste0(interfileFolder, "/GAMLSS.Flankerset2.sum.rds"))
Flanker_data1 <- read_csv(paste0(interfileFolder, "/Flanker_data.set1.csv"))
Flanker_data2 <- read_csv(paste0(interfileFolder, "/Flanker_data.set2.csv"))
derivative_Flanker <- read_csv(paste0(interfileFolder, "/derivative_summary_Flanker.csv"))

sumFlanker_deviation$Sex <- as.factor(sumFlanker_deviation$Sex)
fillcolor <- c("#F5B7BF", "#91ACE0")
# bar plot
Flanker_plot <- ggplot(data = sumFlanker_deviation, aes(x = Age_year, y = after_stat(count), fill = Sex)) +
  geom_histogram(binwidth = 1, color = "black", boundary = 11, position = "stack", linewidth = 0.5) + 
  labs(x = "Age (year)", y = NULL, title = paste0("Go/No-go, N=", nrow(sumFlanker_deviation))) +
  scale_fill_manual(values = fillcolor, name = "Sex") +
  scale_x_continuous(breaks = seq(8, max(sumFlanker_deviation$Age_year, na.rm = TRUE), by = 1)) + 
  theme_classic() +
  theme(
    aspect.ratio = 0.8,
    plot.title = element_text(color = "black", size = 8.5, hjust = 0.5),
    axis.title = element_text(color = "black", size = 8.5),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    axis.text = element_text(color = "black", size = 8.5)
  )
print(Flanker_plot)
# save as PDF 
ggsave(filename = paste0(FigureFolder, "/flankernumber.pdf"), plot = Flanker_plot, width = 9, height=7, units="cm")
#####
###plot centile
# hist(sumFlanker_deviation$nihtbx_flanker_uncorrected, breaks = 30, main = "Histogram of nihtbx_flanker_uncorrected", xlab = "nihtbx_flanker_uncorrected")
# hist(sumFlanker_deviation$nihtbx_flanker_uncorrected_centile, main = "Centile Distribution", xlab = "nihtbx_flanker_uncorrected_centile", ylab = "Frequency", col = "blue")
# 
# 
# ggplot(sumFlanker_deviation, aes(x = nihtbx_flanker_uncorrected_centile)) +
#   geom_histogram(binwidth = 0.1, fill = "blue", color = "black") +
#   theme_minimal()



## 1. deviation proportion
sumFlanker_deviation[["nihtbx_flanker_uncorrected_deviationZ_extreme_p"]] <- sumFlanker_deviation[["nihtbx_flanker_uncorrected_deviationZ"]] > 1.96
sumFlanker_deviation[["nihtbx_flanker_uncorrected_deviationZ_extreme_n"]] <- sumFlanker_deviation[["nihtbx_flanker_uncorrected_deviationZ"]] < -1.96
proportion_Flanker_extreme_p <- mean(sumFlanker_deviation[["nihtbx_flanker_uncorrected_deviationZ_extreme_p"]], na.rm = TRUE)
print(paste("Proportion of positive extreme deviations is", round(proportion_Flanker_extreme_p * 100, 2), "%"))

proportion_Flanker_extreme_n <- mean(sumFlanker_deviation[["nihtbx_flanker_uncorrected_deviationZ_extreme_n"]], na.rm = TRUE)
print(paste("Proportion of negative extreme deviations is", round(proportion_Flanker_extreme_n * 100, 2), "%"))

saveRDS(sumFlanker_deviation, paste0(interfileFolder, "/processed_sum_Flanker_deviation.rds"))

# Plotting extreme points
pFlankerextreme <- ggplot(sumFlanker_deviation, aes(x = Age_year, y = nihtbx_flanker_uncorrected)) +
  geom_point(aes(color = nihtbx_flanker_uncorrected_deviationZ_extreme_p | nihtbx_flanker_uncorrected_deviationZ_extreme_n), alpha = 0.6) +
  scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "blue"), name = "Extreme Deviation") +
  labs(title = "Deviation Z with Extreme Points", x = "Age", y = "Deviation Z") +
  theme_minimal()
# Show plot
print(pFlankerextreme)

# plotting extreme line
ggplot(sumFlanker_deviation, aes(x = nihtbx_flanker_uncorrected)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 1.96, linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = -1.96, linetype = "dashed", color = "red", size = 1) +
  labs(title = "Distribution of Flanker Deviation Z-Scores",
       x = "Flanker Deviation Z-Scores",
       y = "Count") +
  theme_minimal()


###2 plot normitive trajectory
# Calculate centiles
n_points <- 1000
centiles.tmp.set1 <- modFlanker.set1.sum$centiles_strat
centiles.tmp.set2 <- modFlanker.set2.sum$centiles_strat

centiles.F.tmp <- (centiles.tmp.set1[[2]] + centiles.tmp.set2[[2]]) / 2
centiles.M.tmp <- (centiles.tmp.set1[[1]] + centiles.tmp.set2[[1]]) / 2
Centiles <- (centiles.F.tmp + centiles.M.tmp) / 2

X  <- seq(min(sumFlanker_deviation$Age_year), max(sumFlanker_deviation$Age_year), length.out=1000)
sumFlanker_deviation$Sex <- factor(sumFlanker_deviation$Sex, levels = c("F", "M"))
# Plotting
p_main <- ggplot() +
  #geom_hline(yintercept = seq(70, 120, by = 10), linetype = "solid", color = "gray70", size = 0.5) + 
  #geom_point(data=sumFlanker_deviation, aes(x=Age_year, y=nihtbx_flanker_uncorrected), color= "#fac858", size = 3, alpha=0.25, show.legend=FALSE, stroke = 0) +
  geom_jitter(data=sumFlanker_deviation, aes(x=Age_year, y=nihtbx_flanker_uncorrected), 
              color= "#A4C5DF", size = 2, alpha=0.25,show.legend=FALSE, stroke = 0, width = 0, height = 1) +
  geom_line(aes(x=X, y=Centiles[3,]), linetype="dashed", color = "gray10", size = 0.8) +
  geom_line(aes(x=X, y=Centiles[4,]), linetype="dashed", color = "gray10", size = 0.8) +
  geom_line(aes(x=X, y=Centiles[5,]), linetype="solid", color = "black", size = 0.8) +
  geom_line(aes(x=X, y=Centiles[6,]), linetype="dashed", color = "gray10", size = 0.8) +
  geom_line(aes(x=X, y=Centiles[7,]), linetype="dashed", color = "gray10", size = 0.8) +
  scale_y_continuous(name="Score",limits = c(60, 120),breaks = seq(60, 120, by = 15))+
  scale_x_continuous(name="", limits = c(8.8, 16),breaks = seq(9, 16, by = 1))+
  labs(x="", title="Flanker") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 8.5, hjust = 0.5, color = "black"),
    axis.text = element_text(colour = "black",size=8.5), 
    axis.title = element_text(size=8.5,face = "plain"),
    axis.line = element_line(colour = "black", size = 0.5),
    axis.ticks = element_line(colour = "black", size = 0.5), 
    axis.ticks.length = unit(0.2, "cm"),
    axis.text.y = element_text(colour = "black",size = 8.5),
    plot.background=element_rect(fill="white",color = NA),
    legend.position = "none" ,
    panel.grid = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm"))
x_length <- 9  
y_length <- 7
print(p_main)
ggsave(
  filename = paste0(FigureFolder, "/NormativeDevCurve_Flanker_blue.pdf"),
  plot = p_main,
  width = x_length,  
  height = y_length,  
  units = "cm"
)

###bootsrap p
age_range <- range(X)
# remove empty data derivative_Flanker 中的缺失值
derivative_Flanker <- derivative_Flanker %>% drop_na(P50_lower, P50_upper)
derivative_Flanker <- derivative_Flanker %>%
  mutate(significance = ifelse(P50_lower > 0, "Increasing",
                               ifelse(P50_upper < 0, "Decreasing", "Non-significant")))


derivative_Flanker$color_group <- ifelse(derivative_Flanker$significance == "Non-significant", NA, derivative_Flanker$P50_mean)
p_bar <- ggplot(derivative_Flanker) +
  geom_bar(aes(x=Age, y=0.5, fill=color_group, colour=color_group), stat="identity", position="stack") +
  scale_fill_gradient2(high="#c24e44", low="white", midpoint=0, na.value="white", labels=NULL) +
  scale_color_gradient2(high="#c24e44", low="white", midpoint=0, na.value="white", labels=NULL) +
  # 使用 annotate() 绘制矩形边框
  annotate("rect", xmin=min(derivative_Flanker$Age), xmax=max(derivative_Flanker$Age), ymin=0, ymax=0.5, 
           color="black", fill=NA, size=0.5) +
  scale_y_continuous(breaks=NULL) +
  ylab(NULL) + xlab("Age (years)") +
  scale_x_continuous(breaks=NULL) +
  labs(fill="mu Mean") +
  theme_void() +
  theme(
    axis.text=element_text(size=8.5, color='black'),
    legend.position="none",
    legend.title=element_text(size=8.5),
    legend.text=element_text(size=8.5),
    plot.margin=unit(c(0, 0, 0, 0), "cm")
  )

x_length <- 9 
y_length_main <- 7 
y_length_bar <- 0.5 
final_plot <- p_main / p_bar + plot_layout(heights = c(y_length_main, y_length_bar))
print(final_plot)
ggsave(
  filename = paste0(FigureFolder, "/NormativeDevCurve_Flanker_Addbar.pdf"), 
  plot = final_plot, 
  width = x_length,  
  height = y_length_main + y_length_bar,  
  units = "cm"
)

lower_CI <- derivative_Flanker$sigma_pred_lower
upper_CI <- derivative_Flanker$sigma_pred_upper

p_sigma <- ggplot(derivative_Flanker, aes(x = Age, y = sigma_pred)) +
  #geom_hline(yintercept = seq(0.055, 0.075, by = 0.005), linetype = "solid", color = "gray90", size = 0.5) + 
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.3, fill = "gray80") +
  geom_line(size = 0.8, color = "black") +
  labs(x = "Age (years)", y = "Standard Deviation", title = "Flanker") +
  theme_minimal() +
  scale_y_continuous(name="SD", limits = c(0.05, 0.08), breaks = seq(0.05, 0.08, by = 0.01)) +
  scale_x_continuous(name="", limits = c(8.8, 16), breaks = seq(9, 16, by = 1)) +
  theme(
    plot.title = element_text(size = 8.5, hjust = 0.5),
    axis.text = element_text(colour = "black",size = 8.5),
    axis.title = element_text(size = 8.5, hjust = 0.5),
    axis.line = element_line(colour = "black", size = 0.5),
    axis.ticks = element_line(colour = "black", size = 0.5), 
    axis.ticks.length = unit(0.2, "cm"),
    axis.text.y = element_text(colour = "black",size = 8.5),   
    panel.grid = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm"))
x_length <- 9  
y_length <- 7
p_sigma
ggsave(
  filename = paste0(FigureFolder, "/NormativesigmaCurve_Flanker.pdf"),
  plot = p_sigma,
  width = x_length, 
  height = y_length,  
  units = "cm"
)
###bootsrap p
age_range <- range(X)
# remove empty data from derivative_Flanker 
derivative_Flanker <- derivative_Flanker %>% drop_na(sigma_lower, sigma_upper)
derivative_Flanker <- derivative_Flanker %>%
  mutate(significance = ifelse(sigma_lower > 0, "Increasing",
                               ifelse(sigma_upper < 0, "Decreasing", "Non-significant")))


derivative_Flanker$color_group <- ifelse(derivative_Flanker$significance == "Non-significant", NA, derivative_Flanker$sigma_mean)
p_sigmabar <- ggplot(derivative_Flanker) +
  geom_bar(aes(x=Age, y=0.5, fill=color_group, colour=color_group), stat="identity", position="stack") +
  scale_fill_gradient2(high="white", low="#c24e44", midpoint=0, na.value="white", labels=NULL) +
  scale_color_gradient2(high="white", low="#c24e44", midpoint=0, na.value="white", labels=NULL) +
  # 使用 annotate() 绘制矩形边框
  annotate("rect", xmin=min(derivative_Flanker$Age), xmax=max(derivative_Flanker$Age), ymin=0, ymax=0.5, 
           color="black", fill=NA, size=0.5) +
  scale_y_continuous(breaks=NULL) +
  ylab(NULL) + xlab("Age (years)") +
  scale_x_continuous(breaks=NULL) +
  labs(fill="P50 mean") +
  theme_void() +
  theme(
    axis.text=element_text(size=8.5, color='black'),
    legend.position="none",
    legend.title=element_text(size=8.5),
    legend.text=element_text(size=8.5),
    plot.margin=unit(c(0, 0, 0, 0), "cm")
  )
x_length <- 9 
y_length_main <- 7 
y_length_bar <- 0.5
sigma_plot <- p_sigma / p_sigmabar + plot_layout(heights = c(y_length_main, y_length_bar))
print(sigma_plot)
ggsave(
  filename = paste0(FigureFolder, "/NormativesigmaCurve_back1Acc_Addbar.pdf"), 
  plot = sigma_plot, 
  width = x_length,  
  height = y_length_main + y_length_bar,  
  units = "cm"
)



