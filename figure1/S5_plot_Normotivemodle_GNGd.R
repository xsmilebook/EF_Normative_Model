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
library(moments)
library(showtext)
library(patchwork)
# input directory
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- '/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/data0715'
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/interfileFolder_Acc"
  functionFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/Rcode_EFnorms/functions"
  resultFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/results_Acc"
}else{
  datapath <- '/Users/tanlirou/Documents/YF_EF_psy/EF_psy/interfileFolder/Gonogo'
  FigureFolder <- '/Users/tanlirou/Documents/YF_EF_psy/EF_psy/figureFolder/Gonogo'
  interfileFolder <- "/Users/tanlirou/Documents/YF_EF_psy/EF_psy/interfileFolder/Gonogo"
  functionFolder <- "/Users/tanlirou/Documents/YF_EF_psy/EF_psy/code/functions"
  resultFolder <- "/Users/tanlirou/Documents/YF_EF_psy/EF_psy/figureFolder/Gonogo"
}
sumGNGd_prime_deviation <- readRDS(paste0(interfileFolder,"/GNGd_prime.deviations.rds"))
modGNGd_prime.set1.sum <- readRDS(paste0(interfileFolder, "/GAMLSS.GNGd_primeset1.sum.rds"))
modGNGd_prime.set2.sum <- readRDS(paste0(interfileFolder, "/GAMLSS.GNGd_primeset2.sum.rds"))
GNGd_prime_data1 <- read_csv(paste0(interfileFolder, "/GNGd_prime.data1.csv"))
GNGd_prime_data2 <- read_csv(paste0(interfileFolder, "/GNGd_prime.data2.csv"))
derivative_GNGd <- read_csv(paste0(interfileFolder, "/derivative_summary_GNGd.csv"))

fillcolor <- c("#F5B7BF", "#91ACE0")
# 绘制柱状图
gng_plot <- ggplot(data = sumGNGd_prime_deviation, aes(x = Age_year, y = after_stat(count), fill = Sex)) +
  geom_histogram(binwidth = 1, color = "black", boundary = 11, position = "stack", linewidth = 0.5) + 
  labs(x = "Age (year)", y = NULL, title = paste0("Go/No-go, N=", nrow(sumGNGd_prime_deviation))) +
  scale_fill_manual(values = fillcolor, name = "Sex") +
  scale_x_continuous(breaks = seq(11, max(sumGNGd_prime_deviation$Age_year, na.rm = TRUE), by = 1)) + 
  theme_classic() +
  theme(
    aspect.ratio = 0.8,
    plot.title = element_text(color = "black", size = 8.5, hjust = 0.5),
    axis.title = element_text(color = "black", size = 8.5),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    axis.text = element_text(color = "black", size = 8.5)
  )
print(gng_plot)
# 保存图像为 PDF 文件
ggsave(filename = paste0(FigureFolder, "/gngnumber.pdf"), plot = gng_plot, width = 9, height=7, units="cm")



modGNGd_prime.set1 <- modGNGd_prime.set1.sum$performance.tb
modGNGd_prime.set2 <- modGNGd_prime.set2.sum$performance.tb

# 打印出性能表格中的 R² 或相关指标
print(modGNGd_prime.set1)
print(modGNGd_prime.set2)

# 假设性能表格中有 partialRsq（部分决定系数），可以直接提取
partial_Rsq_set1 <- modGNGd_prime.set1$partialRsq
partial_Rsq_set2 <- modGNGd_prime.set2$partialRsq

print(paste("R² for Set1: ", partial_Rsq_set1))
print(paste("R² for Set2: ", partial_Rsq_set2))

# 计算偏度（Skewness）
skewness_value <- skewness(sumGNGd_prime_deviation$d_prime_deviationZ)

# 打印偏度值
print(paste("Skewness: ", skewness_value))

# 计算峰度（Kurtosis）
kurtosis_value <- kurtosis(sumGNGd_prime_deviation$d_prime_deviationZ)

# 打印峰度值
print(paste("Kurtosis: ", kurtosis_value))
# 设置随机种子，确保可重复性
set.seed(123)
# 从大数据集中随机抽取 500 个样本
sample_size <- 5000
sample_data <- sample(sumGNGd_prime_deviation$d_prime_deviationZ, size = sample_size, replace = FALSE)
# 进行 Shapiro-Wilk 检验
shapiro_test_result <- shapiro.test(sample_data)
print(shapiro_test_result)
#####
# ###plot centile
# hist(sumGNGd_prime_deviation$d_prime, breaks = 30, main = "Histogram of d_prime", xlab = "d_prime")
# hist(sumGNGd_prime_deviation$d_prime_centile, main = "Centile Distribution", xlab = "d_prime_centile", ylab = "Frequency", col = "blue")
# 
# 
# ggplot(sumGNGd_prime_deviation, aes(x = d_prime_centile)) +
#   geom_histogram(binwidth = 0.1, fill = "blue", color = "black") +
#   theme_minimal()



## 1. deviation proportion
sumGNGd_prime_deviation[["d_prime_deviationZ_extreme_p"]] <- sumGNGd_prime_deviation[["d_prime_deviationZ"]] > 1.96
sumGNGd_prime_deviation[["d_prime_deviationZ_extreme_n"]] <- sumGNGd_prime_deviation[["d_prime_deviationZ"]] < -1.96
proportion_GNGd_prime_extreme_p <- mean(sumGNGd_prime_deviation[["d_prime_deviationZ_extreme_p"]], na.rm = TRUE)
print(paste("Proportion of positive extreme deviations is", round(proportion_GNGd_prime_extreme_p * 100, 2), "%"))

proportion_GNGd_prime_extreme_n <- mean(sumGNGd_prime_deviation[["d_prime_deviationZ_extreme_n"]], na.rm = TRUE)
print(paste("Proportion of negative extreme deviations is", round(proportion_GNGd_prime_extreme_n * 100, 2), "%"))

#saveRDS(sumGNGd_prime_deviation, paste0(interfileFolder, "/processed_sum_GNGd_prime_deviation.rds"))

# # Plotting extreme points
# pGNGd_primeextreme <- ggplot(sumGNGd_prime_deviation, aes(x = Age_year, y = d_prime)) +
#   geom_point(aes(color = d_prime_deviationZ_extreme_p | d_prime_deviationZ_extreme_n), alpha = 0.6) +
#   scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "blue"), name = "Extreme Deviation") +
#   labs(title = "Deviation Z with Extreme Points", x = "Age", y = "Deviation Z") +
#   theme_minimal()
# # Show plot
# print(pGNGd_primeextreme)
# 
# # plotting extreme line
# ggplot(sumGNGd_prime_deviation, aes(x = d_prime)) +
#   geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
#   geom_vline(xintercept = 1.96, linetype = "dashed", color = "red", size = 1) +
#   geom_vline(xintercept = -1.96, linetype = "dashed", color = "red", size = 1) +
#   labs(title = "Distribution of GNGd_prime Deviation Z-Scores",
#        x = "GNGd_prime Deviation Z-Scores",
#        y = "Count") +
#   theme_minimal()


###2 plot normitive trajectory
# Calculate centiles
n_points <- 1000
centiles.tmp.set1 <- modGNGd_prime.set1.sum$centiles_strat
centiles.tmp.set2 <- modGNGd_prime.set2.sum$centiles_strat

centiles.F.tmp <- (centiles.tmp.set1[[2]] + centiles.tmp.set2[[2]]) / 2
centiles.M.tmp <- (centiles.tmp.set1[[1]] + centiles.tmp.set2[[1]]) / 2
Centiles <- (centiles.F.tmp + centiles.M.tmp) / 2

X  <- seq(min(sumGNGd_prime_deviation$Age_year), max(sumGNGd_prime_deviation$Age_year), length.out=1000)
sumGNGd_prime_deviation$Sex <- factor(sumGNGd_prime_deviation$Sex, levels = c("F", "M"))
# Plotting
#6299ca
p_main <- ggplot() +
  #geom_hline(yintercept = seq(-1, 5, by = 1), linetype = "solid", color = "gray70", size = 0.5) + 
  #geom_point(data=sumGNGd_prime_deviation, aes(x=Age_year, y=d_prime),color = "#fac858",size = 3, alpha=0.25, show.legend=FALSE, stroke = 0) +
  geom_jitter(data=sumGNGd_prime_deviation, aes(x=Age_year, y=d_prime), 
              color = "#A4C5DF", size = 2, alpha = 0.25, show.legend=FALSE, stroke = 0, width = 0, height = 0.1) +
  # stat_density2d(data = sumGNGd_prime_deviation,
  #                aes(x = Age_year, y = d_prime, fill = ..level..),
  #                geom = "polygon",  alpha = 0.1, contour = TRUE, bins = 15, adjust = 0.5) +
  geom_line(aes(x=X, y=Centiles[3,]), linetype="dashed", color = "gray10", size = 0.8) +
  geom_line(aes(x=X, y=Centiles[4,]), linetype="dashed", color = "gray10", size = 0.8) +
  geom_line(aes(x=X, y=Centiles[5,]), linetype="solid", color = "black", size = 0.8) +
  geom_line(aes(x=X, y=Centiles[6,]), linetype="dashed", color = "gray10", size = 0.8) +
  geom_line(aes(x=X, y=Centiles[7,]), linetype="dashed", color = "gray10", size = 0.8) +
  scale_y_continuous(name="d'", limits = c(-1, 5.5),breaks = seq(-1, 5, by = 1))+
  scale_x_continuous(name="", limits = c(11, 18),breaks = seq(11, 18, by = 1))+
  labs(x="", title="Go/No-go") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 8.5, hjust = 0.5, color = "black"),
    axis.text = element_text(colour = "black",size=8.5), 
    axis.title = element_text(size = 8.5,face = "plain"),
    axis.line = element_line(colour = "black", size = 0.5),
    axis.ticks = element_line(colour = "black", size = 0.5), 
    axis.ticks.length = unit(0.2, "cm"),
    axis.text.y = element_text(colour = "black",size = 8.5),  
    axis.line.y = element_line(colour = "black", size = 0.5),
    plot.background=element_rect(fill="white",color = NA),
    legend.position = "none" ,
    panel.grid = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm"))
x_length <- 9  # x 轴的目标长度
y_length <- 7
print(p_main)
ggsave(
  filename = paste0(FigureFolder, "/NormativeDevCurve_GNGd_prime_blue.pdf"),
  plot = p_main,
  width = x_length,  # x轴长度 + 左侧文字的留白
  height = y_length,  # 根据宽高比调整
  units = "cm"
)
#ggsave(paste0(FigureFolder, "/NormativeDevCurve_GNGd_prime_blue.pdf"), width = 20, height=15, units="cm")
###bootsrap p
age_range <- range(X)
# 移除 derivative_GNGd 中的缺失值
derivative_GNGd <- derivative_GNGd %>% drop_na(P50_lower, P50_upper)
derivative_GNGd <- derivative_GNGd %>%
  mutate(significance = ifelse(P50_lower > 0, "Increasing",
                               ifelse(P50_upper < 0, "Decreasing", "Non-significant")))

# # 添加显著性列：判断导数置信区间的方向
derivative_GNGd$color_group <- ifelse(derivative_GNGd$significance == "Non-significant", NA, derivative_GNGd$P50_mean)

p_bar <- ggplot(derivative_GNGd) +
  # 绘制条形图的内容，不添加边框
  geom_bar(aes(x=Age, y=0.5, fill=color_group, colour=color_group), stat="identity", position="stack") +
  scale_fill_gradient2(high="#c24e44", low="white", midpoint=0, na.value="white", labels=NULL) +
  scale_color_gradient2(high="#c24e44", low="white", midpoint=0, na.value="white", labels=NULL) +
  # 使用 annotate() 绘制矩形边框
  annotate("rect", xmin=min(derivative_GNGd$Age), xmax=max(derivative_GNGd$Age), ymin=0, ymax=0.5, 
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
# 合并主图和条形图
#final_plot <- p_main / p_bar + plot_layout(heights=c(10, 1))
# 显示最终图像
#print(final_plot)
#ggsave(paste0(FigureFolder, "/NormativeDevCurve_GNGd_prime_Addbar_blue.pdf"), width = 20, height=15, units="cm")
x_length <- 9 
y_length_main <- 7
y_length_bar <- 0.5 
final_plot <- p_main / p_bar + plot_layout(heights = c(y_length_main, y_length_bar))
print(final_plot)
ggsave(
  filename = paste0(FigureFolder, "/NormativeDevCurve_GNGd_prime_Addbar_blue.pdf"), 
  plot = final_plot, 
  width = x_length,  
  height = y_length_main + y_length_bar,  
  units = "cm"
)

# ggplot() +
#   # 原始数据的散点图，按性别颜色区分
#   #geom_point(data = sumGNGd_prime_deviation, aes(x = Age_year, y = d_prime, color = Sex),size = 3, alpha = 0.5, show.legend = FALSE, stroke = 0) +
#   #M
#   geom_line(aes(x = X, y = centiles.M.tmp[2,]), linetype = "dashed", color = "#3288bd", size = 1) +
#   geom_line(aes(x = X, y = centiles.M.tmp[5,]), linetype = "solid", color = "#3288bd", size = 1.5) +
#   geom_line(aes(x = X, y = centiles.M.tmp[8,]), linetype = "dashed", color = "#3288bd", size = 1) +
#   #F
#   geom_line(aes(x = X, y = centiles.F.tmp[2,]), linetype = "dashed", color = "#d53e4f", size = 1) +
#   geom_line(aes(x = X, y = centiles.F.tmp[5,]), linetype = "solid", color = "#d53e4f", size = 1.5) +
#   geom_line(aes(x = X, y = centiles.F.tmp[8,]), linetype = "dashed", color = "#d53e4f", size = 1) +
#   #scale_color_manual(values = c("M" = "#7fa2ff", "F" = "#f4a3b5")) +  
#   scale_y_continuous(name = "") +  
#   labs(x = "", title = "Go No go trajectory by Sex") +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(size = 22, hjust = 0.5, color = "black"),
#     axis.text = element_text(size = 18), 
#     axis.title = element_text(size = 20, face = "plain"),
#     axis.line = element_line(colour = "black", size = 1),
#     plot.background = element_rect(fill = "white", color = NA),
#     panel.grid = element_blank()
#   )
# ggsave(paste0(FigureFolder, "/NormativeDevCurve_GNGd_prime_bysex.pdf"), width = 20, height=15, units="cm")
# 
# 
# #####sumGNGd_prime_deviation$predicted_sigma <- predict(sumGNGd_prime_deviation, newdata = sumGNGd_prime_deviation, type = "response")
# ####sigma
# ggplot(sumGNGd_prime_deviation, aes(x = Age_year, y = d_prime_sigma)) +
#   geom_smooth(method = "gam", formula = y ~ s(x), size = 1.5, color = "steelblue") +
#   labs(x = "Age(years)", y = "Standard Error", title = "Standard Error trajectory of Go No go d_prime") +
#   theme_minimal() +
#   #scale_y_continuous(breaks = seq(0.74,0.92, by = 0.02)) + 
#   scale_y_continuous(labels = function(y) y) +
#   theme(
#     plot.title = element_text(size = 26, hjust = 0.5),
#     axis.text = element_text(size = 20),
#     axis.title = element_text(size = 22),
#     axis.line = element_line(colour = "black", size = 1),
#     panel.grid = element_blank(),
#     legend.position = "none"
#   )
# ggsave(paste0(FigureFolder, "/NormativesigmaCurve_GNGd_prime.pdf"), width = 20, height=15, units="cm")

# predict_sigma <- data.frame(Age_year = X, Sex = "F")
# # 获取 sigma 参数的预测值
# sigma_pred1 <- predict(modGNGd_prime.set1.sum$mod.tmp, what = "sigma", newdata = predict_sigma, type = "response",se.fit = T,data = GNGd_prime_data1)
# sigma_pred2 <- predict(modGNGd_prime.set2.sum$mod.tmp, what = "sigma", newdata = predict_sigma, type = "response",se.fit = T,data = GNGd_prime_data2)
# # 计算均值和置信区间
# Sigma_mean <- (sigma_pred1 + sigma_pred2) / 2
# # 获取标准误差
# se1 <- sigma_pred1
# se2 <- sigma_pred2
# n1 <- length(se1)
# n2 <- length(se2)
# # 计算合并标准误差
# Sigma_se <- sqrt(((n1 - 1) * se1^2 + (n2 - 1) * se2^2) / (n1 + n2 - 2)) 
# lower_CI <- Sigma_mean - 1.96 * Sigma_se
# upper_CI <- Sigma_mean + 1.96 * Sigma_se
# # 整理数据
# sigma_data <- data.frame(
#   Age_year = X,
#   Sigma = Sigma_mean,
#   lower_CI = lower_CI,
#   upper_CI = upper_CI
# )
lower_CI <- derivative_GNGd$sigma_pred_lower
upper_CI <- derivative_GNGd$sigma_pred_upper
#5470c6
# 绘制图形
p_sigma <- ggplot(derivative_GNGd, aes(x = Age, y = sigma_pred)) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.3, fill = "gray80") +
  #geom_hline(yintercept = seq(0.8, 1, by = 0.05), linetype = "solid", color = "gray90", size = 0.5) + 
  geom_line(size = 0.8, color = "black") +
  labs(x = "Age (years)", y = "Standard Deviation", title = "Go/No-go") +
  theme_minimal() +
  scale_y_continuous(name="SD", limits = c(0.7, 1.02), breaks = seq(0.7, 1, by = 0.1)) +
  scale_x_continuous(name="", limits = c(11, 18), breaks = seq(11, 18, by = 1)) +
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
#ggsave(paste0(FigureFolder, "/NormativesigmaCurve_GNGd_prime.pdf"), width = 20, height=15, units="cm")
ggsave(
  filename = paste0(FigureFolder, "/NormativesigmaCurve_GNGd_prime.pdf"),
  plot = p_sigma,
  width = x_length, 
  height = y_length,  
  units = "cm"
)
###bootsrap p
age_range <- range(X)
# 移除 derivative_GNGd 中的缺失值
derivative_GNGd <- derivative_GNGd %>% drop_na(sigma_lower, sigma_upper)
derivative_GNGd <- derivative_GNGd %>%
  mutate(significance = ifelse(sigma_lower > 0, "Increasing",
                               ifelse(sigma_upper < 0, "Decreasing", "Non-significant")))

# # 添加显著性列：判断导数置信区间的方向
derivative_GNGd$color_group <- ifelse(derivative_GNGd$significance == "Non-significant", NA, derivative_GNGd$sigma_mean)

p_sigmabar <- ggplot(derivative_GNGd) +
  # 绘制条形图的内容，不添加边框
  #313695
  #d7301f
  geom_bar(aes(x=Age, y=0.5, fill=color_group, colour=color_group), stat="identity", position="stack") +
  scale_fill_gradient2(high="white", low="#c24e44", midpoint=0, na.value="white", labels=NULL) +
  scale_color_gradient2(high="white", low="#c24e44", midpoint=0, na.value="white", labels=NULL) +
  # 使用 annotate() 绘制矩形边框
  annotate("rect", xmin=min(derivative_GNGd$Age), xmax=max(derivative_GNGd$Age), ymin=0, ymax=0.5, 
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
# 合并主图和条形图
#sigma_plot <- p_sigma / p_sigmabar + plot_layout(heights=c(10, 1))
# 显示最终图像
#print(sigma_plot)
#ggsave(paste0(FigureFolder, "/NormativesigmaCurve_GNGd_prime_Addbar.pdf"), width = 20, height=15, units="cm")
x_length <- 9 
y_length_main <- 7 
y_length_bar <- 0.5
sigma_plot <- p_sigma / p_sigmabar + plot_layout(heights = c(y_length_main, y_length_bar))
print(sigma_plot)
ggsave(
  filename = paste0(FigureFolder, "/NormativesigmaCurve_GNGd_prime_Addbar.pdf"), 
  plot = sigma_plot, 
  width = x_length,  
  height = y_length_main + y_length_bar,  
  units = "cm"
)



## 筛选0-0.1百分位数的数据
# GNG_bottom10_percent <- sumGNGd_prime_deviation %>%
#   filter(d_prime_centile >= 0 & d_prime_centile <= 0.1) %>%
#   write.xlsx(paste0(interfileFolder,"/GNG_bottom10%.xlsx"))
## 筛选0.45-0.55百分位数的数据
# GNG_middle10_percent <- sumGNGd_prime_deviation %>%
#   filter(d_prime_centile >= 0.45 & d_prime_centile <= 0.55) %>%
#   write.xlsx(paste0(interfileFolder,"/GNG_middle10%.xlsx"))
## 筛选0.9-1百分位数的数据
# GNG_top10_percent <- sumGNGd_prime_deviation %>%
#   filter(d_prime_centile >= 0.9 & d_prime_centile <= 1) %>%
#   write.xlsx(paste0(interfileFolder,"/GNG_top10%.xlsx"))
# GNG_bottom10 <- read.xlsx((paste0(interfileFolder,"/GNG_bottom10%.xlsx")))
# GNG_middle10 <- read.xlsx((paste0(interfileFolder,"/GNG_middle10%.xlsx")))
# GNG_top10 <- read.xlsx((paste0(interfileFolder,"/GNG_top10%.xlsx")))
#####不同人群的gam模型
####top 10%
# source("~/Documents/EF_yunfu_check/code/functions/gamsmooth.R")
# GNG_top10$Sex <- as.factor(GNG_top10$Sex)
# dependentvar <- "d_prime"
# dataname <- "GNG_top10"
# smooth_var <- "Age_year"
# covariates <- "Sex"
# knots <- 3
# set_fx = FALSE
# stats_only = FALSE
# mod_only = FALSE
# Gam.GNG_top10.mod <- gam.fit.smooth(dependentvar, dataname, smooth_var, covariates, knots, set_fx = FALSE, stats_only = FALSE, mod_only = T)
# # 使用 predict 获取平滑曲线的预测值
# GNG_top10$Sex <- as.factor(GNG_top10$Sex)  # 确保 Sex 是因子类型
# GNG_top10_pred <- data.frame(
#   Age_year = seq(min(GNG_top10$Age_year), max(GNG_top10$Age_year), length.out = 1000),
#   Sex = levels(GNG_top10$Sex)[1]  # 可以选择 "M" 或 "F"，这里以 "M" 为例
# )
# GNG_top_predictions <- predict(Gam.GNG_top10.mod, newdata = GNG_top10_pred, type = "response", se.fit = TRUE)
# GNG_top10_pred$d_prime <- GNG_top_predictions$fit
# GNG_top10_pred$ymin <- GNG_top_predictions$fit - 1.96 * GNG_top_predictions$se.fit
# GNG_top10_pred$ymax <- GNG_top_predictions$fit + 1.96 * GNG_top_predictions$se.fit
# write_xlsx(GNG_top10_pred,"/Users/tanlirou/Documents/EF_yunfu_check/results/GNG_top_pred.xlsx")
# 
# ####middle 10%
# GNG_middle10$Sex <- as.factor(GNG_middle10$Sex)
# dependentvar <- "d_prime"
# dataname <- "GNG_middle10"
# smooth_var <- "Age_year"
# covariates <- "Sex"
# knots <- 3
# set_fx = FALSE
# stats_only = FALSE
# mod_only = FALSE
# Gam.GNG_middle10.mod <- gam.fit.smooth(dependentvar, dataname, smooth_var, covariates, knots, set_fx = FALSE, stats_only = FALSE, mod_only = T)
# # 使用 predict 获取平滑曲线的预测值
# GNG_middle10$Sex <- as.factor(GNG_middle10$Sex)  # 确保 Sex 是因子类型
# GNG_middle10_pred <- data.frame(
#   Age_year = seq(min(GNG_middle10$Age_year), max(GNG_middle10$Age_year), length.out = 1000),
#   Sex = levels(GNG_middle10$Sex)[1]  # 可以选择 "M" 或 "F"，这里以 "M" 为例
# )
# GNG_middle_predictions <- predict(Gam.GNG_middle10.mod, newdata = GNG_middle10_pred, type = "response", se.fit = TRUE)
# GNG_middle10_pred$d_prime <- GNG_middle_predictions$fit
# GNG_middle10_pred$ymin <- GNG_middle_predictions$fit - 1.96 * GNG_middle_predictions$se.fit
# GNG_middle10_pred$ymax <- GNG_middle_predictions$fit + 1.96 * GNG_middle_predictions$se.fit
# write_xlsx(GNG_middle10_pred,"/Users/tanlirou/Documents/EF_yunfu_check/results/GNG_middle_pred.xlsx")
# 
# ####bottom 10%
# GNG_bottom10$Sex <- as.factor(GNG_bottom10$Sex)
# dependentvar <- "d_prime"
# dataname <- "GNG_bottom10"
# smooth_var <- "Age_year"
# covariates <- "Sex"
# knots <- 3
# set_fx = FALSE
# stats_only = FALSE
# mod_only = FALSE
# Gam.GNG_bottom10.mod <- gam.fit.smooth(dependentvar, dataname, smooth_var, covariates, knots, set_fx = FALSE, stats_only = FALSE, mod_only = T)
# # 使用 predict 获取平滑曲线的预测值
# GNG_bottom10$Sex <- as.factor(GNG_bottom10$Sex)  # 确保 Sex 是因子类型
# GNG_bottom10_pred <- data.frame(
#   Age_year = seq(min(GNG_bottom10$Age_year), max(GNG_bottom10$Age_year), length.out = 1000),
#   Sex = levels(GNG_bottom10$Sex)[1]  # 可以选择 "M" 或 "F"，这里以 "M" 为例
# )
# # 获取预测值
# GNG_bottom_predictions <- predict(Gam.GNG_bottom10.mod, newdata = GNG_bottom10_pred, type = "response", se.fit = TRUE)
# GNG_bottom10_pred$d_prime <- GNG_bottom_predictions$fit
# GNG_bottom10_pred$ymin <- GNG_bottom_predictions$fit - 1.96 * GNG_bottom_predictions$se.fit
# GNG_bottom10_pred$ymax <- GNG_bottom_predictions$fit + 1.96 * GNG_bottom_predictions$se.fit
# write_xlsx(GNG_bottom10_pred,"/Users/tanlirou/Documents/EF_yunfu_check/results/GNG_bottom_pred.xlsx")


# #load data
# GNG_top10_pred <- read_xlsx("/Users/tanlirou/Documents/EF_yunfu_check/results/GNG_top_pred.xlsx")
# GNG_middle10_pred <- read_xlsx("/Users/tanlirou/Documents/EF_yunfu_check/results/GNG_middle_pred.xlsx")
# GNG_bottom10_pred <- read_xlsx("/Users/tanlirou/Documents/EF_yunfu_check/results/GNG_bottom_pred.xlsx")

# pGNG3percentile <- ggplot() +
#   geom_point(data = GNG_bottom10, aes(x = Age_year, y = d_prime), color = "gray70", size = 1, alpha = 0.5, stroke = 0) +
#   geom_line(data = GNG_bottom10_pred, aes(x = Age_year, y = d_prime), color = "#7b92c7", fill = "#d2d6f5", size = 1) +
#   geom_point(data = GNG_middle10, aes(x = Age_year, y = d_prime), color = "gray70", size = 1, alpha = 0.5, stroke = 0) +
#   geom_line(data = GNG_middle10_pred, aes(x = Age_year, y = d_prime), color = "#37939a", fill = "#d3f0f2", size = 1) +
#   geom_point(data = GNG_top10, aes(x = Age_year, y = d_prime), color = "gray70", size = 1, alpha = 0.5, stroke = 0) +
#   geom_line(data = GNG_top10_pred, aes(x = Age_year, y = d_prime), color = "#304e7e", fill = "#aed4e5", size = 1) +
#   theme(panel.background = element_blank(),
#         axis.line = element_line(colour = "black", size = 1)) +
#   theme(axis.text.y = element_text(size = 16)) +
#   theme(axis.text.x = element_text(size = 16)) +
#   theme(axis.title = element_text(size = 14)) +
#   theme(plot.title = element_text(size = 16, hjust = 0.5)) +
#   scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.1),
#                      labels = function(y) y) +
#   scale_x_continuous(breaks = seq(from = 12, to = 18, by = 2),
#                      labels = function(x) x) +
#   labs(title = "Go No go", x = "Age", y = "Mean accuracy") +
#   scale_color_manual(values = c("#7b92c7", "#37939a", "#304e7e")) +
#   scale_fill_manual(values = c("#d2d6f5", "#d3f0f2", "#aed4e5"))
# print(pGNG3percentile)
# ggsave(paste0(FigureFolder, "/GNG3percentile_plot_point.pdf"), plot = pGNG3percentile, width = 20, height = 15, units = "cm")
# 
# 
# GNG_main_plot <- ggplot() + 
#   # 分别绘制不同组的密度图
#   stat_density2d(data = GNG_bottom10, aes(x = Age_year, y = d_prime, fill = "Bottom 10%"), 
#                  geom = "polygon", alpha = 0.2, contour = TRUE, bins = 8, adjust = 1) +
#   stat_density2d(data = GNG_middle10, aes(x = Age_year, y = d_prime, fill = "Middle 10%"), 
#                  geom = "polygon", alpha = 0.2, contour = TRUE, bins = 8, adjust = 1) +
#   stat_density2d(data = GNG_top10, aes(x = Age_year, y = d_prime, fill = "Top 10%"), 
#                  geom = "polygon", alpha = 0.2, contour = TRUE, bins = 8, adjust = 1) +
#   geom_ribbon(data = GNG_bottom10_pred, aes(x = Age_year, ymin = ymin, ymax = ymax), fill = "#807dba", alpha = 0.5) +
#   geom_line(data = GNG_bottom10_pred, aes(x = Age_year, y = d_prime), color = "#807dba",  size = 1) +
#   geom_ribbon(data = GNG_middle10_pred, aes(x = Age_year, ymin = ymin, ymax = ymax), fill = "#045a8d", alpha = 0.5) +
#   geom_line(data = GNG_middle10_pred, aes(x = Age_year, y = d_prime), color = "#045a8d",  size = 1) +
#   geom_ribbon(data = GNG_top10_pred, aes(x = Age_year, ymin = ymin, ymax = ymax), fill = "#016c59", alpha = 0.5) +
#   geom_line(data = GNG_top10_pred, aes(x = Age_year, y = d_prime), color = "#016c59",  size = 1) +
#   theme(panel.background = element_blank(),
#         axis.line = element_line(colour = "black", size = 1),
#         axis.text.y = element_text(size = 18),
#         axis.text.x = element_text(size = 18),
#         axis.title = element_text(size = 20),
#         plot.title = element_text(size = 22, hjust = 0.5),
#         legend.position = "none") +  # 去掉图例
#   # 设置坐标轴
#   scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2), limits = c(0.2, 1), labels = function(y) y) +
#   scale_x_continuous(breaks = seq(from = 12, to = 18, by = 2), labels = function(x) x) +
#   # 设置标题和标签
#   labs(title = "Percentile Trajectories in Go No go", x = "Age(years)", y = "Accuracy") +
#   # 使用 scale_fill_manual 设置不同的颜色
#   scale_fill_manual(values = c("Bottom 10%" = "#9e9ac8", "Middle 10%" = "#0570b0", "Top 10%" = "#02818a"), name = "Group")
# print(GNG_main_plot)
# ggsave(paste0(FigureFolder, "/GNG3percentile_plot.pdf"), plot = GNG_main_plot, width = 20, height = 15, units = "cm")
# 
# GNG_legend <- ggplot() + 
#   # 分别绘制不同组的密度图
#   stat_density2d(data = GNG_bottom10, aes(x = Age_year, y = d_prime, fill = "Bottom 10%"), 
#                  geom = "polygon", alpha = 1, contour = TRUE, bins = 8, adjust = 1) +
#   stat_density2d(data = GNG_middle10, aes(x = Age_year, y = d_prime, fill = "Middle 10%"), 
#                  geom = "polygon", alpha = 1, contour = TRUE, bins = 8, adjust = 1) +
#   stat_density2d(data = GNG_top10, aes(x = Age_year, y = d_prime, fill = "Top 10%"), 
#                  geom = "polygon", alpha = 1, contour = TRUE, bins = 8, adjust = 1) +
#   geom_ribbon(data = GNG_bottom10_pred, aes(x = Age_year, ymin = ymin, ymax = ymax), fill = "#807dba", alpha = 0.5) +
#   geom_line(data = GNG_bottom10_pred, aes(x = Age_year, y = d_prime), color = "#807dba",  size = 1) +
#   geom_ribbon(data = GNG_middle10_pred, aes(x = Age_year, ymin = ymin, ymax = ymax), fill = "#045a8d", alpha = 0.5) +
#   geom_line(data = GNG_middle10_pred, aes(x = Age_year, y = d_prime), color = "#045a8d",  size = 1) +
#   geom_ribbon(data = GNG_top10_pred, aes(x = Age_year, ymin = ymin, ymax = ymax), fill = "#016c59", alpha = 0.5) +
#   geom_line(data = GNG_top10_pred, aes(x = Age_year, y = d_prime), color = "#016c59",  size = 1) +
#   theme(panel.background = element_blank(),
#         axis.line = element_line(colour = "black", size = 1),
#         axis.text.y = element_text(size = 18),
#         axis.text.x = element_text(size = 18),
#         axis.title = element_text(size = 20),
#         plot.title = element_text(size = 22, hjust = 0.5),
#         legend.position = "bottom") +  
#   # 设置坐标轴
#   scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2), limits = c(0.2, 1), labels = function(y) y) +
#   scale_x_continuous(breaks = seq(from = 12, to = 18, by = 2), labels = function(x) x) +
#   # 设置标题和标签
#   labs(title = "Percentile Trajectories in Go No go", x = "Age(years)", y = "Accuracy") +
#   # 使用 scale_fill_manual 设置不同的颜色
#   scale_fill_manual(values = c("Bottom 10%" = "#9e9ac8", "Middle 10%" = "#0570b0", "Top 10%" = "#02818a"), name = "Group")
# print(GNG_legend)
# ggsave(paste0(FigureFolder, "/GNG3percentile_plot_legend.pdf"), plot = GNG_legend, width = 20, height = 15, units = "cm")


