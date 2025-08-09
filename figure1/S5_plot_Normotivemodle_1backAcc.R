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
# input directory
wd <- getwd()
if (str_detect(wd, "cuizaixu_lab")){
  datapath <- '/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/data0715'
  demopath <- '/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/data0715'
  interfileFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/interfileFolder_Acc"
  functionFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/Rcode_EFnorms/functions"
  resultFolder <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/results_Acc"
}else{
  datapath <- 'D:/datasets/yunfu/raw_data'
  FigureFolder <- 'D:/datasets/yunfu/figures/fig1'
  interfileFolder <- "D:/datasets/yunfu/interfile_folder/Normative_Model"
  functionFolder <- "D:/code/EF_Normative_Model/functions"
  resultFolder <- "D:/datasets/yunfu/results/gamlss/Normative_Model"
}

taskname <- "1-back"
sum1backAcc_deviation <- readRDS(file.path(interfileFolder, taskname, "back1Acc.deviations.rds"))
mod1backAcc.set1.sum <- readRDS(file.path(interfileFolder, taskname, "GAMLSS.back1Accset1.sum.rds"))
mod1backAcc.set2.sum <- readRDS(file.path(interfileFolder, taskname, "GAMLSS.back1Accset2.sum.rds"))
back1Acc_data1 <- read_csv(file.path(interfileFolder, taskname, "back1Acc.data1.csv"))
back1Acc_data2 <- read_csv(file.path(interfileFolder, taskname, "back1Acc.data2.csv"))
derivative_back1 <- read_csv(file.path(interfileFolder, taskname, "derivative_summary_1back.csv"))

fillcolor <- c("#F5B7BF", "#91ACE0")
back1_plot <- ggplot(data = sum1backAcc_deviation, aes(x = Age_year, y = after_stat(count), fill = Sex)) +
  geom_histogram(binwidth = 1, color = "black", boundary = 11, position = "stack", linewidth = 0.5) + 
  labs(x = "Age (year)", y = NULL, title = paste0("1-back, N=", nrow(sum1backAcc_deviation))) +
  scale_fill_manual(values = fillcolor, name = "Sex") +
  scale_x_continuous(breaks = seq(11, max(sum1backAcc_deviation$Age_year, na.rm = TRUE), by = 1)) + 
  theme_classic() +
  theme(
    aspect.ratio = 0.8,
    plot.title = element_text(color = "black", size = 8.5, hjust = 0.5),
    axis.title = element_text(color = "black", size = 8.5),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    axis.text = element_text(color = "black", size = 8.5),
    legend.text = element_text(size = 8.5),
    legend.title = element_text(size = 8.5)
  )
print(back1_plot)
# 保存图像为 PDF 文件
ggsave(filename = paste0(FigureFolder, "/onebacknumber.pdf"), plot = back1_plot, width = 9, height=7, units="cm")


mod1backAcc.set1 <- mod1backAcc.set1.sum$performance.tb
mod1backAcc.set2 <- mod1backAcc.set2.sum$performance.tb

# 打印出性能表格中的 R² 或相关指标
print(mod1backAcc.set1)
print(mod1backAcc.set2)

# 假设性能表格中有 partialRsq（部分决定系数），可以直接提取
partial_Rsq_set1 <- mod1backAcc.set1$partialRsq
partial_Rsq_set2 <- mod1backAcc.set2$partialRsq

print(paste("R² for Set1: ", partial_Rsq_set1))
print(paste("R² for Set2: ", partial_Rsq_set2))

# 计算偏度（Skewness）
skewness_value <- skewness(sum1backAcc_deviation$Oneback_acc_deviationZ)

# 打印偏度值
print(paste("Skewness: ", skewness_value))

# 计算峰度（Kurtosis）
kurtosis_value <- kurtosis(sum1backAcc_deviation$Oneback_acc_deviationZ)

# 打印峰度值
print(paste("Kurtosis: ", kurtosis_value))
# 设置随机种子，确保可重复性
set.seed(123)
 # 从大数据集中随机抽取 500 个样本
sample_size <- 5000
sample_data <- sample(sum1backAcc_deviation$Oneback_acc_deviationZ, size = sample_size, replace = FALSE)
# 进行 Shapiro-Wilk 检验
shapiro_test_result <- shapiro.test(sample_data)
print(shapiro_test_result)


# ###plot centiles
# ###plot
# hist(sum1backAcc_deviation$Oneback_acc, breaks = 30, main = "Histogram of 1back_acc", xlab = "1back_acc")
# hist(sum1backAcc_deviation$Oneback_acc_deviationZ, breaks = 30, main = "Histogram of 1back_accZ", xlab = "1back_acc")
# hist(sum1backAcc_deviation$Oneback_acc_centile, 
#      main = "Centile Distribution", xlab = "Oneback_acc_centile", ylab = "Frequency", col = "blue")


## 1. deviation proportion
sum1backAcc_deviation[["Oneback_acc_deviationZ_extreme_p"]] <- sum1backAcc_deviation[["Oneback_acc_deviationZ"]] > 1.96
sum1backAcc_deviation[["Oneback_acc_deviationZ_extreme_n"]] <- sum1backAcc_deviation[["Oneback_acc_deviationZ"]] < -1.96
proportion_1backAcc_extreme_p <- mean(sum1backAcc_deviation[["Oneback_acc_deviationZ_extreme_p"]], na.rm = TRUE)
print(paste("Proportion of positive extreme deviations is", round(proportion_1backAcc_extreme_p * 100, 2), "%"))

proportion_1backAcc_extreme_n <- mean(sum1backAcc_deviation[["Oneback_acc_deviationZ_extreme_n"]], na.rm = TRUE)
print(paste("Proportion of negative extreme deviations is", round(proportion_1backAcc_extreme_n * 100, 2), "%"))

saveRDS(sum1backAcc_deviation, paste0(interfileFolder, "/processed_sum_1backAcc_deviation.rds"))

# # Plotting extreme points
# p1backAccextreme <- ggplot(sum1backAcc_deviation, aes(x = Age_year, y = Oneback_acc)) +
#   geom_point(aes(color = Oneback_acc_deviationZ_extreme_p | Oneback_acc_deviationZ_extreme_n), alpha = 0.6) +
#   scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "blue"), name = "Extreme Deviation") +
#   labs(title = "Deviation Z with Extreme Points", x = "Age", y = "Deviation Z") +
#   theme_minimal()
# # Show plot
# print(p1backAccextreme)
# # plotting extreme line
# ggplot(sum1backAcc_deviation, aes(x = Oneback_acc)) +
#   geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
#   geom_vline(xintercept = 1.96, linetype = "dashed", color = "red", size = 1) +
#   geom_vline(xintercept = -1.96, linetype = "dashed", color = "red", size = 1) +
#   labs(title = "Distribution of 1backAcc Deviation Z-Scores",
#        x = "1backAcc Deviation Z-Scores",
#        y = "Count") +
#   theme_minimal()


###2 plot normitive trajectory
# Calculate centiles
n_points <- 1000
centiles.tmp.set1 <- mod1backAcc.set1.sum$centiles_strat
centiles.tmp.set2 <- mod1backAcc.set2.sum$centiles_strat

centiles.F.tmp <- (centiles.tmp.set1[[2]] + centiles.tmp.set2[[2]]) / 2
centiles.M.tmp <- (centiles.tmp.set1[[1]] + centiles.tmp.set2[[1]]) / 2
Centiles <- (centiles.F.tmp + centiles.M.tmp) / 2

X  <- seq(min(sum1backAcc_deviation$Age_year), max(sum1backAcc_deviation$Age_year), length.out=1000)

# Plotting
p_main  <- ggplot() +
  #geom_hline(yintercept = seq(0.2, 1, by = 0.2), linetype = "solid", color = "gray70", size = 0.5) + 
  geom_point(data=sum1backAcc_deviation, aes(x=Age_year, y=Oneback_acc), 
             color= "#A4C5DF", size = 2, alpha=0.25,show.legend=FALSE, stroke = 0) +
  # stat_density2d(data = sumGNGd_prime_deviation,
  #                aes(x = Age_year, y = d_prime, fill = ..level..),
  #                geom = "polygon",  alpha = 0.1, contour = TRUE, bins = 15, adjust = 0.5) +
  geom_line(aes(x=X, y=Centiles[3,]), linetype="dashed", color = "gray10", size = 0.8) +
  geom_line(aes(x=X, y=Centiles[4,]), linetype="dashed", color = "gray10", size = 0.8) +
  geom_line(aes(x=X, y=Centiles[5,]), linetype="solid", color = "black", size = 0.8) +
  geom_line(aes(x=X, y=Centiles[6,]), linetype="dashed", color = "gray10", size = 0.8) +
  geom_line(aes(x=X, y=Centiles[7,]), linetype="dashed", color = "gray10", size = 0.8) +
  scale_y_continuous(name="Acc",limits = c(0.1, 1),breaks = seq(0.2, 1, by = 0.2))+
  scale_x_continuous(name="", limits = c(11, 18),breaks = seq(11, 18, by = 1))+
  labs(x="", title="1-back") +
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
x_length <- 9  # x 轴的目标长度
y_length <- 7
print(p_main)
ggsave(
  filename = paste0(FigureFolder, "/NormativeDevCurve_1backAcc_blue.pdf"),
  plot = p_main,
  width = x_length,  # x轴长度 + 左侧文字的留白
  height = y_length,  # 根据宽高比调整
  units = "cm"
)
# print(p_main)
# ggsave(paste0(FigureFolder, "/NormativeDevCurve_1backAcc_blue.pdf"), width = 20, height=15, units="cm")

###bootsrap p
age_range <- range(X)
# 移除 derivative_back1 中的缺失值
derivative_back1 <- derivative_back1 %>% drop_na(P50_lower, P50_upper)
derivative_back1 <- derivative_back1 %>%
  mutate(significance = ifelse(P50_lower > 0, "Increasing",
                               ifelse(P50_upper < 0, "Decreasing", "Non-significant")))

# # 添加显著性列：判断导数置信区间的方向
derivative_back1$color_group <- ifelse(derivative_back1$significance == "Non-significant", NA, derivative_back1$P50_mean)

p_bar <- ggplot(derivative_back1) +
  # 绘制条形图的内容，不添加边框
  geom_bar(aes(x=Age, y=0.5, fill=color_group, colour=color_group), stat="identity", position="stack") +
  scale_fill_gradient2(high="#c24e44", low="white", midpoint=0, na.value="white", labels=NULL) +
  scale_color_gradient2(high="#c24e44", low="white", midpoint=0, na.value="white", labels=NULL) +
  # 使用 annotate() 绘制矩形边框
  annotate("rect", xmin=min(derivative_back1$Age), xmax=max(derivative_back1$Age), ymin=0, ymax=0.5, 
           color="black", fill=NA, size=0.5) +
  scale_y_continuous(breaks=NULL) +
  ylab(NULL) + xlab("Age (years)") +
  scale_x_continuous(breaks=NULL) +
  labs(fill="mu Mean") +
  theme_void() +
  theme(
    axis.text=element_text(size=18.5, color='black'),
    legend.position="none",
    legend.title=element_text(size=12),
    legend.text=element_text(size=10),
    plot.margin=unit(c(0, 0, 0, 0), "cm")
  )
# # 合并主图和条形图
# final_plot <- p_main / p_bar + plot_layout(heights=c(10, 1))
# # 显示最终图像
# print(final_plot)
# ggsave(paste0(FigureFolder, "/NormativeDevCurve_1backAcc_Addbar_blue.pdf"), width = 20, height=15, units="cm")
x_length <- 9
y_length_main <- 7 
y_length_bar <- 0.5 
final_plot <- p_main / p_bar + plot_layout(heights = c(y_length_main, y_length_bar))
print(final_plot)
ggsave(
  filename = paste0(FigureFolder, "/NormativeDevCurve_1backAcc_Addbar_blue.pdf"), 
  plot = final_plot, 
  width = x_length,  
  height = y_length_main + y_length_bar,  
  units = "cm"
)


# ggplot() +
#   #M
#   geom_line(aes(x = X, y = centiles.M.tmp[2,]), linetype = "dashed", color = "#3288bd", size = 1) +
#   geom_line(aes(x = X, y = centiles.M.tmp[5,]), linetype = "solid", color = "#3288bd", size = 1.5) +
#   geom_line(aes(x = X, y = centiles.M.tmp[8,]), linetype = "dashed", color = "#3288bd", size = 1) +
#   #F
#   geom_line(aes(x = X, y = centiles.F.tmp[2,]), linetype = "dashed", color = "#d53e4f", size = 1) +
#   geom_line(aes(x = X, y = centiles.F.tmp[5,]), linetype = "solid", color = "#d53e4f", size = 1.5) +
#   geom_line(aes(x = X, y = centiles.F.tmp[8,]), linetype = "dashed", color = "#d53e4f", size = 1) +
#   #scale_color_manual(values = c("M" = "#7fa2ff", "F" = "#f4a3b5")) +  
#   scale_y_continuous((name = "Accuracy"), breaks = seq(from = 0, to = 1, by = 0.2) )+  
#   labs(x = "", title = "1back trajectory by Sex") +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(size = 22, hjust = 0.5, color = "black"),
#     axis.text = element_text(size = 18), 
#     axis.title = element_text(size = 20, face = "plain"),
#     axis.line = element_line(colour = "black", size = 1),
#     plot.background = element_rect(fill = "white", color = NA),
#     panel.grid = element_blank()
#   )
# ggsave(paste0(FigureFolder, "/NormativeDevCurve_1backAcc_bysex.pdf"), width = 20, height=15, units="cm")

####sigma
# ggplot(sum1backAcc_deviation, aes(x = Age_year, y = Oneback_acc_sigma) )+
#   geom_smooth(method = "gam", formula = y ~ s(x), size = 1.5,color = "steelblue") +
#   labs(x = "", y = "Standard Error", title = "Standard Error trajectory of 1back") +
#   theme_minimal() +
#   scale_y_continuous(breaks = seq(from = 0, to = 0.28, by = 0.02),  labels = function(y) y) +
#   theme(
#     plot.title = element_text(size = 28, hjust = 0.5),
#     axis.text = element_text(size = 20),
#     axis.title = element_text(size = 22),
#     axis.line = element_line(colour = "black", size = 1),
#     panel.grid = element_blank(),
#     legend.position = "none"
#   )
# ggsave(paste0(FigureFolder, "/NormativesigmaCurve_1backAcc.pdf"), width = 20, height=15, units="cm")



# predict_sigma <- data.frame(Age_year = X, Sex = "F")
# # 获取 sigma 参数的预测值
# sigma_pred1 <- predict(mod1backAcc.set1.sum$mod.tmp, what = "sigma", newdata = predict_sigma, type = "response",se.fit = TRUE,data = back1Acc_data1)
# sigma_pred2 <- predict(mod1backAcc.set2.sum$mod.tmp, what = "sigma", newdata = predict_sigma, type = "response",se.fit = TRUE,data = back1Acc_data2)
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
lower_CI <- derivative_back1$sigma_pred_lower
upper_CI <- derivative_back1$sigma_pred_upper
# 绘制图形
p_sigma <- ggplot(derivative_back1, aes(x = Age, y = sigma_pred)) +
  #geom_hline(yintercept = seq(0.15, 0.27, by = 0.03), linetype = "solid", color = "gray90", size = 0.5) + 
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.3, fill = "gray80") +
  geom_line(size = 0.8, color = "black") +
  labs(x = "Age (years)", y = "Standard Deviation", title = "1-back") +
  theme_minimal() +
  scale_y_continuous(name="SD", limits = c(0.1, 0.28), breaks = seq(0.1, 0.28, by = 0.05)) +
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
ggsave(
  filename = paste0(FigureFolder, "/NormativesigmaCurve_back1Acc.pdf"),
  plot = p_sigma,
  width = x_length, 
  height = y_length,  
  units = "cm"
)
#ggsave(paste0(FigureFolder, "/NormativesigmaCurve_back1Acc.pdf"), width = 20, height=15, units="cm")
###bootsrap p
age_range <- range(X)
# 移除 derivative_back1 中的缺失值
derivative_back1 <- derivative_back1 %>% drop_na(sigma_lower, sigma_upper)
derivative_back1 <- derivative_back1 %>%
  mutate(significance = ifelse(sigma_lower > 0, "Increasing",
                               ifelse(sigma_upper < 0, "Decreasing", "Non-significant")))

# # 添加显著性列：判断导数置信区间的方向
derivative_back1$color_group <- ifelse(derivative_back1$significance == "Non-significant", NA, derivative_back1$sigma_mean)

p_sigmabar <- ggplot(derivative_back1) +
  # 绘制条形图的内容，不添加边框
  geom_bar(aes(x=Age, y=0.5, fill=color_group, colour=color_group), stat="identity", position="stack") +
  scale_fill_gradient2(high="white", low="#c24e44", midpoint=0, na.value="white", labels=NULL) +
  scale_color_gradient2(high="white", low="#c24e44", midpoint=0, na.value="white", labels=NULL) +
  # 使用 annotate() 绘制矩形边框
  annotate("rect", xmin=min(derivative_back1$Age), xmax=max(derivative_back1$Age), ymin=0, ymax=0.5, 
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
# # 合并主图和条形图
# sigma_plot <- p_sigma / p_sigmabar + plot_layout(heights=c(10, 1))
# # 显示最终图像
# print(sigma_plot)
# ggsave(paste0(FigureFolder, "/NormativesigmaCurve_back1Acc_Addbar.pdf"), width = 20, height=15, units="cm")
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



predict_mu <- data.frame(Age_year = X, Sex = "F")
# 获取 sigma 参数的预测值
mu_pred1 <- predict(mod1backAcc.set1.sum$mod.tmp, what = "mu", newdata = predict_mu, type = "response",se.fit = TRUE,data = back1Acc_data1)
mu_pred2 <- predict(mod1backAcc.set2.sum$mod.tmp, what = "mu", newdata = predict_mu, type = "response",se.fit = TRUE,data = back1Acc_data2)
# 计算均值和置信区间
mu_mean <- (mu_pred1 + mu_pred2) / 2
# 获取标准误差
se1 <- mu_pred1
se2 <- mu_pred2
n1 <- length(se1)
n2 <- length(se2)
# 计算合并标准误差
mu_se <- sqrt(((n1 - 1) * se1^2 + (n2 - 1) * se2^2) / (n1 + n2 - 2)) 
lower_CI <- mu_mean - 1.96 * mu_se
upper_CI <- mu_mean + 1.96 * mu_se
# 整理数据
mu_data <- data.frame(
  Age_year = X,
  mu = P50_mean,
  lower_CI = lower_CI,
  upper_CI = upper_CI
)

# 绘制图形
ggplot(mu_data, aes(x = Age_year, y = mu)) +
  geom_line(size = 1.5, color = "steelblue") +
  #geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.3, fill = "lightblue") +
  labs(x = "Age (years)", y = "Standard Error", title = "Standard Error Trajectories of 1back Accuracy") +
  theme_minimal() +
  #scale_y_continuous(breaks = seq(0,0.32, by = 0.02)) + 
  theme(
    plot.title = element_text(size = 26, hjust = 0.5),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 22),
    axis.line = element_line(colour = "black", size = 1),
    panel.grid = element_blank()
  )



# ggplot(sum1backAcc_deviation, aes(x = Age_year, y = Oneback_acc_sigma, color = Sex)) +
#   geom_smooth(method = "gam", formula = y ~ s(x), size = 1.5) +
#   labs(x = "Age(years)", y = "Sigma", title = "Sigma of 1 back") +
#   theme_minimal() +
#   scale_color_manual(values = c("M" = "#3288bd", "F" = "#d53e4f")) +  # 指定颜色
#   scale_y_continuous(breaks = seq(from = 0, to = 0.25, by = 0.02),  labels = function(y) y) +
#   theme(
#     plot.title = element_text(size = 22, hjust = 0.5),
#     axis.text = element_text(size = 18),
#     axis.title = element_text(size = 20),
#     axis.line = element_line(colour = "black", size = 1),
#     panel.grid = element_blank(),
#     legend.position = "bottom"
#   )
# ggsave(paste0(FigureFolder, "/Sex_legend.pdf"), width = 20, height=15, units="cm")

#####画不同水平人群的线
# # 筛选0-0.1百分位数的数据
# back1_bottom10_percent <- sum1backAcc_deviation %>%
#   filter(Oneback_acc_centile >= 0 & Oneback_acc_centile <= 0.1) %>%
#   write.xlsx(paste0(interfileFolder,"/Onebackbottom10%.xlsx"))
# # 筛选0.45-0.55百分位数的数据
# back1_middle10_percent <- sum1backAcc_deviation %>%
#   filter(Oneback_acc_centile >= 0.45 & Oneback_acc_centile <= 0.55) %>%
#   write.xlsx(paste0(interfileFolder,"/Onebackmiddle10%.xlsx"))
# # 筛选0.9-1百分位数的数据
# back1_top10_percent <- sum1backAcc_deviation %>%
#   filter(Oneback_acc_centile >= 0.9 & Oneback_acc_centile <= 1) %>%
#   write.xlsx(paste0(interfileFolder,"/Onebacktop10%.xlsx"))

back1_top10 <- read.xlsx((paste0(interfileFolder,"/Onebacktop10%.xlsx")))
back1_bottom10 <- read.xlsx((paste0(interfileFolder,"/Onebackbottom10%.xlsx")))
back1_middle10 <- read.xlsx((paste0(interfileFolder,"/Onebackmiddle10%.xlsx")))

# source("~/Documents/EF_yunfu_check/code/functions/gamsmooth.R")
# back1_top10$Sex <- as.factor(back1_top10$Sex)
# dependentvar <- "Oneback_acc"
# dataname <- "back1_top10"
# smooth_var <- "Age_year"
# covariates <- "Sex"
# knots <- 3
# set_fx = FALSE
# stats_only = FALSE
# mod_only = FALSE
# Gam.back1_top10.mod <- gam.fit.smooth(dependentvar, dataname, smooth_var, covariates, knots, set_fx = FALSE, stats_only = FALSE, mod_only = T)
# # 使用 predict 获取平滑曲线的预测值
# back1_top10$Sex <- as.factor(back1_top10$Sex)  # 确保 Sex 是因子类型
# back1_top10_pred <- data.frame(
#   Age_year = seq(min(back1_top10$Age_year), max(back1_top10$Age_year), length.out = 1000),
#   Sex = levels(back1_top10$Sex)[1]  # 可以选择 "M" 或 "F"，这里以 "M" 为例
# )
# back1_top_predictions <- predict(Gam.back1_top10.mod, newdata = back1_top10_pred, type = "response", se.fit = TRUE)
# back1_top10_pred$Oneback_acc <- back1_top_predictions$fit
# back1_top10_pred$ymin <- back1_top_predictions$fit - 1.96 * back1_top_predictions$se.fit
# back1_top10_pred$ymax <- back1_top_predictions$fit + 1.96 * back1_top_predictions$se.fit
# write_xlsx(back1_top10_pred,"/Users/tanlirou/Documents/EF_yunfu_check/results/back1_top_pred.xlsx")
# 
# ####middle 10%
# back1_middle10$Sex <- as.factor(back1_middle10$Sex)
# dependentvar <- "Oneback_acc"
# dataname <- "back1_middle10"
# smooth_var <- "Age_year"
# covariates <- "Sex"
# knots <- 3
# set_fx = FALSE
# stats_only = FALSE
# mod_only = FALSE
# Gam.back1_middle10.mod <- gam.fit.smooth(dependentvar, dataname, smooth_var, covariates, knots, set_fx = FALSE, stats_only = FALSE, mod_only = T)
# # 使用 predict 获取平滑曲线的预测值
# back1_middle10$Sex <- as.factor(back1_middle10$Sex)  # 确保 Sex 是因子类型
# back1_middle10_pred <- data.frame(
#   Age_year = seq(min(back1_middle10$Age_year), max(back1_middle10$Age_year), length.out = 1000),
#   Sex = levels(back1_middle10$Sex)[1]  # 可以选择 "M" 或 "F"，这里以 "M" 为例
# )
# back1_middle_predictions <- predict(Gam.back1_middle10.mod, newdata = back1_middle10_pred, type = "response", se.fit = TRUE)
# back1_middle10_pred$Oneback_acc <- back1_middle_predictions$fit
# back1_middle10_pred$ymin <- back1_middle_predictions$fit - 1.96 * back1_middle_predictions$se.fit
# back1_middle10_pred$ymax <- back1_middle_predictions$fit + 1.96 * back1_middle_predictions$se.fit
# write_xlsx(back1_middle10_pred,"/Users/tanlirou/Documents/EF_yunfu_check/results/back1_middle_pred.xlsx")
# 
# ####bottom 10%
# back1_bottom10$Sex <- as.factor(back1_bottom10$Sex)
# dependentvar <- "Oneback_acc"
# dataname <- "back1_bottom10"
# smooth_var <- "Age_year"
# covariates <- "Sex"
# knots <- 3
# set_fx = FALSE
# stats_only = FALSE
# mod_only = FALSE
# Gam.back1_bottom10.mod <- gam.fit.smooth(dependentvar, dataname, smooth_var, covariates, knots, set_fx = FALSE, stats_only = FALSE, mod_only = T)
# # 使用 predict 获取平滑曲线的预测值
# back1_bottom10$Sex <- as.factor(back1_bottom10$Sex)  # 确保 Sex 是因子类型
# back1_bottom10_pred <- data.frame(
#   Age_year = seq(min(back1_bottom10$Age_year), max(back1_bottom10$Age_year), length.out = 1000),
#   Sex = levels(back1_bottom10$Sex)[1]  # 可以选择 "M" 或 "F"，这里以 "M" 为例
# )
# # 获取预测值
# back1_bottom_predictions <- predict(Gam.back1_bottom10.mod, newdata = back1_bottom10_pred, type = "response", se.fit = TRUE)
# back1_bottom10_pred$Oneback_acc <- back1_bottom_predictions$fit
# back1_bottom10_pred$ymin <- back1_bottom_predictions$fit - 1.96 * back1_bottom_predictions$se.fit
# back1_bottom10_pred$ymax <- back1_bottom_predictions$fit + 1.96 * back1_bottom_predictions$se.fit
# write_xlsx(back1_bottom10_pred,"/Users/tanlirou/Documents/EF_yunfu_check/results/back1_bottom_pred.xlsx")


##load
back1_top10_pred <- read_xlsx("/Users/tanlirou/Documents/EF_yunfu_check/results/back1_top_pred.xlsx")
back1_middle10_pred <- read_xlsx("/Users/tanlirou/Documents/EF_yunfu_check/results/back1_middle_pred.xlsx")
back1_bottom10_pred <- read_xlsx("/Users/tanlirou/Documents/EF_yunfu_check/results/back1_bottom_pred.xlsx")

back1_main_plot <- ggplot() + 
  # 分别绘制不同组的密度图
  stat_density2d(data = back1_bottom10, aes(x = Age_year, y = Oneback_acc, fill = "Bottom 10%"), 
                 geom = "polygon", alpha = 0.2, contour = TRUE, bins = 8, adjust = 1) +
  stat_density2d(data = back1_middle10, aes(x = Age_year, y = Oneback_acc, fill = "Middle 10%"), 
                 geom = "polygon", alpha = 0.2, contour = TRUE, bins = 8, adjust = 1) +
  stat_density2d(data = back1_top10, aes(x = Age_year, y = Oneback_acc, fill = "Top 10%"), 
                 geom = "polygon", alpha = 0.2, contour = TRUE, bins = 8, adjust = 1) +
  geom_ribbon(data = back1_bottom10_pred, aes(x = Age_year, ymin = ymin, ymax = ymax), fill = "#807dba", alpha = 0.5) +
  geom_line(data = back1_bottom10_pred, aes(x = Age_year, y = Oneback_acc), color = "#807dba",  size = 1) +
  geom_ribbon(data = back1_middle10_pred, aes(x = Age_year, ymin = ymin, ymax = ymax), fill = "#045a8d", alpha = 0.5) +
  geom_line(data = back1_middle10_pred, aes(x = Age_year, y = Oneback_acc), color = "#045a8d",  size = 1) +
  geom_ribbon(data = back1_top10_pred, aes(x = Age_year, ymin = ymin, ymax = ymax), fill = "#016c59", alpha = 0.5) +
  geom_line(data = back1_top10_pred, aes(x = Age_year, y = Oneback_acc), color = "#016c59",  size = 1) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 1),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 22, hjust = 0.5),
        legend.position = "none") +  # 去掉图例
  # 设置坐标轴
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2), limits = c(0.2, 1), labels = function(y) y) +
  scale_x_continuous(breaks = seq(from = 12, to = 18, by = 2), labels = function(x) x) +
  # 设置标题和标签
  labs(title = "Percentile Trajectories in 1back", x = "Age(years)", y = "Accuracy") +
  # 使用 scale_fill_manual 设置不同的颜色
  scale_fill_manual(values = c("Bottom 10%" = "#9e9ac8", "Middle 10%" = "#0570b0", "Top 10%" = "#02818a"), name = "Group")
print(back1_main_plot)
ggsave(paste0(FigureFolder, "/back13percentile_plot.pdf"), plot = back1_main_plot, width = 20, height = 15, units = "cm")

