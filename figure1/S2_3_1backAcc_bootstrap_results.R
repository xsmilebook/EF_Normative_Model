# 加载必要的库
library(tidyverse)
library(ggplot2)

# 设置文件目录
bootstrap_dir <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psy/interfileFolder/bootstrap/1backAcc"
output_dir <- "/ibmgpfs/cuizaixu_lab/tanlirou1/Yunfu/YF_EF_psy/interfileFolder/bootstrap"
file_list <- list.files(path = bootstrap_dir, pattern = "*.rds", full.names = TRUE)

# 初始化列表以保存所有导数结果
mu_derivatives_all <- list()
sigma_derivatives_all <- list()
sigma_all <- list()
# 读取并提取每个文件的导数结果
for (file in file_list) {
  bootstrap_result <- readRDS(file)
  mu_derivatives_all[[length(mu_derivatives_all) + 1]] <- bootstrap_result$P50_derivative
  sigma_derivatives_all[[length(sigma_derivatives_all) + 1]] <- bootstrap_result$sigma_derivative
  sigma_all[[length(sigma_all) + 1]] <- bootstrap_result$sigma_pred
}

# 合并导数数据到矩阵形式，每列代表一次bootstrap的结果
mu_derivatives_matrix <- do.call(cbind, mu_derivatives_all)
sigma_derivatives_matrix <- do.call(cbind, sigma_derivatives_all)
sigma_matrix <- do.call(cbind, sigma_all)

# 计算 mu 的 P50（中位数）及置信区间
P50_mean <- rowMeans(mu_derivatives_matrix)
P50_ci <- apply(mu_derivatives_matrix, 1, quantile, probs = c(0.025, 0.975))

sigma_mean <- rowMeans(sigma_derivatives_matrix)
sigma_ci <- apply(sigma_derivatives_matrix, 1, quantile, probs = c(0.025, 0.975))

sigma_pred <- rowMeans(sigma_matrix)
sigma_pred_ci <- apply(sigma_matrix, 1, quantile, probs = c(0.025, 0.975))

# 从第一个文件中提取年龄点信息
age_points <- readRDS(file_list[1])$Age_points

# 将结果保存为数据框
derivative_summary <- data.frame(
  Age = age_points,
  P50_mean = P50_mean,
  P50_lower = P50_ci[1, ],
  P50_upper = P50_ci[2, ],
  sigma_mean = sigma_mean,
  sigma_lower = sigma_ci[1, ],
  sigma_upper = sigma_ci[2, ],
  sigma_pred = sigma_pred,
  sigma_pred_lower = sigma_pred_ci[1, ],
  sigma_pred_upper = sigma_pred_ci[2, ]
)

# 绘制 mu 的导数随年龄变化的均值和置信区间
p1 <- ggplot(derivative_summary, aes(x = Age, y = P50_mean)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = P5_lower, ymax = P50_upper), alpha = 0.2, fill = "blue") +
  labs(title = "μ的导数随年龄变化的均值和置信区间", x = "年龄", y = "μ的导数") +
  theme_minimal()

# 绘制 sigma 的导数随年龄变化的均值和置信区间
p2 <- ggplot(derivative_summary, aes(x = Age, y = sigma_mean)) +
  geom_line(color = "red") +
  geom_ribbon(aes(ymin = sigma_lower, ymax = sigma_upper), alpha = 0.2, fill = "red") +
  labs(title = "σ的导数随年龄变化的均值和置信区间", x = "年龄", y = "σ的导数") +
  theme_minimal()

# 保存结果数据和图像
write.csv(derivative_summary, file = file.path(output_dir, "derivative_summary_1back.csv"), row.names = FALSE)
ggsave(file.path(output_dir, "mu_derivative_age_1back.png"), plot = p1, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir, "sigma_derivative_age_1back.png"), plot = p2, width = 8, height = 6, dpi = 300)