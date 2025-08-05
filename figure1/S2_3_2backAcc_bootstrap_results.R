# import library
library(tidyverse)
library(ggplot2)

# Set file directories
bootstrap_dir <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/datasets/yunfu/interfile_folder/bootstrap/2backAcc"
output_dir <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/datasets/yunfu/interfile_folder/bootstrap"
figure_dir <- "/ibmgpfs/cuizaixu_lab/xuhaoshu/datasets/yunfu/figures/fig1"
file_list <- list.files(path = bootstrap_dir, pattern = "*.rds", full.names = TRUE)

# Initialize lists to store all derivative results
mu_derivatives_all <- list()
sigma_derivatives_all <- list()
sigma_all <- list()

# read derivative results
for (file in file_list) {
  bootstrap_result <- readRDS(file)
  mu_derivatives_all[[length(mu_derivatives_all) + 1]] <- bootstrap_result$P50_derivative
  sigma_derivatives_all[[length(sigma_derivatives_all) + 1]] <- bootstrap_result$sigma_derivative
  sigma_all[[length(sigma_all) + 1]] <- bootstrap_result$sigma_pred
}

# Combine the derivative data into matrices
mu_derivatives_matrix <- do.call(cbind, mu_derivatives_all)
sigma_derivatives_matrix <- do.call(cbind, sigma_derivatives_all)
sigma_matrix <- do.call(cbind, sigma_all)

# compute the median number of mu
P50_mean <- rowMeans(mu_derivatives_matrix)
P50_ci <- apply(mu_derivatives_matrix, 1, quantile, probs = c(0.025, 0.975))

# compute the mean and confidence interval of sigma derivative
sigma_mean <- rowMeans(sigma_derivatives_matrix)
sigma_ci <- apply(sigma_derivatives_matrix, 1, quantile, probs = c(0.025, 0.975))

# compute the mean and confidence interval of sigma prediction
sigma_pred <- rowMeans(sigma_matrix)
sigma_pred_ci <- apply(sigma_matrix, 1, quantile, probs = c(0.025, 0.975))

# Extract age points from the first file
age_points <- readRDS(file_list[1])$Age_points

# Create a summary data frame
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

# Plot the mean and confidence interval of the derivative of mu as a function of age
p1 <- ggplot(derivative_summary, aes(x = Age, y = P50_mean)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = P50_lower, ymax = P50_upper), alpha = 0.2, fill = "blue") +
  labs(title = "Mean and confidence interval of the derivative of μ varying with age", x = "age", y = "derivative of μ") +
  theme_minimal()

# Plot the mean and confidence interval of the derivative of sigma as a function of age
p2 <- ggplot(derivative_summary, aes(x = Age, y = sigma_mean)) +
  geom_line(color = "red") +
  geom_ribbon(aes(ymin = sigma_lower, ymax = sigma_upper), alpha = 0.2, fill = "red") +
  labs(title = "Mean and confidence interval of the derivative of σ varying with age", x = "age", y = "derivative of σ") +
  theme_minimal()

# Save the summary data and the plots
write.csv(derivative_summary, file = file.path(output_dir, "2backAcc", "derivative_summary_2backAcc.csv"), row.names = FALSE)
ggsave(file.path(figure_dir, "mu_derivative_age_2backAcc.png"), plot = p1, width = 8, height = 6, dpi = 300)
ggsave(file.path(figure_dir, "sigma_derivative_age_2backAcc.png"), plot = p2, width = 8, height = 6, dpi = 300)