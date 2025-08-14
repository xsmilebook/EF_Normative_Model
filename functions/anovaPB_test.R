library(mgcv)
library(tidyr)
library(dplyr)

source("/ibmgpfs/cuizaixu_lab/xuhaoshu/EF_Normative_Model/functions/gamcog_withsmoothvar_deviation.R")

# --- 1. 模拟测试数据 ---
set.seed(123)
n <- 400
test_data <- tibble(
  y = rnorm(n),                           # 响应变量
  x1 = rnorm(n),                          # 线性协变量
  x2 = runif(n, 0, 10),                   # 非线性变量（用于样条）
  group = factor(sample(c("A", "B"), n, replace = TRUE))
)

# 生成真实的非线性关系
test_data$y <- with(test_data, 
                    2 + 0.5 * x1 + 
                      sin(x2) +  # 真实非线性
                      as.numeric(group) + 
                      rnorm(n, sd = 0.5))

# --- 2. 定义变量（对应 %s 占位符） ---
response_var   <- "y"      # 第一个 %s
linear_pred1   <- "x1"     # 第二个 %s
linear_pred2   <- "group"  # 第三个 %s  
smooth_pred    <- "x2"     # 第四个 %s（要检验的平滑项）

# --- 3. 构建 Null 和 Alternative 模型公式 ---
# Null 模型：不含平滑项（可选：移除或改为线性）
null_formula_str <- sprintf("%s ~ %s + %s", response_var, linear_pred1, linear_pred2)
alt_formula_str  <- sprintf("%s ~ %s + %s + s(%s, k=3, fx=FALSE)", 
                            response_var, linear_pred1, linear_pred2, smooth_pred)

null_formula <- as.formula(null_formula_str)
alt_formula  <- as.formula(alt_formula_str)

cat("Null model: ", null_formula_str, "\n")
cat("Alt model:  ", alt_formula_str,  "\n\n")

# --- 4. 拟合模型 ---
# 注意：gam 默认用 GCV，但 LRT 检验建议用 ML 或 REML
m0 <- gam(null_formula, data = test_data, method = "REML")
m1 <- gam(alt_formula,  data = test_data, method = "REML")

# --- 5. 使用 anovaPB_ext 进行参数自助检验 ---
# 确保 anovaPB_ext 支持 gam 模型（基于 gam 或 gamObject）

# 如果 anovaPB_ext 不支持 gam，可尝试：
# - 改用 logLik + 手动 LRT
# - 或确保 anovaPB_ext 能处理 mgcv::gam 对象

# 假设 anovaPB_ext 支持 gam：
set.seed(456)
pb_test <- anovaPB_ext(
  objectNull = m0,
  object     = m1,
  n.sim      = 99,  
  test = 'Chisq'  
)

# --- 6. 输出结果 ---
print(pb_test)

cat("\n=== PB-ANOVA Results ===\n")
cat("Observed test stat:  ", attr(pb_test, "observed_stat"), "\n")
cat("Simulated stats (first 5): ", head(attr(pb_test, "simulated_stats"), 5), "\n")
cat("P-value:             ", pb_test$`Pr(>Chi)`[2], "\n")
cat("Simulated stats length:", length(attr(pb_test, "simulated_stats")), "\n")

# --- 7. 可视化模拟分布 ---
sim_stats <- attr(pb_test, "simulated_stats")
obs_stat  <- attr(pb_test, "observed_stat")

hist(sim_stats, breaks = 20, 
     main = "PB Null Distribution vs Observed Statistic",
     xlab = "Likelihood Ratio Statistic", 
     col = "lightblue", border = "white")
abline(v = obs_stat, col = "red", lwd = 2, lty = 2)
legend("topright", 
       legend = sprintf("Observed: %.3f", obs_stat), 
       col = "red", lwd = 2, lty = 2, bty = "n")