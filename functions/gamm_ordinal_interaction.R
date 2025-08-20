# Filename: gamm_ordinal_interaction.R (Final Version)

library(mgcv)
library(tidyverse)
library(gamm4)

gamm.ordinal.interaction <- function(dependentvar, dataname, smooth_var, int_var, covariates, knots, draws = 1000) {
  
  cat("\n--- Running Ordinal Interaction Analysis for:", dependentvar, " & ", int_var, "---\n")
  
  # --- 1. Data Preparation ---
  gam.data <- get(dataname)
  
  if (!"ID" %in% names(gam.data)) {
    stop("Error: Dataframe must contain a column named 'ID' for random effects.")
  }
  
  # Filter NAs for the specific variables in this model run
  original_rows <- nrow(gam.data)
  gam.data <- gam.data %>%
    filter(!is.na(.data[[dependentvar]]) & !is.na(.data[[int_var]]))
  filtered_rows <- nrow(gam.data)
  cat(sprintf("Data filtered from %d to %d rows after removing NAs.\n", original_rows, filtered_rows))
  
  # --- 【关键修复】: 丢弃未使用的因子水平 ---
  gam.data <- droplevels(gam.data)
  
  # --- 增强的诊断检查 ---
  unique_ids <- length(unique(gam.data$ID))
  cat(sprintf("Diagnostics: Found %d unique IDs after filtering and droplevels().\n", unique_ids))
  
  # 检查协变量的水平
  if (!is.null(covariates) && covariates != "" && covariates %in% names(gam.data)) {
    cov_levels <- length(levels(as.factor(gam.data[[covariates]])))
    cat(sprintf("Diagnostics: Covariate '%s' has %d levels in the current data.\n", covariates, cov_levels))
    if (cov_levels <= 1) {
      cat("!!! WARNING: Covariate has only one level after filtering. This might cause model issues. Proceeding with caution.\n")
    }
  }
  
  if (nrow(gam.data) < 20 || unique_ids <= 1) {
    cat(sprintf("!!! WARNING: Insufficient data. Found %d rows and %d unique IDs. Skipping this analysis. !!!\n", 
                nrow(gam.data), unique_ids))
    return(NULL)
  }
  
  # Prepare the ordinal dependent variable
  gam.data[[dependentvar]] <- as.integer(gam.data[[dependentvar]])
  unique_responses <- sort(unique(gam.data[[dependentvar]]))
  R <- length(unique_responses)
  
  if (min(unique_responses, na.rm = TRUE) < 1) {
    gam.data[[dependentvar]] <- gam.data[[dependentvar]] - min(unique_responses, na.rm = TRUE) + 1
  }
  
  # --- 2. Model Fitting ---
  gam.family <- ocat(R = R)
  
  formula_full <- as.formula(sprintf(
    "%s ~ s(%s, k=%d) + s(%s, by = %s, k=%d) + %s + %s",
    dependentvar, smooth_var, knots, smooth_var, int_var, knots, int_var, covariates
  ))
  
  formula_null <- as.formula(sprintf(
    "%s ~ s(%s, k=%d) + %s + %s",
    dependentvar, smooth_var, knots, int_var, covariates
  ))
  
  model_full <- tryCatch({
    gamm4(formula_full, random = ~(1|ID), data = gam.data, family = gam.family)
  }, error = function(e) {
    cat("!!! ERROR during full model fitting: ", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (is.null(model_full)) return(NULL)
  
  model_null <- tryCatch({
    gamm4(formula_null, random = ~(1|ID), data = gam.data, family = gam.family)
  }, error = function(e) {
    cat("!!! ERROR during null model fitting: ", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (is.null(model_null)) return(NULL)
  
  # --- 3. Statistical Inference & 4. Effect Curve Calculation ---
  # ... (函数的其余部分保持不变) ...
  
  model_comparison <- anova(model_null$gam, model_full$gam, test = "Chisq")
  interaction_p_value <- model_comparison$`Pr(>Chi)`[2]
  
  deviance_full <- summary(model_full$gam)$dev.expl
  deviance_null <- summary(model_null$gam)$dev.expl
  interaction_pseudo_rsq <- deviance_full - deviance_null
  
  stats_results <- data.frame(
    dependentvar = dependentvar,
    int_var = int_var,
    p_value = interaction_p_value,
    pseudo_rsq = interaction_pseudo_rsq
  )
  
  smooth_min <- min(gam.data[[smooth_var]])
  smooth_max <- max(gam.data[[smooth_var]])
  
  pred_grid <- expand.grid(
    smooth_val = seq(smooth_min, smooth_max, length.out = 100),
    int_val = quantile(gam.data[[int_var]], probs = c(0.1, 0.9), na.rm = TRUE)
  )
  names(pred_grid) <- c(smooth_var, int_var)
  
  if (!is.null(covariates) && covariates != "") {
    pred_grid[[covariates]] <- levels(as.factor(gam.data[[covariates]]))[1]
  }
  
  coefs <- coef(model_full$gam)
  vcov_matrix <- vcov(model_full$gam, unconditional = TRUE)
  
  simulated_coefs <- suppressWarnings(MASS::mvrnorm(draws, mu = coefs, Sigma = vcov_matrix))
  
  lp_matrix <- predict(model_full$gam, newdata = pred_grid, type = "lpmatrix")
  
  posterior_preds <- lp_matrix %*% t(simulated_coefs)
  
  n_points <- nrow(pred_grid)
  expected_values <- matrix(NA, nrow = n_points, ncol = draws)
  
  for (d in 1:draws) {
    pred_draw <- matrix(posterior_preds[, d], nrow = n_points, ncol = R - 1, byrow = TRUE)
    cum_probs <- plogis(pred_draw)
    cum_probs <- cbind(cum_probs, 1)
    cat_probs <- t(apply(cum_probs, 1, function(x) diff(c(0, x))))
    expected_values[, d] <- cat_probs %*% 1:R
  }
  
  pred_grid_with_exp <- cbind(pred_grid, expected_values)
  
  effect_over_age <- pred_grid_with_exp %>%
    pivot_longer(cols = starts_with("V"), names_to = "draw", values_to = "expected_phq") %>%
    mutate(ef_level = ifelse(!!sym(int_var) == min(!!sym(int_var)), "low_ef", "high_ef")) %>%
    select(!!sym(smooth_var), draw, ef_level, expected_phq) %>%
    pivot_wider(names_from = "ef_level", values_from = "expected_phq") %>%
    mutate(effect = high_ef - low_ef) %>%
    group_by(!!sym(smooth_var)) %>%
    summarise(
      median_effect = median(effect, na.rm = TRUE),
      lower_ci = quantile(effect, probs = 0.025, na.rm = TRUE),
      upper_ci = quantile(effect, probs = 0.975, na.rm = TRUE)
    )
  
  return(list(stats = stats_results, effect_data = effect_over_age))
}