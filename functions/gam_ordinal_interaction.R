# Filename: gam_ordinal_interaction.R
# Description: 
#   Final version. Returns effect curves in BOTH the expected value space (for accuracy)
#   and the latent linear predictor space (for smooth visualization).

# --- Required Libraries ---
library(mgcv)
library(tidyverse)
library(ecostats)

# --- Function Definition ---

gam.ordinal.interaction <- function(
    dependentvar, 
    dataname, 
    smooth_var, 
    int_var, 
    covariates, 
    knots, 
    draws = 1000, 
    bootstrap_sims = 1000
) {
  
  # --- Steps 1, 2, 3 (No major changes) ---
  cat("\n--- Running Ordinal Varying-Coefficient GAM for:", dependentvar, "&", int_var, "---\n")
  gam.data <- get(dataname)
  gam.data <- gam.data %>% filter(!is.na(.data[[dependentvar]]) & !is.na(.data[[int_var]]))
  if (nrow(gam.data) < 50) { cat("!!! WARNING: Insufficient data. Skipping. !!!\n"); return(NULL) }
  gam.data[[dependentvar]] <- as.integer(gam.data[[dependentvar]])
  unique_responses <- sort(unique(gam.data[[dependentvar]]))
  if (min(unique_responses, na.rm = TRUE) == 0) {
    gam.data[[dependentvar]] <- gam.data[[dependentvar]] + 1
    cat("Adjusted dependent variable to start from 1.\n")
    unique_responses <- unique_responses + 1
  }
  R <- length(unique_responses)
  gam.family <- ocat(R = R)
  formula_full <- as.formula(sprintf("%s ~ s(%s, k=%d) + s(%s, by = %s, k=%d) + %s + %s", dependentvar, smooth_var, knots, smooth_var, int_var, knots, int_var, covariates))
  formula_null <- as.formula(sprintf("%s ~ s(%s, k=%d) + %s + %s", dependentvar, smooth_var, knots, int_var, covariates))
  cat("Fitting full & null models...\n")
  model_full <- gam(formula_full, data = gam.data, family = gam.family, method = "REML")
  model_null <- gam(formula_null, data = gam.data, family = gam.family, method = "REML")
  model_comparison <- anova(model_null, model_full, test = "Chisq")
  anova_p_value <- model_comparison$`Pr(>Chi)`[2]
  cat("Attempting bootstrap for p-value...\n")
  bootstrap_p_value <- NA
  tryCatch({
    bootstrap_p_value <- anovaPB(model_null, model_full, n.sim = bootstrap_sims, test = 'Chisq')$`Pr(>Chi)`[2]
  }, warning = function(w) { cat("--- INFO: anovaPB not implemented. Using standard ANOVA p-value.\n") }, error = function(e) { cat("--- INFO: anovaPB not implemented. Using standard ANOVA p-value.\n") })
  deviance_full <- summary(model_full)$dev.expl
  deviance_null <- summary(model_null)$dev.expl
  interaction_pseudo_rsq <- deviance_full - deviance_null
  stats_results <- data.frame(dependentvar = dependentvar, int_var = int_var, anova_p = anova_p_value, bootstrap_p = bootstrap_p_value, p_value_final = ifelse(!is.na(bootstrap_p_value), bootstrap_p_value, anova_p_value), pseudo_rsq = interaction_pseudo_rsq, n_obs = nrow(gam.data))
  
  # --- 4. Calculate Effect Curves (for both latent and expected value space) ---
  cat("Calculating effect curves...\n")
  
  # Create a common prediction grid
  smooth_min <- min(gam.data[[smooth_var]])
  smooth_max <- max(gam.data[[smooth_var]])
  pred_grid <- expand.grid(smooth_val = seq(smooth_min, smooth_max, length.out = 100), int_val = quantile(gam.data[[int_var]], probs = c(0.1, 0.9), na.rm = TRUE))
  names(pred_grid) <- c(smooth_var, int_var)
  if (!is.null(covariates) && covariates != "") {
    gam.data[[covariates]] <- as.factor(gam.data[[covariates]])
    mode_level <- names(which.max(table(gam.data[[covariates]])))
    pred_grid[[covariates]] <- factor(mode_level, levels = levels(gam.data[[covariates]]))
  }
  
  # Posterior simulation
  coefs <- coef(model_full)
  vcov_matrix <- vcov(model_full, unconditional = TRUE)
  simulated_coefs <- suppressWarnings(MASS::mvrnorm(draws, mu = coefs, Sigma = vcov_matrix))
  
  # This gets predictions in the linear predictor (latent) space for the *first category*
  lp_matrix <- predict(model_full, newdata = pred_grid, type = "lpmatrix")
  posterior_preds <- lp_matrix %*% t(simulated_coefs)
  
  # --- A. Calculate SMOOTH effect curve in the LATENT SPACE ---
  
  latent_df <- cbind(pred_grid, as.data.frame(posterior_preds))
  
  effect_data_latent <- latent_df %>%
    pivot_longer(cols = starts_with("V"), names_to = "draw", values_to = "latent_value") %>%
    mutate(ef_level = ifelse(!!sym(int_var) == min(!!sym(int_var)), "low_ef", "high_ef")) %>%
    select(!!sym(smooth_var), draw, ef_level, latent_value) %>%
    pivot_wider(names_from = "ef_level", values_from = "latent_value") %>%
    mutate(effect = high_ef - low_ef) %>% # The effect is the difference on the linear predictor scale
    group_by(!!sym(smooth_var)) %>%
    summarise(
      median_effect = median(effect, na.rm = TRUE),
      lower_ci = quantile(effect, probs = 0.025, na.rm = TRUE),
      upper_ci = quantile(effect, probs = 0.975, na.rm = TRUE)
    )
  
  # --- B. Calculate NON-SMOOTH effect curve in the EXPECTED VALUE SPACE ---
  
  n_points <- nrow(pred_grid)
  expected_values <- matrix(NA, nrow = n_points, ncol = draws)
  for (d in 1:draws) {
    pred_draw <- matrix(posterior_preds[, d], nrow = n_points, ncol = R - 1, byrow = TRUE)
    cum_probs <- plogis(pred_draw)
    cum_probs <- cbind(cum_probs, 1)
    cat_probs <- t(apply(cum_probs, 1, function(x) diff(c(0, x))))
    expected_values[, d] <- cat_probs %*% unique_responses
  }
  
  expected_values_df <- as.data.frame(expected_values)
  names(expected_values_df) <- paste0("draw_", 1:draws)
  pred_grid_with_exp <- cbind(pred_grid, expected_values_df)
  
  effect_data_expected <- pred_grid_with_exp %>%
    pivot_longer(cols = starts_with("draw_"), names_to = "draw", values_to = "expected_phq") %>%
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
  
  # --- 5. Return Results ---
  cat("--- Analysis Complete ---\n")
  return(list(
    stats = stats_results, 
    effect_data_expected = effect_data_expected, # The "bumpy" curve
    effect_data_latent = effect_data_latent     # The SMOOTH curve
  ))
}