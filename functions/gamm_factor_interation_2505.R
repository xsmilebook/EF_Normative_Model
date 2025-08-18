library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse) # 包含 dplyr, stringr 等
library(lme4)
library(gamm4)
library(pbkrtest)
library(parallel)  # 用于并行处理
library(psych)     # for corr.test in gamm.smooth.predict.interaction
# library(tableone) # for CreateTableOne in your main script
# library(writexl)  # for write.xlsx in your main script


# HELPER FUNCTION (if needed, based on your str_split_i)
# If str_split_i is from a package, ensure that package is loaded.
# If it's custom, define it here. For example:
# str_split_i_custom <- function(string, pattern, i) {
#   return(str_split(string, pattern, simplify = TRUE)[,i])
# }
# And then use str_split_i_custom in pbootint.
# OR, more directly, modify pbootint:

pbootint <- function(modelobj, int_var=NA){
  numsims <- 1000 # Can be reduced for speed, e.g., 199 or 249 for permutations
  set.seed(925) 
  df <- modelobj$gam$model
  thisResp <- as.character(modelobj$gam$terms[[2]])
  f1 <- modelobj$gam$formula
  theseVars <- attr(terms(f1),"term.labels")
  
  if (!is.na(int_var)){ # This block seems for a specific interaction structure
    # Ensure stringr is loaded for str_split and str_detect
    if (!requireNamespace("stringr", quietly = TRUE)) {
      stop("Package 'stringr' needed for pbootint. Please install it.", call. = FALSE)
    }
    indexcomma <- gregexpr(",", theseVars[1])[[1]]
    addsmooth <- theseVars[1]
    addsmooth <- paste0(substr(addsmooth, 1, indexcomma[1]), substr(addsmooth, indexcomma[2]+1, nchar(addsmooth)))
    theseVars <-  c(theseVars, addsmooth)
  }
  
  if (sum(stringr::str_detect(theseVars, "by ="))==1){
    # Using stringr::str_split and then selecting the part
    split_by <- stringr::str_split(theseVars[1], "by = ")[[1]]
    if (length(split_by) > 1) {
      int_var_detected <- stringr::str_split(split_by[2], ", ")[[1]][1]
      f2 <- reformulate(c(theseVars[2:(length(theseVars))], int_var_detected), response = thisResp)
    } else {
      # Fallback or error if "by =" is present but structure is not as expected
      warning("pbootint: 'by =' detected but structure for f2 is unclear. Using original logic for f2.")
      f2 <- reformulate(theseVars[2:(length(theseVars))], response = thisResp) # Original fallback
    }
  } else {
    f2 <- reformulate(theseVars[2:(length(theseVars))], response = thisResp)
  }
  
  g1_res <- tryCatch(gam(f1, data = df), error = function(e) NULL)
  g2_res <- tryCatch(gam(f2, data = df), error = function(e) NULL)
  
  if(is.null(g1_res) || is.null(g2_res)) {
    warning("pbootint: gam model fitting failed for g1 or g2.")
    return(NA)
  }
  g1 <- g1_res
  g2 <- g2_res
  
  mat1 <- model.matrix(g1)
  mat2 <- model.matrix(g2)
  
  if (!"src_subject_id" %in% names(df)) {
    warning("pbootint: src_subject_id not found in the model frame. Returning NA.")
    return(NA)
  }
  src_subject_id_vec <- df$src_subject_id
  
  if (length(unique(src_subject_id_vec)) < 2) {
    warning("pbootint: Not enough unique levels in src_subject_id for lmer. Returning NA.")
    return(NA)
  }
  
  y_val <- df[,thisResp] # y was defined earlier. Re-check usage. Using y_val to avoid conflict.
  
  m1_res <- tryCatch({
    lmer(y_val ~ -1 + mat1 + (1|src_subject_id_vec), control = lmerControl(check.nobs.vs.nRE = "ignore"))
  }, error = function(e) {
    warning(paste("pbootint: lmer for m1 failed:", e$message))
    return(NULL)
  })
  if (is.null(m1_res)) return(NA)
  m1 <- m1_res
  
  m2_res <- tryCatch({
    lmer(y_val ~ -1 + mat2 + (1|src_subject_id_vec), control = lmerControl(check.nobs.vs.nRE = "ignore"))
  }, error = function(e) {
    warning(paste("pbootint: lmer for m2 failed:", e$message))
    return(NULL)
  })
  if (is.null(m2_res)) return(NA)
  m2 <- m2_res
  
  pb_res <- tryCatch({
    # PBrefdist can be slow. For permutations, consider reducing nsim if this is a bottleneck.
    refdist <- PBrefdist(m1, m2, nsim = max(20, round(numsims/10))) # Reduced nsim for speed in PBrefdist
    PBmodcomp(m1, m2, ref = refdist)
  }, error = function(e) {
    warning(paste("pbootint: PBmodcomp/PBrefdist failed:", e$message))
    return(NULL)
  })
  
  if (is.null(pb_res) || !("PBtest" %in% rownames(pb_res$test)) ) return(NA)
  
  int_pval <- pb_res$test["PBtest","p.value"]
  return(int_pval)
}


gamm.smooth.predict.interaction <- function(dependentvar, data_input, smoothvar, interest.indep.var, covariates){
  gam.data <- data_input
  
  required_cols_gamm <- c(dependentvar, interest.indep.var, "src_subject_id")
  if (!is.na(smoothvar)) required_cols_gamm <- c(required_cols_gamm, smoothvar)
  if (!is.na(covariates) && all(nchar(covariates) > 0)) required_cols_gamm <- c(required_cols_gamm, covariates)
  
  if (!all(required_cols_gamm %in% names(gam.data))) {
    missing_cols <- required_cols_gamm[!required_cols_gamm %in% names(gam.data)]
    warning(paste("gamm.smooth.predict.interaction: Missing required columns:", paste(missing_cols, collapse=", ")))
    return(data.frame( parcel = dependentvar, interest.indep.var = interest.indep.var, partialRsq = NA, bettermodel=NA,
                       gamm.independent.t = NA, gamm.independent.pvalue = NA, boots.pvalues = NA,
                       correstimate = NA, corrp = NA, slope = NA ))
  }
  
  gam.data <- gam.data[stats::complete.cases(gam.data[, required_cols_gamm]), ]
  
  if (nrow(gam.data) < 20 || length(unique(gam.data$src_subject_id)) < 2) {
    warning(paste0("gamm.smooth.predict.interaction: Not enough data (", nrow(gam.data), 
                   " rows, ", length(unique(gam.data$src_subject_id)), " unique IDs) after NA removal. Returning NA results."))
    return(data.frame( parcel = dependentvar, interest.indep.var = interest.indep.var, partialRsq = NA, bettermodel=NA,
                       gamm.independent.t = NA, gamm.independent.pvalue = NA, boots.pvalues = NA,
                       correstimate = NA, corrp = NA, slope = NA ))
  }
  
  # Outlier removal for dependent variable
  tmp <- gam.data[[dependentvar]] 
  mean_tmp <- mean(tmp, na.rm=TRUE)
  sd_tmp <- sd(tmp, na.rm=TRUE)
  if (!is.na(sd_tmp) && sd_tmp > 0) { # only if sd is calculable and > 0
    outlierindx<-which(tmp < mean_tmp - 3*sd_tmp | tmp > mean_tmp + 3*sd_tmp)
    if (length(outlierindx)>0){
      gam.data<-gam.data[-outlierindx, ]
    }
  }
  if (nrow(gam.data) < 20) { 
    warning("gamm.smooth.predict.interaction: Not enough data after dependent var outlier removal. Returning NA.")
    return(data.frame(parcel = dependentvar, interest.indep.var = interest.indep.var, partialRsq = NA, bettermodel=NA,
                      gamm.independent.t = NA, gamm.independent.pvalue = NA, boots.pvalues = NA,
                      correstimate = NA, corrp = NA, slope = NA))
  }
  
  # Outlier removal for independent variable
  Independent.var <- gam.data[[interest.indep.var]]
  mean_indep <- mean(Independent.var, na.rm = TRUE)
  sd_indep <- sd(Independent.var, na.rm = TRUE)
  if (!is.na(sd_indep) && sd_indep > 0) {
    outlierindx1 <- which(Independent.var < mean_indep - 3 * sd_indep | Independent.var > mean_indep + 3 * sd_indep)
    if (length(outlierindx1)>0){
      gam.data[outlierindx1, interest.indep.var] <- NA
      gam.data <- gam.data[stats::complete.cases(gam.data[, interest.indep.var]), ] # Remove rows with new NAs
    }
  }
  if (nrow(gam.data) < 20 || length(unique(gam.data$src_subject_id)) < 2) {
    warning("gamm.smooth.predict.interaction: Not enough data or src_subject_id levels after indep var outlier removal. Returning NA.")
    return(data.frame(parcel = dependentvar, interest.indep.var = interest.indep.var, partialRsq = NA, bettermodel=NA,
                      gamm.independent.t = NA, gamm.independent.pvalue = NA, boots.pvalues = NA,
                      correstimate = NA, corrp = NA, slope = NA))
  }
  
  # Formula construction
  covariates_string <- ""
  if (!is.na(covariates) && all(nchar(covariates) > 0)) {
    # Check if covariates have variance
    valid_covariates <- c()
    for(cv in covariates){
      if(cv %in% names(gam.data) && length(unique(na.omit(gam.data[[cv]]))) > 1){
        valid_covariates <- c(valid_covariates, cv)
      } else {
        warning(paste("Covariate", cv, "has no variance or is missing in subset, removing from formula for this run."))
      }
    }
    if(length(valid_covariates) > 0) {
      covariates_string <- paste(valid_covariates, collapse=" + ")
    } else {
      covariates_string <- "" # All covariates removed
    }
  }
  
  # Base terms (interest var and covariates)
  base_terms_model <- interest.indep.var
  if (nchar(covariates_string) > 0) {
    base_terms_model <- paste(covariates_string, interest.indep.var, sep=" + ")
    base_terms_null <- covariates_string
  } else {
    base_terms_null <- "1" # Null model is intercept only if no covariates
  }
  
  # Smooth term
  smooth_term_str <- ""
  smooth_term_present <- !is.na(smoothvar) && nchar(smoothvar) > 0 && (smoothvar %in% names(gam.data)) && (length(unique(na.omit(gam.data[[smoothvar]]))) >= 3) # k=3 needs at least 3 unique points
  
  if (smooth_term_present) {
    smooth_term_str <- paste0("s(", smoothvar, ", k=3, fx=T)") # fx=T means fixed df (not estimated)
  } else {
    if(!is.na(smoothvar) && nchar(smoothvar) > 0) warning(paste("Smooth variable", smoothvar, "not usable (missing, no variance, or too few unique points). Omitting smooth term."))
  }
  
  # Construct formulas for GAMM (nonlinear path)
  if (nchar(smooth_term_str) > 0) {
    form_gamm_nonlinear_model_str <- paste(dependentvar, "~", base_terms_model, "+", smooth_term_str)
    form_gamm_nonlinear_null_str  <- paste(dependentvar, "~", base_terms_null, "+", smooth_term_str)
    # For LM path, smoothvar is treated as linear predictor
    form_lm_linear_model_str <- paste(dependentvar, "~", base_terms_model, "+", smoothvar, "+ (1|src_subject_id)")
    form_lm_linear_null_str  <- paste(dependentvar, "~", base_terms_null, "+", smoothvar, "+ (1|src_subject_id)")
  } else { # No smooth term
    form_gamm_nonlinear_model_str <- paste(dependentvar, "~", base_terms_model) # gamm4 needs random effect
    form_gamm_nonlinear_null_str  <- paste(dependentvar, "~", base_terms_null)
    form_lm_linear_model_str <- paste(dependentvar, "~", base_terms_model, "+ (1|src_subject_id)")
    form_lm_linear_null_str  <- paste(dependentvar, "~", base_terms_null, "+ (1|src_subject_id)")
    if (base_terms_null == "1" && !grepl("1|src_subject_id", form_lm_linear_null_str) ) { # Ensure random effect for intercept-only null LM
      form_lm_linear_null_str <- paste(dependentvar, "~ 1 + (1|src_subject_id)")
    }
    if (base_terms_null == "1" && !grepl("1|src_subject_id", form_gamm_nonlinear_null_str) ) { # Ensure random effect for intercept-only null GAMM
      # gamm4 handles this by adding random effect specified in random= argument
    }
  }
  
  modelformula.nonlinear <- as.formula(form_gamm_nonlinear_model_str)
  modelformula.nonlinear.null <- as.formula(form_gamm_nonlinear_null_str)
  modelformula.linear <- as.formula(form_lm_linear_model_str)
  modelformula.linear.null <- as.formula(form_lm_linear_null_str)
  
  # Initialize results to NA
  partialRsq <- NA; bettermodel_out <- NA; gamm.independent.t <- NA; gamm.independent.pvalue <- NA;
  boots.pvalues <- NA; correstimate <- NA; corrp <- NA; slope <- NA
  
  # Fit models with tryCatch
  gamm.model <- tryCatch(gamm4(modelformula.nonlinear, random=~(1|src_subject_id), REML=T, data = gam.data), error=function(e){warning(paste("gamm4 full failed:",e$message));NULL})
  lm.model   <- tryCatch(lmer(modelformula.linear, REML = T, data=gam.data, control = lmerControl(check.nobs.vs.nRE = "ignore")), error=function(e){warning(paste("lmer full failed:",e$message));NULL})
  
  if(is.null(gamm.model) && is.null(lm.model)) {
    warning("Both GAMM and LM full models failed to fit.")
    return(data.frame( parcel = dependentvar, interest.indep.var = interest.indep.var, partialRsq, bettermodel=bettermodel_out,
                       gamm.independent.t, gamm.independent.pvalue, boots.pvalues, correstimate, corrp, slope))
  }
  
  gamm_AIC <- if(!is.null(gamm.model)) AIC(gamm.model$mer) else Inf
  lm_AIC   <- if(!is.null(lm.model)) AIC(lm.model) else Inf
  
  if (is.infinite(gamm_AIC) && is.infinite(lm_AIC)) { # Both failed
    bettermodel_out <- "None"
  } else if (gamm_AIC < lm_AIC) {
    bettermodel_out <- "GAMM"
  } else {
    bettermodel_out <- "LM"
  }
  
  if (bettermodel_out == "GAMM" && !is.null(gamm.model)) {
    gamm.results <- summary(gamm.model$gam)
    gamm.model.null <- tryCatch(gamm4(modelformula.nonlinear.null, random=~(1|src_subject_id), REML=T, data = gam.data), error=function(e){warning(paste("gamm4 null failed:",e$message));NULL})
    
    if (!is.null(gamm.model.null)) {
      sse.model <- sum((gamm.model$gam$y - gamm.model$gam$fitted.values)^2)
      sse.nullmodel <- sum((gamm.model.null$gam$y - gamm.model.null$gam$fitted.values)^2)
      if (sse.nullmodel > 0 && !is.na(sse.nullmodel) && !is.na(sse.model)) partialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel else partialRsq <- NA
      boots.pvalues <- pbootint(gamm.model) # pbootint itself handles its specific null model logic
    }
    
    # Parametric coefficient for interest.indep.var
    # Check if interest.indep.var is in the parametric part of the summary
    if (interest.indep.var %in% rownames(gamm.results$p.table)) {
      gamm.independent.t <- gamm.results$p.table[interest.indep.var, "t value"] # or "z value" depending on gam version/summary method
      gamm.independent.pvalue <- gamm.results$p.table[interest.indep.var, "Pr(>|t|)"] # or "Pr(>|z|)"
      slope <- gamm.results$p.table[interest.indep.var, "Estimate"]
    } else if (nrow(gamm.results$p.table) >= 2 && rownames(gamm.results$p.table)[2] == interest.indep.var) { 
      # Fallback to second row if it's the interest var (common for simple models)
      gamm.independent.t <- gamm.results$p.table[2, "t value"] 
      gamm.independent.pvalue <- gamm.results$p.table[2, "Pr(>|t|)"]
      slope <- gamm.results$p.table[2, "Estimate"]
    } else {
      warning(paste("Interest variable", interest.indep.var, "not found as expected in GAMM parametric summary."))
    }
    
    
  } else if (bettermodel_out == "LM" && !is.null(lm.model)) {
    lm.results <- summary(lm.model)
    lm.model.null <- tryCatch(lmer(modelformula.linear.null, REML = T, data=gam.data, control = lmerControl(check.nobs.vs.nRE = "ignore")), error=function(e){warning(paste("lmer null failed:",e$message));NULL})
    
    if (!is.null(lm.model.null)) {
      # Using PBmodcomp for LM p-value (model comparison)
      # Reduce nsim for PBmodcomp here for speed, e.g., 199
      pb_comp_res <- tryCatch(PBmodcomp(lm.model, lm.model.null, nsim=199), error = function(e) {warning(paste("PBmodcomp for LM failed:", e$message)); NULL})
      if (!is.null(pb_comp_res) && "PBtest" %in% rownames(pb_comp_res$test)) boots.pvalues <- pb_comp_res$test["PBtest","p.value"] else boots.pvalues <- NA
      
      sse.model <- sum(residuals(lm.model)^2)
      sse.nullmodel <- sum(residuals(lm.model.null)^2)
      if (sse.nullmodel > 0 && !is.na(sse.nullmodel) && !is.na(sse.model)) partialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel else partialRsq <- NA
    }
    
    coef_table <- lm.results$coefficients
    if (interest.indep.var %in% rownames(coef_table)) {
      gamm.independent.t <- coef_table[interest.indep.var, "t value"]
      # p-values from lmer summary (lmerTest) are Satterthwaite's approximation
      if ("Pr(>|t|)" %in% colnames(coef_table)) {
        gamm.independent.pvalue <- coef_table[interest.indep.var, "Pr(>|t|)"]
      } else { # if using base lme4, p-values are not standard output
        gamm.independent.pvalue <- NA 
      }
      slope <- coef_table[interest.indep.var, "Estimate"]
    } else if (nrow(coef_table) >= 2 && rownames(coef_table)[2] == interest.indep.var) {
      gamm.independent.t <- coef_table[2, "t value"]
      if ("Pr(>|t|)" %in% colnames(coef_table)) {
        gamm.independent.pvalue <- coef_table[2, "Pr(>|t|)"]
      } else { gamm.independent.pvalue <- NA }
      slope <- coef_table[2, "Estimate"]
    } else {
      warning(paste("Interest variable", interest.indep.var, "not found as expected in LM summary."))
    }
  } else {
    warning(paste("No better model determined or chosen model failed. Model was:", bettermodel_out))
  }
  
  # Partial correlation part - simplified, ensure residuals are valid
  res1 <- NA; res2 <- gam.data[[interest.indep.var]]
  if (bettermodel_out == "GAMM" && !is.null(gamm.model)) {
    # Residuals of DV after accounting for covariates and smooth (but not interest.indep.var)
    # This requires fitting a model Y ~ covariates + smooth
    # Your original varcorformula1 implies regressing out covariates and smooth from DV
    # And varcorformula2 implies regressing out covariates and smooth from IV (interest.indep.var)
    # Let's stick to your original logic for res1 (Y | Zs), res2 (raw X or X | Zs)
    # Original res1: Y ~ covariates + s(smoothvar)
    # Original res2: gam.data[[interest.indep.var]] (raw X)
    if (nchar(covariates_string) > 0 || nchar(smooth_term_str) > 0) {
      form_res1_str <- dependentvar %s+% " ~ "
      terms_res1 <- c(if(nchar(covariates_string)>0) covariates_string else NULL, 
                      if(nchar(smooth_term_str)>0) smooth_term_str else NULL)
      if(length(terms_res1) == 0) terms_res1 <- "1" # Intercept only if no covs/smooth
      form_res1_str <- paste0(dependentvar, " ~ ", paste(terms_res1, collapse=" + "))
      
      res1_model_gamm <- tryCatch(gamm4(as.formula(form_res1_str), random=~(1|src_subject_id), data=gam.data)$gam, error=function(e) NULL)
      if(!is.null(res1_model_gamm)) res1 <- residuals(res1_model_gamm)
    } else { # No covariates or smooth, res1 is just DV
      res1 <- gam.data[[dependentvar]]
    }
    
  } else if (bettermodel_out == "LM" && !is.null(lm.model)) {
    if (nchar(covariates_string) > 0 || (nchar(smoothvar) > 0 && smoothvar %in% names(gam.data))) { # smoothvar used linearly
      form_res1_str <- dependentvar %s+% " ~ "
      terms_res1_lm <- c(if(nchar(covariates_string)>0) covariates_string else NULL,
                         if(nchar(smoothvar)>0 && smoothvar %in% names(gam.data)) smoothvar else NULL)
      if(length(terms_res1_lm) == 0) terms_res1_lm <- "1"
      form_res1_str <- paste0(dependentvar, " ~ ", paste(terms_res1_lm, collapse=" + "), " + (1|src_subject_id)")
      
      res1_model_lm <- tryCatch(lmer(as.formula(form_res1_str), data=gam.data, control = lmerControl(check.nobs.vs.nRE = "ignore")), error=function(e) NULL)
      if(!is.null(res1_model_lm)) res1 <- residuals(res1_model_lm)
    } else {
      res1 <- gam.data[[dependentvar]]
    }
  }
  
  if (all(!is.na(res1)) && all(!is.na(res2)) && length(res1) > 5 && length(res2) > 5 && sd(res1, na.rm=T) > 1e-6 && sd(res2, na.rm=T) > 1e-6) {
    # Check normality for correlation method (ensure psych package is loaded for corr.test)
    if (!requireNamespace("psych", quietly = TRUE)) { stop("Package 'psych' needed for corr.test.") }
    # KS test can be sensitive to large N, consider shapiro.test for N < 5000
    res1_clean <- na.omit(res1)
    res2_clean <- na.omit(res2) # res2 should not have NAs at this stage based on earlier processing
    
    # Ensure enough distinct points for KS test
    if(length(unique(res1_clean)) > 1 && length(unique(res2_clean)) > 1 && length(res1_clean) > 3 && length(res2_clean) > 3){
      # ks.test is for continuous variables.
      # Defaulting to spearman if unsure or tests fail.
      corrmethod <- "spearman" # Default
      res1.normtest_p <- tryCatch(ks.test(res1_clean, "pnorm", mean(res1_clean), sd(res1_clean))$p.value, error = function(e) 0)
      res2.normtest_p <- tryCatch(ks.test(res2_clean, "pnorm", mean(res2_clean), sd(res2_clean))$p.value, error = function(e) 0)
      
      if (res1.normtest_p > 0.01 && res2.normtest_p > 0.01) {
        corrmethod <- "pearson"
      }
      
      PCorr_Test <- tryCatch(psych::corr.test(res1_clean, res2_clean, method=corrmethod, adjust="none"), error=function(e) {warning(paste("corr.test failed:", e$message)); NULL})
      if(!is.null(PCorr_Test)){
        correstimate <- as.numeric(PCorr_Test$r)
        corrp <- as.numeric(PCorr_Test$p)
      }
    } else {
      warning("Not enough distinct data points in residuals for normality test or correlation.")
    }
  } else {
    warning("Residuals for correlation are problematic (NA, no variance, or too few).")
  }
  
  stats.results <- data.frame(
    parcel = dependentvar,
    interest.indep.var = interest.indep.var,
    partialRsq = partialRsq,
    bettermodel = bettermodel_out,
    gamm.independent.t = gamm.independent.t,
    gamm.independent.pvalue = gamm.independent.pvalue,
    boots.pvalues = boots.pvalues, # This is the model comparison p-value (from pbootint or PBmodcomp for LM)
    correstimate = correstimate,
    corrp = corrp,
    slope = slope
  )
  return(stats.results)
}


run_robustness_analysis_parallel <- function(dependentvar, full_data, smoothvar, interest.indep.var, covariates,
                                             stratify_by_var = "Sex", 
                                             n_splithalf = 101, n_permutations = 1000,
                                             statistic_to_permute = "partialRsq",
                                             seed = 123,
                                             num_cores = max(1, detectCores() - 1)) {
  main_seed <- seed # Save main seed for reference
  set.seed(main_seed) 
  
  # 1. Calculate observed median statistic on original data
  cat("Calculating observed median statistic using seed:", main_seed, "for split-half sampling...\n")
  if(!is.numeric(full_data[[dependentvar]])) {
    stop(paste("Dependent variable", dependentvar, "is not numeric."))
  }
  if (!is.null(stratify_by_var) && !(stratify_by_var %in% names(full_data))) {
    warning(paste("Stratification variable", stratify_by_var, "not found in full_data. Proceeding without stratification for observed statistic."))
    stratify_by_var_obs <- NULL
  } else {
    stratify_by_var_obs <- stratify_by_var
  }
  # Check for NAs in stratification variable if used for observed
  if (!is.null(stratify_by_var_obs) && any(is.na(full_data[[stratify_by_var_obs]]))) {
    warning(paste("Stratification variable '", stratify_by_var_obs, "' contains NAs. Rows with NA in this variable will affect stratification."))
  }
  
  # Define the split-half function to be used for observed and within permutation_task
  # This function will be defined inside permutation_task for workers, but we also need it here for observed.
  .get_median_stat_from_splithalf_local <- function(data_for_splithalf, dep_var, sm_var, int_indep_var, covs, 
                                                    n_sh_internal, stat_col_internal, strat_var_internal, current_seed) {
    set.seed(current_seed) # Ensure reproducibility for this set of split-halves
    split_half_stats_collected <- numeric(0)
    
    if (!is.null(strat_var_internal) && !(strat_var_internal %in% names(data_for_splithalf))) {
      warning(paste("Stratification var", strat_var_internal, "not in data. No stratification."))
      strat_var_internal <- NULL
    }
    if (!is.null(strat_var_internal) && length(stats::na.omit(unique(data_for_splithalf[[strat_var_internal]]))) < 1) { # or <2 if expecting multiple groups
      warning(paste("Stratification var", strat_var_internal, "has <1 unique non-NA levels. No stratification."))
      strat_var_internal <- NULL
    }
    
    for (i_sh in 1:n_sh_internal) {
      if (nrow(data_for_splithalf) < 40) { next } # Skip if data too small
      
      data_s1 <- NULL; data_s2 <- NULL
      
      if (!is.null(strat_var_internal)) {
        # Stratified sampling logic (simplified from previous, using dplyr if available)
        # Ensure strat_var_internal is a factor
        data_for_splithalf[[strat_var_internal]] <- as.factor(data_for_splithalf[[strat_var_internal]])
        
        # Try using dplyr for cleaner stratified sampling if available
        if (requireNamespace("dplyr", quietly = TRUE)) {
          # dplyr approach: group by strata, then sample half from each.
          # This needs careful handling of small strata.
          # A simpler approach for balanced splits:
          temp_data_strat <- data_for_splithalf %>% 
            dplyr::mutate(.original_row_id_strat = dplyr::row_number()) %>%
            dplyr::group_by_at(dplyr::vars(one_of(strat_var_internal))) %>%
            dplyr::mutate(.group_row_id = dplyr::row_number(),
                          .group_size = dplyr::n())
          
          # Try to assign roughly half of each stratum to each split
          # This is tricky to get perfectly balanced 50/50 splits while maintaining strata proportions
          # A common method is to sample n/2 from each stratum for one half.
          
          # Simplification: create an index for each half within strata
          # This is a basic attempt; more robust methods exist (e.g. 'splitstackshape' package)
          
          # Fallback to simpler non-dplyr way if dplyr is complex here
          # Manual stratified sampling:
          all_indices <- 1:nrow(data_for_splithalf)
          s1_indices_list <- list()
          s2_indices_list <- list()
          valid_stratification = TRUE
          
          unique_strata_vals <- unique(stats::na.omit(data_for_splithalf[[strat_var_internal]]))
          
          if(length(unique_strata_vals) == 0) { # No valid strata levels
            valid_stratification = FALSE
          } else {
            for(str_val in unique_strata_vals){
              str_indices <- all_indices[which(data_for_splithalf[[strat_var_internal]] == str_val & !is.na(data_for_splithalf[[strat_var_internal]]))]
              n_s <- length(str_indices)
              if(n_s == 0) next
              if(n_s == 1){ # Assign singletons randomly or consistently
                if(runif(1) < 0.5) s1_indices_list[[as.character(str_val)]] <- str_indices
                else s2_indices_list[[as.character(str_val)]] <- str_indices
                next
              }
              half_n_s <- sample(str_indices, floor(n_s/2))
              s1_indices_list[[as.character(str_val)]] <- half_n_s
              s2_indices_list[[as.character(str_val)]] <- setdiff(str_indices, half_n_s)
            }
          }
          
          s1_final_indices <- unlist(s1_indices_list)
          s2_final_indices <- unlist(s2_indices_list)
          
          if(length(s1_final_indices) < 10 || length(s2_final_indices) < 10 || !valid_stratification){ # Min N for a split
            warning("Stratified sampling resulted in very small splits or failed. Falling back to non-stratified for this SH iter.")
            strat_var_internal <- NULL # Fallback for this SH iter
          } else {
            data_s1 <- data_for_splithalf[s1_final_indices, ]
            data_s2 <- data_for_splithalf[s2_final_indices, ]
          }
        } else { # dplyr not available
          warning("dplyr package not available for robust stratified sampling. Falling back to non-stratified.")
          strat_var_internal <- NULL
        }
      } # End if !is.null(strat_var_internal)
      
      if (is.null(strat_var_internal) || is.null(data_s1) || is.null(data_s2)) { # Non-stratified or fallback
        sample_indices_all <- 1:nrow(data_for_splithalf)
        size_half1_overall <- floor(nrow(data_for_splithalf) / 2)
        if (size_half1_overall < 10 || (nrow(data_for_splithalf) - size_half1_overall) < 10) {next}
        original_row_ids_half1 <- sample(sample_indices_all, size = size_half1_overall)
        data_s1 <- data_for_splithalf[original_row_ids_half1, ]
        data_s2 <- data_for_splithalf[setdiff(sample_indices_all, original_row_ids_half1), ]
      }
      
      if(is.null(data_s1) || is.null(data_s2) || nrow(data_s1) < 20 || nrow(data_s2) < 20){next}
      
      res1 <- gamm.smooth.predict.interaction(dependentvar = dep_var, data_input = data_s1, smoothvar = sm_var, interest.indep.var = int_indep_var, covariates = covs)
      if (!is.null(res1) && nrow(res1) > 0 && !is.na(res1[[stat_col_internal]])) { split_half_stats_collected <- c(split_half_stats_collected, res1[[stat_col_internal]])}
      
      res2 <- gamm.smooth.predict.interaction(dependentvar = dep_var, data_input = data_s2, smoothvar = sm_var, interest.indep.var = int_indep_var, covariates = covs)
      if (!is.null(res2) && nrow(res2) > 0 && !is.na(res2[[stat_col_internal]])) { split_half_stats_collected <- c(split_half_stats_collected, res2[[stat_col_internal]])}
    }
    if (length(split_half_stats_collected) > 0) { return(median(split_half_stats_collected, na.rm = TRUE))
    } else { return(NA_real_) }
  } # End of .get_median_stat_from_splithalf_local definition
  
  # Generate a seed for the observed statistic's split-half process
  observed_splithalf_seed <- sample.int(1e6, 1) 
  observed_median_stat <- .get_median_stat_from_splithalf_local(
    data_for_splithalf = full_data, 
    dep_var = dependentvar, 
    sm_var = smoothvar, 
    int_indep_var = interest.indep.var, 
    covs = covariates, 
    n_sh_internal = n_splithalf, 
    stat_col_internal = statistic_to_permute,
    strat_var_internal = stratify_by_var_obs, # Use the potentially nulled version
    current_seed = observed_splithalf_seed
  )
  cat("Observed median", statistic_to_permute, ":", observed_median_stat, "\n")
  
  if (is.na(observed_median_stat)) {
    stop("Failed to calculate observed median statistic. Check warnings.")
  }
  
  # 2. Permutation loop (parallelized)
  y_original <- full_data[[dependentvar]]
  
  cat("Starting permutations with", num_cores, "cores...\n")
  cl <- makeCluster(num_cores)
  
  # Export necessary variables and functions to the cluster
  clusterExport(cl, c("pbootint", "gamm.smooth.predict.interaction"), envir = .GlobalEnv)
  # If gamm.smooth.predict.interaction or pbootint uses any custom helper functions not from packages, export them too.
  # e.g. clusterExport(cl, "my_helper_function_name")
  
  # Load packages on workers
  clusterEvalQ(cl, {
    library(mgcv)
    library(lme4)
    library(gamm4)
    library(pbkrtest)
    library(tidyverse) # Includes dplyr, stringr
    library(psych)
    # Define .get_median_stat_from_splithalf_local on workers (copied from above)
    # This is crucial because parLapply needs this function definition on each worker.
    # (Definition is identical to the one above)
    assign(".get_median_stat_from_splithalf_on_worker", 
           eval(body(.get_median_stat_from_splithalf_local)), 
           envir = .GlobalEnv)
    # Make sure the arguments are correctly captured for the worker version
    formals(.get_median_stat_from_splithalf_on_worker) <- formals(.get_median_stat_from_splithalf_local)
    
  })
  
  # Generate seeds for each permutation's split-half process
  # These seeds are for the random sampling within each permutation's 101 split-halves
  permutation_splithalf_seeds <- sample.int(1e6, n_permutations, replace = FALSE) # Unique seed for each perm's SH series
  
  permutation_task_for_parLapply <- function(perm_index, .full_data_task, .y_original_task, .dep_var_task, 
                                             .sm_var_task, .int_indep_var_task, .covs_task, .n_sh_task, 
                                             .stat_col_task, .strat_var_task, .iter_seed) {
    # Note: .iter_seed is the seed for THIS permutation's split-half sequence.
    # The permutation of Y itself is done based on the main R session's RNG state before parLapply.
    # To make Y permutation itself reproducible across parallel runs if needed,
    # one would set seed for `sample(y_original)` using a seed derived from perm_index.
    # However, standard permutation tests usually rely on the overall RNG sequence.
    
    # The actual permutation (shuffling of Y) happens here, inside the worker for this specific iteration.
    # This ensures each worker shuffles Y independently based on its own RNG state,
    # which is good for permutation tests. If we want the *exact same set* of permuted Ys
    # across different runs of the entire script (given the same main seed),
    # we'd need to pre-generate permuted Ys or use seeds for shuffling.
    # For now, standard `sample()` behavior is fine.
    
    current_perm_data <- .full_data_task
    current_perm_data[[.dep_var_task]] <- sample(.y_original_task) 
    
    # Stratification variable for permuted data run
    strat_var_perm_run <- .strat_var_task
    if (!is.null(strat_var_perm_run) && !(strat_var_perm_run %in% names(current_perm_data))) {
      strat_var_perm_run <- NULL
    }
    if (!is.null(strat_var_perm_run) && any(is.na(current_perm_data[[strat_var_perm_run]]))) {
      # Silently proceed, NAs handled by split-half func
    }
    
    
    # Call the split-half function (which should now be defined on the worker)
    # Pass .iter_seed as current_seed to .get_median_stat_from_splithalf_on_worker
    median_stat <- .GlobalEnv$.get_median_stat_from_splithalf_on_worker( 
      data_for_splithalf = current_perm_data, 
      dep_var = .dep_var_task, 
      sm_var = .sm_var_task, 
      int_indep_var = .int_indep_var_task, 
      covs = .covs_task, 
      n_sh_internal = .n_sh_task, 
      stat_col_internal = .stat_col_task,
      strat_var_internal = strat_var_perm_run,
      current_seed = .iter_seed 
    )
    return(median_stat)
  }
  
  permuted_median_stats_list <- parLapply(cl, 1:n_permutations, 
                                          permutation_task_for_parLapply,
                                          .full_data_task = full_data, 
                                          .y_original_task = y_original,
                                          .dep_var_task = dependentvar,
                                          .sm_var_task = smoothvar,
                                          .int_indep_var_task = interest.indep.var,
                                          .covs_task = covariates,
                                          .n_sh_task = n_splithalf,
                                          .stat_col_task = statistic_to_permute,
                                          .strat_var_task = stratify_by_var,
                                          .iter_seed = permutation_splithalf_seeds # Pass the vector of seeds
  )
  
  stopCluster(cl)
  cat("\nPermutations finished.\n")
  
  permuted_median_stats <- unlist(permuted_median_stats_list)
  permuted_median_stats_valid <- permuted_median_stats[!is.na(permuted_median_stats)]
  
  if(length(permuted_median_stats_valid) < n_permutations * 0.8 && n_permutations > 0) {
    warning(paste0("More than 20% of permutations resulted in NA. ",
                   length(permuted_median_stats_valid), " valid out of ", n_permutations, "."))
  }
  if(length(permuted_median_stats_valid) == 0 && n_permutations > 0) {
    permutation_p_value <- NA
    warning("All permutations resulted in NA. Cannot calculate p-value.")
  } else if (n_permutations == 0) {
    permutation_p_value <- NA # No permutations run
  } else {
    if (statistic_to_permute %in% c("gamm.independent.t", "slope")) {
      permutation_p_value <- (sum(abs(permuted_median_stats_valid) >= abs(observed_median_stat), na.rm = TRUE) + 1) / (length(permuted_median_stats_valid) + 1)
    } else { # For Rsq (non-negative)
      permutation_p_value <- (sum(permuted_median_stats_valid >= observed_median_stat, na.rm = TRUE) + 1) / (length(permuted_median_stats_valid) + 1)
    }
  }
  
  return(list(
    dependent_var = dependentvar,
    interest_var = interest.indep.var,
    statistic_tested = statistic_to_permute,
    observed_median_statistic = observed_median_stat,
    permuted_median_statistics_dist = permuted_median_stats_valid,
    permutation_p_value = permutation_p_value,
    n_permutations_requested = n_permutations,
    n_permutations_completed_valid = length(permuted_median_stats_valid),
    main_seed_used = main_seed,
    observed_splithalf_seed = observed_splithalf_seed,
    permutation_splithalf_seeds_summary = summary(permutation_splithalf_seeds) # Just to confirm they were generated
  ))
}