library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(lme4)
library(gamm4)
library(pbkrtest)

pbootint <- function(modelobj, int_var = NA) {
  numsims <- 1000
  set.seed(925)
  
  df <- modelobj$gam$model
  thisResp <- as.character(modelobj$gam$terms[[2]])
  f1 <- modelobj$gam$formula
  theseVars <- unique(attr(terms(f1), "term.labels"))
  
  print("Original formula (f1):")
  print(f1)
  
  if (sum(str_detect(theseVars, "by ="))==1){
    int_var <- str_split_i(theseVars[1], "by = ", 2)
    int_var <- str_split_i(int_var, ", ", 1)
    f2 <- reformulate(c(theseVars[2:(length(theseVars))], int_var),response = thisResp)
  }else{
    f2 <- reformulate(theseVars[2:(length(theseVars))],response = thisResp)
  }
  
  print("Modified formula (f2):")
  print(f2)
  
  g1 <- gam(f1, data = df)
  g2 <- gam(f2, data = df)
  
  mat1 <- model.matrix(g1)
  mat2 <- model.matrix(g2)
  
  ID <- df$ID
  y <- df[, thisResp]
  
  m1 <- lmer(y ~ -1 + mat1 + (1 | ID))
  m2 <- lmer(y ~ -1 + mat2 + (1 | ID))
  
  refdist <- PBrefdist(m1, m2, nsim = numsims)
  pb <- PBmodcomp(m1, m2, ref = refdist)
  int_pval <- pb$test["PBtest", "p.value"]
  
  return(int_pval)
}
# pbootint <- function(modelobj, int_var=NA){
#   numsims <- 1000
#   set.seed(925)
#   df <- modelobj$gam$model
#   thisResp <- as.character(modelobj$gam$terms[[2]])
#   f1 <- modelobj$gam$formula
#   theseVars <- attr(terms(f1),"term.labels")
#   if (!is.na(int_var)){
#     indexcomma <- gregexpr(",", theseVars[1])[[1]]
#     addsmooth <- theseVars[1]
#     addsmooth <- paste0(substr(addsmooth, 1, indexcomma[1]), substr(addsmooth, indexcomma[2]+1, nchar(addsmooth)))
#     theseVars <-  c(theseVars, addsmooth)
#   }
#   # if (sum(str_detect(theseVars, "by ="))==1){
#   #   int_var <- str_split_i(theseVars[1], "by = ", 2)
#   #   int_var <- str_split_i(int_var, ", ", 1)
#   #   f2 <- reformulate(c(theseVars[2:(length(theseVars))], int_var),response = thisResp)
#   # }else{
#   #   f2 <- reformulate(theseVars[2:(length(theseVars))],response = thisResp)
#   # }
#   if (any(str_detect(theseVars, "by ="))) {
#     f2_terms <- theseVars[!str_detect(theseVars, "by =") & !str_detect(theseVars, "s\\(Age_year")]
#     f2 <- reformulate(c(f2_terms, "s(Age_year, k = 3, fx = TRUE)"), response = thisResp) # 添加主要平滑项
#   } else {
#     f2 <- reformulate(theseVars, response = thisResp) # 如果没有交互项，保留其他项
#   }
# 
# 
#   g1 <- gam(f1,data = df)
#   g2 <- gam(f2,data = df)
# 
#   mat1 <- model.matrix(g1)
#   mat2 <- model.matrix(g2)
# 
#   ID<-df$ID
#   y <- df[,thisResp]
# 
#   m1 <- lmer(y ~ -1 + mat1 + (1|ID))
#   m2 <- lmer(y ~ -1 + mat2 + (1|ID))
#   refdist <- PBrefdist(m1, m2, nsim=numsims)
#   pb <- PBmodcomp(m1, m2, ref = refdist)
#   int_pval <- pb$test["PBtest","p.value"]
# 
#   return(int_pval)
# }


# pbootint <- function(modelobj, int_var = NA) {
#   numsims <- 1000
#   set.seed(925)
#   
#   # 获取数据和模型信息
#   df <- modelobj$gam$model
#   thisResp <- as.character(modelobj$gam$terms[[2]])
#   f1 <- modelobj$gam$formula
#   theseVars <- unique(attr(terms(f1), "term.labels")) # 去重
#   
#   
#   # 如果指定了 int_var，添加额外平滑项（暂时不使用这个逻辑，因为 int_var 在公式中已定义）
#   if (!is.na(int_var)) {
#     indexcomma <- gregexpr(",", theseVars[1])[[1]]
#     addsmooth <- theseVars[1]
#     addsmooth <- paste0(
#       substr(addsmooth, 1, indexcomma[1]),
#       substr(addsmooth, indexcomma[2] + 1, nchar(addsmooth))
#     )
#     theseVars <- unique(c(theseVars, addsmooth))
#   }
#   
#   
#   # 构建两个模型
#   g1 <- gam(f1, data = df)
#   g2 <- gam(f2, data = df)
#   
#   # 计算模型矩阵
#   mat1 <- model.matrix(g1)
#   mat2 <- model.matrix(g2)
#   
#   # 获取被试 ID 和响应变量
#   ID <- df$ID
#   y <- df[, thisResp]
#   
#   # 使用线性混合模型拟合
#   m1 <- lmer(y ~ -1 + mat1 + (1 | ID))
#   m2 <- lmer(y ~ -1 + mat2 + (1 | ID))
#   
#   # 使用 PBrefdist 和 PBmodcomp 计算显著性
#   refdist <- PBrefdist(m1, m2, nsim = numsims)
#   pb <- PBmodcomp(m1, m2, ref = refdist)
#   int_pval <- pb$test["PBtest", "p.value"]
#   
#   # 返回交互效应的 p 值
#   return(int_pval)
# }



# gamm.varyingcoefficients <- function(dependentvar, dataname, smooth_var, int_var, covariates, knots, set_fx = FALSE, increments, draws, return_posterior_coefficients = FALSE){
# 
#   #Set parameters
#   npd <- as.numeric(draws) #number of draws from the posterior distribution
#   np <- as.numeric(increments) #number of smooth_var increments to get derivatives at
#   UNCONDITIONAL <- FALSE #should we account for uncertainty when estimating smoothness parameters?
# 
#   #Fit the gam
#   gam.data <- get(dataname)
#   smoothmin <- min(gam.data[ , smooth_var])
#   smoothmax <- max(gam.data[ , smooth_var])
#   cognition<- as.numeric(gam.data[[int_var]])
#   outlierindx1<-which(cognition<mean(cognition,na.rm = T)-3*sd(cognition,na.rm = T) | cognition>mean(cognition,na.rm = T)+3*sd(cognition,na.rm = T))
#   gam.data[outlierindx1 ,int_var]<-NA
#   NonNANIndex <- which(!is.na(gam.data[ ,int_var]) & !is.na(gam.data[ ,dependentvar]))
# 
#   gam.data <- gam.data[NonNANIndex,]
#   parcel <- dependentvar
#   modelformula <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s) + s(%2$s, by=%5$s, k=%3$s, fx=%4$s) + %6$s", dependentvar, smooth_var, knots, set_fx, int_var,covariates))
#   gamm.model <- gamm4(modelformula, random=~(1|ID), REML=TRUE, data = gam.data)
#   # Fit null model for ANOVA test
#   modelformula.null <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s) + %5$s + %6$s", dependentvar, smooth_var, knots, set_fx, int_var,covariates))
#   gamm.model.null <- gamm4(modelformula.null, random = ~(1|ID), REML = TRUE, data = gam.data)
# 
#   # Calculate ANOVA p-value
#   anova_results <- anova(gamm.model$mer, gamm.model.null$mer, test = "Chisq")
#   anova.int.pvalue <- anova_results$`Pr(>Chisq)`[2]
#   boots.pvalues <- pbootint(gamm.model)
# 
#   # Partial R² as effect size
#   sse.model <- sum((gamm.model$gam$y - gamm.model$gam$fitted.values)^2)
#   sse.nullmodel <- sum((gamm.model.null$gam$y - gamm.model.null$gam$fitted.values)^2)
#   IntpartialRsq <- (sse.nullmodel - sse.model) / sse.nullmodel
# 
#   # Z-value for ANOVA p-value
#   anova.int.zvalue <- qnorm(anova.int.pvalue / 2, lower.tail = FALSE)
#   full.results <- cbind(int_var, dependentvar, anova.int.pvalue, IntpartialRsq, boots.pvalues, anova.int.zvalue)
# 
#   #Extract gam input data
#   df <- gamm.model$gam$model #extract the data used to build the gam, i.e., a df of y + predictor values
# 
#   #Create a prediction data frame, used to estimate (posterior) model slopes (varying covariate coefficients)
#   theseVars <- attr(gamm.model$gam$terms,"term.labels")
#   varClasses <- attr(gamm.model$gam$terms,"dataClasses")
# 
#   #prediction df for int_var P10th at np smooth_var increments
#   pred.low <- data.frame(init = rep(0,np))
#   for (v in c(1:length(theseVars))) {
#     thisVar <- theseVars[[v]]
#     thisClass <- varClasses[thisVar]
#     if (thisVar == int_var) {
#       #pred.low[,int_var] <- (min(df[,int_var],na.rm = T))
#       pred.low[,int_var] <- (quantile(df[,int_var], probs=0.1, na.rm = T))
#     } else if (thisVar == smooth_var) {
#       pred.low[,smooth_var] = seq(smoothmin,smoothmax, length.out = np)
#     } else {
#       tab.tmp<-table(df[,thisVar])
#       levelact<-which.max(tab.tmp)
#       switch (thisClass,
#               "numeric" = {pred.low[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
#               "factor" = {pred.low[,thisVar] = levels(df[,thisVar])[[levelact]]}, #make predictions based on first level of factor
#               "ordered" = {pred.low[,thisVar] = levels(df[,thisVar])[[levelact]]} #make predictions based on first level of ordinal variable
#       )}}
# 
#   #prediction df for int_var P90th at np smooth_var increments
#   pred.high <- data.frame(init = rep(0,np))
#   for (v in c(1:length(theseVars))) {
#     thisVar <- theseVars[[v]]
#     thisClass <- varClasses[thisVar]
#     if (thisVar == int_var) {
#       #pred.high[,int_var] <- (max(df[,int_var],na.rm = T))
#       pred.high[,int_var] <- (quantile(df[,int_var], probs=0.9, na.rm = T))
#     } else if (thisVar == smooth_var) {
#       pred.high[,smooth_var] = seq(smoothmin,smoothmax, length.out = np)
#     } else {
#       tab.tmp<-table(df[,thisVar])
#       levelact<-which.max(tab.tmp)
#       switch (thisClass,
#               "numeric" = {pred.high[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
#               "factor" = {pred.high[,thisVar] = levels(df[,thisVar])[[levelact]]}, #make predictions based on first level of factor
#               "ordered" = {pred.high[,thisVar] = levels(df[,thisVar])[[levelact]]} #make predictions based on first level of ordinal variable
#       )}}
# 
#   pred <- rbind(pred.low, pred.high) #complete pred df
#   pred <- pred %>% dplyr::select(-init)
# 
#   #Get effects (slopes) along the smooth function for the true model
#   if(return_posterior_coefficients == FALSE){
#     #varying coefficient slopes
#     predicted.values <- fitted_values(object = gamm.model$gam, data = pred) #predict y at min and mix int_var along the smooth function
#     predicted.values <- predicted.values %>% dplyr::select(fitted, all_of(smooth_var), all_of(int_var))
#     colnames(predicted.values) <- c("fitted", "smooth.var", "int.var")
#     predicted.values$smooth.var <- round(predicted.values$smooth.var, 3)
# 
#     varyingcoeff.slopes <-  predicted.values %>% #calculate the effect of int_var on y  (slope; delta y/delta int_var) along the smooth function
#       group_by(smooth.var) %>%
#       do(slope = diff(.$fitted)/diff(.$int.var)) %>%
#       unnest(cols = c(slope))
#     colnames(varyingcoeff.slopes) <- c(smooth_var, "slope")
#   }
# 
#   varyingcoeff.slopes <- NULL
#   varyingcoeff.slopes.CI <- NULL
# 
#   if (return_posterior_coefficients == FALSE) {
#     # Varying coefficient slopes
#     predicted.values <- fitted_values(object = gamm.model$gam, data = pred)
#     predicted.values <- predicted.values %>% dplyr::select(.fitted, all_of(smooth_var), all_of(int_var))
#     colnames(predicted.values) <- c(".fitted", "smooth.var", "int.var")
#     predicted.values$smooth.var <- round(predicted.values$smooth.var, 3)
# 
#     varyingcoeff.slopes <- predicted.values %>%
#       group_by(smooth.var) %>%
#       summarise(slope = diff(.fitted) / diff(int.var))
#     colnames(varyingcoeff.slopes) <- c(smooth_var, "slope")
#   }
# 
#   # Posterior distribution of slopes
#   if (return_posterior_coefficients == TRUE) {
#     Vb <- vcov(gamm.model$gam, unconditional = UNCONDITIONAL)
#     sims <- MASS::mvrnorm(npd, mu = coef(gamm.model$gam), Sigma = Vb)
#     X0 <- predict(gamm.model$gam, newdata = pred, type = "lpmatrix")
#     predicted.values.posterior <- X0 %*% t(sims)
# 
#     predicted.values.posterior <- as.data.frame(predicted.values.posterior)
#     colnames(predicted.values.posterior) <- sprintf("draw%s", seq_len(npd))
#     predicted.values.posterior <- cbind(pred[, smooth_var], pred[, int_var], predicted.values.posterior)
#     colnames(predicted.values.posterior)[1:2] <- c("smooth.var", "int.var")
#     predicted.values.posterior$smooth.var <- round(predicted.values.posterior$smooth.var, 3)
# 
#     varyingcoeff.slopes.CI <- predicted.values.posterior %>%
#       pivot_longer(cols = starts_with("draw"), names_to = "draw", values_to = "posterior.fitted") %>%
#       mutate(int.var = ifelse(int.var == quantile(df[, int_var], 0.1), "low", "high")) %>%
#       pivot_wider(names_from = int.var, values_from = posterior.fitted) %>%
#       mutate(slope = (high - low) / (quantile(df[, int_var], 0.9) - quantile(df[, int_var], 0.1))) %>%
#       select(smooth.var, draw, slope)
#     varyingcoeff.slopes.CI <- cbind(as.character(dependentvar), varyingcoeff.slopes.CI)
#     colnames(varyingcoeff.slopes.CI) <- c("label", smooth_var, "draw", "slope")
#   }
# 
#   data.results1 <- list(full.results, gamm.model, varyingcoeff.slopes, pred)
#   data.results2 <- list(full.results,  gamm.model, varyingcoeff.slopes.CI, pred)
# 
#   if(return_posterior_coefficients == FALSE)
#     return(data.results1)
#   if(return_posterior_coefficients == TRUE)
#     return(data.results2)
# 
# }

gamm.varyingcoefficients <- function(dependentvar, dataname, smooth_var, int_var, covariates, bestmodel, knots, set_fx = FALSE, increments, draws, return_posterior_coefficients = FALSE){
  
  cat("\n--- Running Analysis for:", dependentvar, " & ", int_var, "---\n")
  
  # Set parameters
  npd <- as.numeric(draws)
  np <- as.numeric(increments)
  UNCONDITIONAL <- FALSE
  
  # Fit the gam
  gam.data.full <- get(dataname)
  
  # --- 调试点 1: 检查进入函数的原始数据 ---
  cat("1. Dimensions of initial data '", dataname, "': ", nrow(gam.data.full), "rows, ", ncol(gam.data.full), "cols\n", sep="")
  cat("   Unique IDs in initial data: ", length(unique(gam.data.full$ID)), "\n")
  
  # 内部数据处理
  cognition <- as.numeric(gam.data.full[[int_var]])
  outlierindx1 <- which(cognition < mean(cognition, na.rm = TRUE) - 3 * sd(cognition, na.rm = TRUE) | cognition > mean(cognition, na.rm = TRUE) + 3 * sd(cognition, na.rm = TRUE))
  if(length(outlierindx1) > 0) gam.data.full[outlierindx1, int_var] <- NA
  
  # 关键的NA筛选
  NonNANIndex <- which(!is.na(gam.data.full[, int_var]) & !is.na(gam.data.full[, dependentvar]))
  gam.data <- gam.data.full[NonNANIndex, ]
  
  # --- 调试点 2: 检查筛选后的数据 ---
  cat("2. Filtering for non-NA in '", dependentvar, "' and '", int_var, "'\n", sep="")
  cat("   Dimensions of data AFTER filtering NAs: ", nrow(gam.data), "rows, ", ncol(gam.data), "cols\n", sep="")
  cat("   Unique IDs in data AFTER filtering NAs: ", length(unique(gam.data$ID)), "\n")
  
  # 检查是否还有足够的数据
  if (nrow(gam.data) < 10 || length(unique(gam.data$ID)) <= 1) {
    cat("   !!! ERROR: Not enough data or unique IDs to run the model. Skipping this combination. !!!\n")
    # 返回一个空结果或者NULL，避免程序崩溃
    return(NULL)
  }
  smoothmin <- min(gam.data[ , smooth_var])
  smoothmax <- max(gam.data[ , smooth_var])
  cognition<- as.numeric(gam.data[[int_var]]) 
  outlierindx1<-which(cognition<mean(cognition,na.rm = T)-3*sd(cognition,na.rm = T) | cognition>mean(cognition,na.rm = T)+3*sd(cognition,na.rm = T))
  gam.data[outlierindx1 ,int_var]<-NA
  NonNANIndex <- which(!is.na(gam.data[ ,int_var]) & !is.na(gam.data[ ,dependentvar]))
  
  gam.data <- gam.data[NonNANIndex,]
  parcel <- dependentvar
  
  if (bestmodel=="GAMM"){
    modelformula <- as.formula(sprintf("%1$s ~ s(%2$s, by=%5$s, k=%3$s, fx=%4$s) + s(%2$s, k=%3$s, fx=%4$s) + %6$s", dependentvar, smooth_var, knots, set_fx, int_var,covariates))
    gamm.model <- gamm4(modelformula, random=~(1|ID), REML=TRUE, data = gam.data)
    # Fit null model for ANOVA test
    modelformula.null <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s) + %5$s + %6$s", dependentvar, smooth_var, knots, set_fx, int_var,covariates))
    gamm.model.null <- gamm4(modelformula.null, random = ~(1|ID), REML = TRUE, data = gam.data)
    anova_results <- anova(gamm.model$mer, gamm.model.null$mer, test = "Chisq")
    anova.int.pvalue <- anova_results$`Pr(>Chisq)`[2]
    if (return_posterior_coefficients == T){
      boots.pvalues <- pbootint(gamm.model, int_var)
    }else{
      boots.pvalues <- NA
    }
    
    # Partial R² as effect size
    sse.model <- sum((gamm.model$gam$y - gamm.model$gam$fitted.values)^2)
    sse.nullmodel <- sum((gamm.model.null$gam$y - gamm.model.null$gam$fitted.values)^2)
    IntpartialRsq <- (sse.nullmodel - sse.model) / sse.nullmodel
    
    stats.results <- cbind(int_var, dependentvar, anova.int.pvalue, IntpartialRsq, boots.pvalues)
    
    #Extract gam input data
    df <- gam.data
    #Create a prediction data frame, used to estimate (posterior) model slopes (varying covariate coefficients)
    theseVars <- attr(gamm.model$gam$terms,"term.labels") 
    varClasses <- unlist(lapply(theseVars, function(x) class(gam.data[[x]])))
    
    model.tmp <- gamm.model$gam
    
  }else if (bestmodel=="LM"){
    modelformula <- as.formula(sprintf("%s ~ %s * %s + %s + (1|ID)", dependentvar, int_var, smooth_var, covariates))
    lm.model <- lmer(modelformula, REML = T, data=gam.data)
    # Fit null model for ANOVA test
    modelformula.null <- as.formula(sprintf("%s ~ %s + %s + %s + (1|ID)", dependentvar, int_var, covariates, smooth_var))
    lm.model.null <- lmer(modelformula.null, REML = T, data=gam.data)
    anova_results <- anova(lm.model, lm.model.null, test = "Chisq")
    anova.int.pvalue <- anova_results$`Pr(>Chisq)`[2]
    set.seed(925)
    if (return_posterior_coefficients == T){
      boots.pvalues <- PBmodcomp(lm.model, lm.model.null, nsim=1000)$test["PBtest","p.value"]
    }else{
      boots.pvalues <- NA
    }
    # effect size
    sse.model <- sum(residuals(lm.model)^2)
    sse.nullmodel <- sum(residuals(lm.model.null)^2)
    IntpartialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
    
    stats.results <- cbind(int_var, dependentvar, anova.int.pvalue, IntpartialRsq, boots.pvalues)
    
    #Extract gam input data
    df <- gam.data
    #Create a prediction data frame, used to estimate (posterior) model slopes (varying covariate coefficients)
    theseVars <- attr(terms(lm.model),"term.labels")[1:3]
    varClasses <- unlist(lapply(theseVars, function(x) class(gam.data[[x]])))
    
    model.tmp <- lm.model
  }
  
  #prediction df for int_var P10th at np smooth_var increments
  pred.low <- data.frame(init = rep(0,np)) 
  for (v in c(1:length(theseVars))) {
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[v]
    if (thisVar == int_var) { 
      #pred.low[,int_var] <- (min(df[,int_var],na.rm = T))
      pred.low[,int_var] <- (quantile(df[,int_var], probs=0.1, na.rm = T))
    } else if (thisVar == smooth_var) {
      pred.low[,smooth_var] = seq(smoothmin,smoothmax, length.out = np)
    } else {
      tab.tmp<-table(df[,thisVar])
      levelact<-which.max(tab.tmp)
      switch (thisClass,
              "numeric" = {pred.low[,thisVar] = median(df[[thisVar]])},
              "factor" = {
                pred.low[,thisVar] = factor(levels(df[[thisVar]])[levelact], levels = levels(df[[thisVar]]))
              },
              "ordered" = {pred.low[,thisVar] = levels(df[[thisVar]])[levelact]}
      )}}
  
  #prediction df for int_var P90th at np smooth_var increments
  pred.high <- data.frame(init = rep(0,np)) 
  for (v in c(1:length(theseVars))) {
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[v]
    if (thisVar == int_var) { 
      #pred.high[,int_var] <- (max(df[,int_var],na.rm = T))
      pred.high[,int_var] <- (quantile(df[,int_var], probs=0.9, na.rm = T))
    } else if (thisVar == smooth_var) {
      pred.high[,smooth_var] = seq(smoothmin,smoothmax, length.out = np)
    } else {
      tab.tmp<-table(df[,thisVar])
      levelact<-which.max(tab.tmp)
      switch (thisClass,
              "numeric" = {pred.high[,thisVar] = median(df[[thisVar]])},
              "factor" = {
                pred.high[,thisVar] = factor(levels(df[[thisVar]])[levelact], levels = levels(df[[thisVar]]))
              },
              "ordered" = {pred.high[,thisVar] = levels(df[[thisVar]])[levelact]}
      )}}
  
  pred <- rbind(pred.low, pred.high) #complete pred df 
  pred <- pred %>% dplyr::select(-init)
  pred[[theseVars[which(varClasses=="factor")]]] <- factor(pred[[theseVars[which(varClasses=="factor")]]], levels = levels(gam.data[[theseVars[which(varClasses=="factor")]]]))
  #Get effects (slopes) along the smooth function for the true model
  if (return_posterior_coefficients == FALSE & bestmodel=="LM"){
    #varying coefficient slopes
    predicted.values <- predict(object = model.tmp, newdata = pred, type="response", re.form=NA) #predict y at min and mix int_var along the smooth function
    predicted.df <- pred
    predicted.df$fitted <- predicted.values
    predicted.values <- predicted.df %>% dplyr::select(fitted, all_of(smooth_var), all_of(int_var))
    colnames(predicted.values) <- c("fitted", "smooth.var", "int.var")
    predicted.values$smooth.var <- round(predicted.values$smooth.var, 3)
    
    varyingcoeff.slopes <-  predicted.values %>% #calculate the effect of int_var on y  (slope; delta y/delta int_var) along the smooth function 
      group_by(smooth.var) %>%
      do(slope = diff(.$fitted)/diff(.$int.var)) %>%
      unnest(cols = c(slope))
    colnames(varyingcoeff.slopes) <- c(smooth_var, "slope")  
  }else if (return_posterior_coefficients == FALSE & bestmodel=="GAMM"){
    #varying coefficient slopes
    predicted.values <- fitted_values(object = gamm.model$gam, data = pred) #predict y at min and mix int_var along the smooth function
    predicted.values <- predicted.values %>% dplyr::select(.fitted, all_of(smooth_var), all_of(int_var))
    colnames(predicted.values) <- c("fitted", "smooth.var", "int.var")
    predicted.values$smooth.var <- round(predicted.values$smooth.var, 3)
    
    varyingcoeff.slopes <-  predicted.values %>% #calculate the effect of int_var on y  (slope; delta y/delta int_var) along the smooth function 
      group_by(smooth.var) %>%
      do(slope = diff(.$fitted)/diff(.$int.var)) %>%
      unnest(cols = c(slope))
    colnames(varyingcoeff.slopes) <- c(smooth_var, "slope")  
  }
  
  
  # Posterior distribution of slopes
  if (return_posterior_coefficients == TRUE & bestmodel=="GAMM") {
    Vb <- vcov(gamm.model$gam, unconditional = UNCONDITIONAL)
    sims <- MASS::mvrnorm(npd, mu = coef(gamm.model$gam), Sigma = Vb)
    X0 <- predict(gamm.model$gam, newdata = pred, type = "lpmatrix")
    predicted.values.posterior <- X0 %*% t(sims)
    
    predicted.values.posterior <- as.data.frame(predicted.values.posterior)
    colnames(predicted.values.posterior) <- sprintf("draw%s", seq_len(npd))
    predicted.values.posterior <- cbind(pred[, smooth_var], pred[, int_var], predicted.values.posterior)
    colnames(predicted.values.posterior)[1:2] <- c("smooth.var", "int.var")
    predicted.values.posterior$smooth.var <- round(predicted.values.posterior$smooth.var, 3)
    
    varyingcoeff.slopes.CI <- predicted.values.posterior %>%
      pivot_longer(cols = starts_with("draw"), names_to = "draw", values_to = "posterior.fitted") %>%
      mutate(int.var = ifelse(int.var == quantile(df[[int_var]], 0.1), "low", "high")) %>%
      pivot_wider(names_from = int.var, values_from = posterior.fitted) %>%
      mutate(slope = (high - low) / (quantile(df[[int_var]], 0.9) - quantile(df[[int_var]], 0.1))) %>%
      select(smooth.var, draw, slope)
    varyingcoeff.slopes.CI <- cbind(as.character(dependentvar), varyingcoeff.slopes.CI)
    colnames(varyingcoeff.slopes.CI) <- c("label", smooth_var, "draw", "slope")
  }else if (return_posterior_coefficients == TRUE & bestmodel=="LM"){
    beta_hat <- lme4::fixef(model.tmp)
    Vb <- vcov(model.tmp)
    sims <- MASS::mvrnorm(npd, mu = beta_hat, Sigma = Vb)
    X0 <- model.matrix(delete.response(terms(model.tmp)), data = pred)
    predicted.values.posterior <- X0 %*% t(sims)
    
    predicted.values.posterior <- as.data.frame(predicted.values.posterior)
    colnames(predicted.values.posterior) <- sprintf("draw%s", seq_len(npd))
    predicted.values.posterior <- cbind(pred[, smooth_var], pred[, int_var], predicted.values.posterior)
    colnames(predicted.values.posterior)[1:2] <- c("smooth.var", "int.var")
    predicted.values.posterior$smooth.var <- round(predicted.values.posterior$smooth.var, 3)
    
    varyingcoeff.slopes.CI <- predicted.values.posterior %>%
      pivot_longer(cols = starts_with("draw"), names_to = "draw", values_to = "posterior.fitted") %>%
      mutate(int.var = ifelse(int.var == quantile(df[[int_var]], 0.1), "low", "high")) %>%
      pivot_wider(names_from = int.var, values_from = posterior.fitted) %>%
      mutate(slope = (high - low) / (quantile(df[[int_var]], 0.9) - quantile(df[[int_var]], 0.1))) %>%
      select(smooth.var, draw, slope)
    varyingcoeff.slopes.CI <- cbind(as.character(dependentvar), varyingcoeff.slopes.CI)
    colnames(varyingcoeff.slopes.CI) <- c("label", smooth_var, "draw", "slope")
  }
  
  
  if(return_posterior_coefficients == FALSE){
    data.results1 <- list(stats.results, varyingcoeff.slopes, pred)
    return(data.results1)
  }
  
  if(return_posterior_coefficients == TRUE){
    data.results2 <- list(stats.results,  varyingcoeff.slopes.CI, pred)
    return(data.results2)
  }
}
