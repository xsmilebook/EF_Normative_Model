library(tidyr)
library(mgcv)
library(psych)
library(dplyr)
library(ecostats)
library(pbkrtest)

#### Fit the linear effects of a continuous variable with a smooth item as covariates ####
## Function to fit evaluate the linear effects of a continuous variable on dependent variables per edge.
gam.fit.Independent.var <- function(dependentvar, dataname, interest.indep.var, covariates, stats_only = FALSE){  
  #Fit the gam
  gam.data <- get(dataname)
  Independent.var<-gam.data[[interest.indep.var]]

  outlierindx1<-which(Independent.var<mean(Independent.var,na.rm = T)-3*sd(Independent.var,na.rm = T) | Independent.var>mean(Independent.var,na.rm = T)+3*sd(Independent.var,na.rm = T))
  gam.data[outlierindx1 ,interest.indep.var]<-NA
  NonNANIndex <- which(!is.na(gam.data[ ,interest.indep.var]) & !is.na(gam.data[ ,dependentvar]))
  gam.data <- gam.data[NonNANIndex,]
  
  # tmp<-gam.data[,dependentvar]
  # outlierindx<-which(tmp<mean(tmp)-3*sd(tmp) | tmp>mean(tmp)+3*sd(tmp))
  # if (length(outlierindx)>0){
  #   gam.data<-gam.data[-outlierindx, ]
  # }
  parcel <- as.character(dependentvar)
  modelformula <- as.formula(sprintf("%s ~ %s + %s",dependentvar, interest.indep.var, covariates))
  modelformula.null<-as.formula(sprintf("%s ~ %s",dependentvar, covariates))
  gam.model <- gam(modelformula, method="REML", data = gam.data)
  gam.results <- summary(gam.model)
  gam.model.null <- gam(modelformula.null, method="REML", data = gam.data)
  
  #GAM statistics
  #F value for the smooth term and GAM-based significance of the smooth term
  gam.indep.t <- gam.results$p.table[2,3]
  gam.indep.pvalue <- gam.results$p.table[2,4]
  
  #Full versus reduced model anova p-value
  anova.cov.pvalue <- anova.gam(gam.model.null,gam.model,test='Chisq')$`Pr(>Chi)`[2]
  # Full versus reduced model anova p-value
  if (stats_only) {
    # Perform ANOVA simulation
    anova_results <- anovaPB(gam.model.null, gam.model, n.sim = 1000, test = 'Chisq')
    anova.pvalues <- anova_results$`Pr(>Chi)`[2]  # Get 1000 p-values
  }
  # if (stats_only){
  #   anova.cov.pvalue <- anovaPB(gam.model.null,gam.model, n.sim = 1000,test='Chisq')$`Pr(>Chi)`[2]
  # }else{anova.cov.pvalue <- NA}
  if(is.na(anova.cov.pvalue)){ #if residual deviance is exactly equal between full and reduced models and p=value = NA, set p = 1
    anova.cov.pvalue <- 1}
  
  ##Full versus reduced model direction-dependent partial R squared
  ### effect size
  sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
  sse.nullmodel <- sum((gam.model.null$y - gam.model.null$fitted.values)^2)
  partialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  ### effect direction
  if(gam.indep.t < 0){ #if the gam t-value for covariate of interest is less than 0, make the partialRsq negative
    partialRsq <- partialRsq*-1}
  
  #residual correlation
  varcorformula1 <- as.formula(sprintf("%s ~  %s", dependentvar, covariates))
  varcorformula2 <- as.formula(sprintf("%s ~  %s", interest.indep.var, covariates))
  res1<-residuals(gam(varcorformula1,method="REML", data = gam.data))
  res2<-residuals(gam(varcorformula2,method="REML", data = gam.data))
  res1.normtest <- ks.test(res1, "pnorm", mean = mean(res1), sd = sd(res1))
  res2.normtest <- ks.test(res2, "pnorm", mean = mean(res2), sd = sd(res2))
  
  if (res1.normtest$p.value > 0.01 & res2.normtest$p.value > 0.01){
    corrmethod <- "pearson"
  }else{
    corrmethod <- "spearman"
  }
  
  PCorr_Test <- corr.test(res1, res2, method=corrmethod)
  correstimate<-as.numeric(PCorr_Test$r)
  corrp <- as.numeric(PCorr_Test$p)
  
  samplesize <- nrow(gam.data)
  stats.results <- cbind(parcel, interest.indep.var, gam.indep.t, gam.indep.pvalue,anova.cov.pvalue, anova.pvalues, partialRsq, corrmethod, correstimate, corrp, samplesize)
  data.results <- list()
  data.results[[1]] <- as.data.frame(stats.results)
  data.results[[2]] <- data.frame(SCres=res1, cogres=res2)
  if(stats_only == TRUE)
    return(stats.results)
  if(stats_only == FALSE)
    return(data.results)
}




