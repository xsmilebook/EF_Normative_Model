# correlation between ordinal variables & continuous variables.
# age effects need to be removed from the continous variables.
library(polycor)
library(tidyverse)
library(mgcv)
library(psych)

ordinalcorr <- function(dependentvar, dataname, interest.indep.var, covariates,smoothvar, knots, set_fx = FALSE, stats_only = FALSE){
  
  #Fit the gam
  gam.data <- get(dataname)
  Independent.var<-gam.data[[interest.indep.var]]
  
  NonNANIndex <- which(!is.na(gam.data[ ,interest.indep.var]) & !is.na(gam.data[ ,dependentvar]))
  gam.data <- gam.data[NonNANIndex,]
  
  parcel <- as.character(dependentvar)
  if (length(levels(gam.data[[interest.indep.var]]))>2){
    modelformula <- as.formula(sprintf("%s ~ as.factor(%s) + %s + s(%s, k=%d, fx=F)",dependentvar, interest.indep.var,  covariates, smoothvar, knots))
  }else{
    modelformula <- as.formula(sprintf("%s ~ %s + %s +  s(%s, k=%d, fx=F)",dependentvar, interest.indep.var,  covariates, smoothvar, knots))
  }
  modelformula.null<-as.formula(sprintf("%s ~  %s + s(%s, k=%d, fx=F)",dependentvar,  covariates, smoothvar, knots))
  gam.model <- gam(modelformula, method="REML", data = gam.data,family=ocat(R=4))
  gam.results <- summary(gam.model)
  gam.model.null <- gam(modelformula.null, method="REML", data = gam.data,family=ocat(R=4))

  
  #GAM statistics
  #F value for the smooth term and GAM-based significance of the smooth term
  gam.smooth.t <- gam.results$p.table[2,3]
  gam.smooth.pvalue <- gam.results$p.table[2,4]
  
  #Full versus reduced model anova p-value
  anova.cov.pvalue <- anova.gam(gam.model.null,gam.model,test='Chisq')$`Pr(>Chi)`[2]
  # if (stats_only){
  #   anova.cov.pvalue <- anovaPB(gam.model.null,gam.model, n.sim = 1000,test='Chisq')$`Pr(>Chi)`[2]
  # }else{anova.cov.pvalue <- NA}
  # if(is.na(anova.cov.pvalue)){ #if residual deviance is exactly equal between full and reduced models and p=value = NA, set p = 1
  #   anova.cov.pvalue <- 1}
  
  ##Full versus reduced model direction-dependent partial R squared
  ### effect size
  sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
  sse.nullmodel <- sum((gam.model.null$y - gam.model.null$fitted.values)^2)
  partialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  
  ### effect direction
  if(gam.smooth.t < 0){ #if the gam t-value for covariate of interest is less than 0, make the partialRsq negative
    partialRsq <- partialRsq*-1}
  
  #residual correlation
  varcorformula1 <- as.formula(sprintf("%s ~  %s +  s(%s, k=%d, fx=F)", dependentvar,  covariates, smoothvar, knots))
  #res1<-residuals(gam(varcorformula1, method="REML", data = gam.data,  family = nb()))
  res1 <- gam.data[[dependentvar]]
  res2<- gam.data[[interest.indep.var]]
  PCorr_Test <-polycor::polyserial(res2, res1, ML=T, std.err=T)
  correstimate <- PCorr_Test$rho[[1]]
  
  samplesize <- nrow(gam.data)
  corrmethod <- "polyserial"
  beta <- gam.results$p.table[2, 1]
  stats.results <- cbind(parcel, interest.indep.var, gam.smooth.t, gam.smooth.pvalue,anova.cov.pvalue, partialRsq, corrmethod, correstimate, samplesize, beta)
  data.results <- list()
  data.results[[1]] <- as.data.frame(stats.results)
  data.results[[2]] <- data.frame(Dependent.res=res1, Independent.var=gam.data[[interest.indep.var]])
  if(stats_only == TRUE)
    return(stats.results)
  if(stats_only == FALSE)
    return(data.results)
  
}



