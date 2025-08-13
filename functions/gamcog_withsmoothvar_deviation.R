library(tidyr)
library(mgcv)
library(psych)
library(dplyr)
library(ecostats)
library(pbkrtest)

# modified anovaPB: return distribution of simulated statistics
anovaPB_ext <- function(objectNull, object, n.sim = 999, 
                        colRef = switch(class(object)[1], "lm" = 5, "lmerMod" = 6, "glmmTMB" = 6, 4),
                        rowRef = 2, ncpus = NULL, ...) {
  
  # check model
  if (length(unlist(coef(objectNull))) > length(unlist(coef(object))))
    stop("The first object (null model) should be smaller than the alternative model.")
  
  # 获取响应变量的维度名称
  respDimnames <- dimnames(model.response(model.frame(object)))
  
  # 设置 CPU 数量
  if (is.null(ncpus)) {
    ncpus <- min(70, parallel::detectCores() - 2)  
  }
  cat(paste("This run will allocate ", ncpus, " cores."))
  # 计算原始 ANOVA 表（真实数据）
  targs <- match.call(expand.dots = FALSE)
  anovaFn <- anova
  statObs <- try(anova(objectNull, object, ...))
  
  if (inherits(statObs, "try-error")) {
    anovaFn <- function(objectNull, object, ...) {
      llAlt  <- logLik(object)
      llNull <- logLik(objectNull)
      table <- data.frame(
        df = c(attr(llNull, "df"), attr(llAlt, "df")),
        deviance = -2 * c(llNull, llAlt),
        LRT = c(NA, -2 * llNull + 2 * llAlt)
      )
      return(table)
    }
    statObs <- anovaFn(objectNull, object)
    statObs$P <- c(NA, NA)
    names(statObs)[4] <- 'Pr(>LRT)'
    colRef <- 3
    modelnamelist <- c(deparse(substitute(objectNull)), deparse(substitute(object)))
    Xnames <- c(paste(deparse(formula(objectNull), width.cutoff = 500), collapse = "\n"),
                paste(deparse(formula(object), width.cutoff = 500), collapse = "\n"))
    topnote <- paste(modelnamelist, ": ", Xnames, sep = "", collapse = "\n")
    title <- "Analysis of Deviance Table\n"
    rownames(statObs) <- modelnamelist
    attr(statObs, "heading") <- c(title, topnote)
  }
  
  # 初始化统计量向量（第1个是真实值）
  stats <- rep(NA, n.sim + 1)
  stats[1] <- statObs[rowRef, colRef]
  
  # 提取模型框架
  if (inherits(object, c("lmerMod", "glmerMod"))) {
    cll <- object@call
    mf <- match.call(call = cll)
    dat <- if (.hasSlot(object, "data")) object@data else NULL
  } else {
    cll <- object$call
    mf <- match.call(call = cll)
    dat <- object$data
  }
  m <- match(c("formula", "data", "subset", "weights", "na.action", "etastart", "mustart", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)
  
  # 构建 model.frame
  modelF <- try(eval(mf, parent.frame()), silent = TRUE)
  if (inherits(modelF, "try-error") || inherits(object, c("lmerMod", "glmmTMB")))
    modelF <- model.frame(object)
  
  respName <- names(model.frame(object))[1]
  whichResp <- 1
  if (!is.null(dat) && is.list(dat)) {
    whichAdd <- which(names(dat) %in% names(modelF) == FALSE)
    for (iAdd in whichAdd) {
      if (!is.list(dat[[iAdd]])) modelF[[names(dat)[iAdd]]] <- dat[[iAdd]]
    }
  }
  
  offs <- try(model.offset(modelF), silent = TRUE)
  if (inherits(offs, "try-error")) offs <- NULL
  
  if (regexpr("(", respName, fixed = TRUE) > 0) {
    newResp <- sprintf("`%s`", respName)
    fm.update <- reformulate(".", response = newResp)
  } else {
    fm.update <- reformulate(".")
  }
  
  is.mva <- ncol(as.matrix(modelF[[whichResp]])) > 1
  
  # simulate response variable
  yNew <- simulate(objectNull, nsim = n.sim)
  
  # define simulate function
  getStat <- function(iSim, yNew, objectNull, object, modelF, anovaFn, is.mva, fm.update, whichResp, respDimnames, rowRef, colRef) {
    modelF[[whichResp]] <- if (is.mva) yNew[,,iSim] else as.matrix(yNew[, iSim], dimnames = respDimnames)
    
    if (is.null(offs)) {
      objectiNull <- update(objectNull, formula = fm.update, data = modelF)
      objecti     <- update(object,     formula = fm.update, data = modelF)
    } else {
      objectiNull <- update(objectNull, formula = fm.update, data = modelF, offset = offs)
      objecti     <- update(object,     formula = fm.update, data = modelF, offset = offs)
    }
    
    return(anovaFn(objectiNull, objecti, ...)[rowRef, colRef])
  }
  
  # simulate
  if (ncpus > 1) {
    # 🔥 关键：使用makeCluster而不是mclapply
    cl <- NULL
    tryCatch({
      cl <- parallel::makeCluster(ncpus)
      # 确保在函数结束时关闭集群
      on.exit(parallel::stopCluster(cl), add = TRUE)
      
      parallel::clusterExport(cl, c("getStat", "anovaFn", "yNew", "objectNull", "object", "modelF", 
                                    "is.mva", "fm.update", "whichResp", "respDimnames", 
                                    "rowRef", "colRef", as.character(cll[[1]])), envir = environment())
      
      statList <- parallel::clusterApplyLB(cl, 1:n.sim, getStat, yNew = yNew, objectNull = objectNull,
                                           object = object, modelF = modelF, anovaFn = anovaFn,
                                           is.mva = is.mva, fm.update = fm.update,
                                           whichResp = whichResp, respDimnames = respDimnames,
                                           rowRef = rowRef, colRef = colRef)
      stats[-1] <- unlist(statList)
    }, error = function(e) {
      # 如果并行失败，回退到串行
      warning("Parallel processing failed, falling back to serial execution: ", e$message)
      for (iBoot in 2:(n.sim + 1)) {
        stats[iBoot] <- getStat(iBoot - 1, yNew = yNew, objectNull = objectNull, object = object,
                                modelF = modelF, anovaFn = anovaFn, is.mva = is.mva,
                                fm.update = fm.update, whichResp = whichResp,
                                respDimnames = respDimnames, rowRef = rowRef, colRef = colRef)
      }
    })
  } else {
    # 串行执行
    for (iBoot in 2:(n.sim + 1)) {
      stats[iBoot] <- getStat(iBoot - 1, yNew = yNew, objectNull = objectNull, object = object,
                              modelF = modelF, anovaFn = anovaFn, is.mva = is.mva,
                              fm.update = fm.update, whichResp = whichResp,
                              respDimnames = respDimnames, rowRef = rowRef, colRef = colRef)
    }
  }
  
  # return anova table
  statReturn <- statObs[, 1:colRef, drop = FALSE]
  statReturn$`P-value` <- NA
  observed_stat <- stats[1]
  simulated_stats <- stats[-1]
  p_value <- mean(simulated_stats >= observed_stat - 1e-8, na.rm = TRUE)  # 单侧检验
  statReturn$`P-value`[rowRef] <- p_value
  
  hasP <- grep("Pr", colnames(statObs))
  if (length(hasP) > 0) colnames(statReturn)[ncol(statReturn)] <- colnames(statObs)[hasP[1]]
  

  attr(statReturn, "heading") <- attr(statObs, "heading")
  class(statReturn) <- c("anovaPB", class(statObs))
  
  attr(statReturn, "simulated_stats") <- simulated_stats
  attr(statReturn, "observed_stat")  <- observed_stat
  attr(statReturn, "n.sim")          <- n.sim
  
  attr(statReturn, "full_stats_vector") <- stats  
  
  return(statReturn)
}

#### Fit the linear effects of a continuous variable with a smooth item as covariates ####
## Function to fit evaluate the linear effects of a continuous variable on dependent variables per edge.
gam.fit.Independent.var <- function(dependentvar, dataname, smoothvar, interest.indep.var, covariates, stats_only = FALSE){  
  #Fit the gam
  gam.data <- get(dataname)
  gam.data <- gam.data %>%
    filter(!is.na(.data[[interest.indep.var]]), !is.na(.data[[dependentvar]]))
  
  Independent.var<-gam.data[[interest.indep.var]]

  indep_outlier_index <- which(Independent.var < mean(Independent.var) - 3 * sd(Independent.var) |
                                 Independent.var > mean(Independent.var) + 3 * sd(Independent.var))
  if (length(indep_outlier_index) > 0) {
    gam.data <- gam.data[-indep_outlier_index, ]
  }

  dependent.var.vec <- gam.data[[dependentvar]]
  dep_outlier_index <- which(dependent.var.vec < mean(dependent.var.vec) - 3 * sd(dependent.var.vec) |
                               dependent.var.vec > mean(dependent.var.vec) + 3 * sd(dependent.var.vec))
  if (length(dep_outlier_index) > 0) {
    gam.data <- gam.data[-dep_outlier_index, ]
  }
  
  parcel <- as.character(dependentvar)
  modelformula <- as.formula(sprintf("%s ~ %s + %s + s(%s, k=3, fx=F)",dependentvar, interest.indep.var, covariates, smoothvar))
  modelformula.null<-as.formula(sprintf("%s ~ %s + s(%s, k=3, fx=F)",dependentvar, covariates, smoothvar))
  gam.model <- gam(modelformula, method="REML", data = gam.data)
  gam.results <- summary(gam.model)
  gam.model.null <- gam(modelformula.null, method="REML", data = gam.data)
  
  #gam statistics
  #F value for the smooth term and gam-based significance of the smooth term
  gam.indep.t <- gam.results$p.table[2,3]
  gam.indep.pvalue <- gam.results$p.table[2,4]
  
  # Full versus reduced model anova p-value
  if (stats_only){
    # Perform ANOVA simulation
    anova_results <- anovaPB_ext(gam.model.null, gam.model, n.sim = 10000, test = 'Chisq', ncpus = 660)
    anova.pvalues <- anova_results$`Pr(>Chi)`[2]  # Get 1000 p-values
    simulated_stats <- attr(anova_results, "simulated_stats")
    observed_stat <- attr(anova_results, "observed_stat")
  }else{
    anova.pvalues <- NA
    simulated_stats <- NULL
    observed_stat <- NULL
  }
  
  ##Full versus reduced model direction-dependent partial R squared
  ### effect size
  sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
  sse.nullmodel <- sum((gam.model.null$y - gam.model.null$fitted.values)^2)
  partialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  ### effect direction
  if(gam.indep.t < 0){ 
    #if the gam t-value for covariate of interest is less than 0, make the partialRsq negative
    partialRsq <- partialRsq*-1
  }
  
  #residual correlation
  varcorformula1 <- as.formula(sprintf("%s ~  %s+s(%s, k=3, fx=F)", dependentvar, covariates, smoothvar))
  res1<-residuals(gam(varcorformula1, method="REML", data = gam.data))
  res2<-gam.data[[interest.indep.var]]
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
  beta = gam.results$p.table[2,1]
  stats.results <- cbind(parcel, interest.indep.var, gam.indep.t, gam.indep.pvalue, anova.pvalues, partialRsq, corrmethod, correstimate, corrp, samplesize, beta)
  sim_results <- list(
    simulated_stats = simulated_stats,
    observed_stat = observed_stat,
    n_sim = if (!is.null(simulated_stats)) length(simulated_stats) else NA
  )
  data.results <- list()
  data.results[[1]] <- as.data.frame(stats.results)
  data.results[[2]] <- data.frame(SCres=res1, cogres=res2)
  data.results[[3]] <- sim_results
  if(stats_only == TRUE)
    return(list(stats = as.data.frame(stats.results), 
                simulation = sim_results))
  if(stats_only == FALSE)
    return(data.results)
}




