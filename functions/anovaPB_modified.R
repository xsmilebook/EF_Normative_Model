library(future)
library(furrr)
library(progressr)

anovaPB_modified_with_stats <- function (objectNull, object, n.sim = 999, 
                                         colRef = switch(class(object)[1], 
                                                         lm = 5, lmerMod = 6, 4), 
                                         rowRef = 2, ncpus = NULL, ...) {
  if (length(unlist(coef(objectNull))) > length(unlist(coef(object)))) 
    stop(paste("Whoah, the first object fitted a larger model than the second object...", 
               "\n", " it should be a smaller 'null' model fit, nested in the second object."))
  respDimnames = dimnames(model.response(model.frame(object)))
  if (is.null(ncpus)) {
    ncpus = max(1, parallel::detectCores() - 2)
    parl = ifelse(ncpus == 1, "no", "snow")
  }
  targs <- match.call(expand.dots = FALSE)
  anovaFn = anova
  statObs = try(anova(objectNull, object, ...))
  if (inherits(statObs, "try-error")) {
    anovaFn = function(objectNull, object, ...) {
      llAlt = logLik(object)
      llNull = logLik(objectNull)
      table = data.frame(df = c(attr(llNull, "df"), attr(llAlt, 
                                                         "df")), deviance = -2 * c(llNull, llAlt), LRT = c(NA, 
                                                                                                           -2 * llNull + 2 * llAlt))
      return(table)
    }
    statObs = anovaFn(objectNull, object)
    statObs$P = c(NA, NA)
    names(statObs)[4] = "Pr(>LRT)"
    colRef = 3
    modelnamelist = c(deparse(substitute(objectNull)), deparse(substitute(object)))
    Xnames <- list(paste(deparse(formula(objectNull), width.cutoff = 500), 
                         collapse = "\n"), paste(deparse(formula(object), 
                                                         width.cutoff = 500), collapse = "\n"))
    topnote <- paste(modelnamelist, ": ", Xnames, sep = "", 
                     collapse = "\n")
    title <- "Analysis of Deviance Table\n"
    rownames(statObs) = modelnamelist
    attr(statObs, "heading") <- c(title, topnote)
  }
  stats = rep(NA, n.sim + 1)
  stats[1] = statObs[rowRef, colRef]
  if (inherits(object, c("lmerMod", "glmerMod"))) {
    cll <- object@call
    mf <- match.call(call = cll)
    if (.hasSlot(object, "data")) 
      dat <- object@data
    else dat <- NULL
  }
  else {
    cll <- object$call
    mf <- match.call(call = cll)
    dat <- object$data
  }
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
               "etastart", "mustart", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)
  modelF <- try(eval(mf, parent.frame()), silent = TRUE)
  if (inherits(modelF, "try-error") | inherits(object, c("lmerMod", 
                                                         "glmerMod", "glmmTMB"))) 
    modelF <- model.frame(object)
  respName <- names(model.frame(object))[1]
  whichResp <- 1
  if (is.null(dat) == FALSE) 
    if (is.list(dat)) {
      whichAdd = which(names(dat) %in% names(modelF) == 
                         FALSE)
      if (length(whichAdd) > 0) 
        for (iAdd in whichAdd) if (is.list(dat[[iAdd]]) == 
                                   FALSE) 
          modelF[[names(dat)[iAdd]]] = dat[[iAdd]]
    }
  offs = NULL
  modelF$offs = try(model.offset(modelF))
  if (regexpr("(", respName, fixed = TRUE) > 0) {
    newResp = sprintf("`%s`", respName)
    fm.update = reformulate(".", response = newResp)
  }
  else fm.update = reformulate(".")
  is.mva = ncol(as.matrix(modelF[[whichResp]])) > 1
  yNew = simulate(objectNull, n.sim)
  
  getStat = function(iSim, yNew, objectNull = objectNull, object = object, 
                     modelF = modelF, anovaFn = anovaFn, is.mva = is.mva, 
                     fm.update = fm.update, whichResp = whichResp, respDimnames = respDimnames, 
                     rowRef = rowRef, colRef = colRef) {
    modelF[[whichResp]] = if (is.mva) 
      yNew[, , iSim]
    else as.matrix(yNew[, iSim], dimnames = respDimnames)
    if (inherits(modelF$offs, "try-error") | is.null(modelF$offs)) {
      objectiNull = update(objectNull, formula = fm.update, 
                           data = modelF)
      objecti = update(object, formula = fm.update, data = modelF)
    }
    else {
      objectiNull = update(object, formula = fm.update, 
                           data = modelF, offset = offs)
      objecti = update(object, formula = fm.update, data = modelF, 
                       offset = offs)
    }
    return(anovaFn(objectiNull, objecti, ...)[rowRef, colRef])
  }
  
  # compute cross node(future)
  if (ncpus > 1) {
    statList <- future_map(1:n.sim, function(i) {
      getStat(i, yNew = yNew, objectNull = objectNull, object = object, 
              modelF = modelF, anovaFn = anovaFn, is.mva = is.mva, 
              fm.update = fm.update, whichResp = whichResp, respDimnames = respDimnames, 
              rowRef = rowRef, colRef = colRef)
    }, .progress = FALSE)
    
    stats[2:(n.sim + 1)] = unlist(statList)
  }
  else {
    for (iBoot in 1:n.sim) {
      stats[iBoot + 1] = getStat(iBoot, yNew = yNew, objectNull = objectNull, object = object, 
                                 modelF = modelF, anovaFn = anovaFn, is.mva = is.mva, 
                                 fm.update = fm.update, whichResp = whichResp, respDimnames = respDimnames, 
                                 rowRef = rowRef, colRef = colRef)
    }
  }
  
  statReturn = statObs[, 1:colRef]
  statReturn$"P-value" = NA
  statReturn$"P-value"[rowRef] = mean(stats > (stats[1] - 1e-08), 
                                      na.rm = TRUE)
  attr(statReturn, "n.sim") = sum(is.na(stats) == FALSE) - 
    1
  hasP = grep("Pr", colnames(statObs))
  colnames(statReturn)[colRef + 1] = colnames(statObs)[hasP[1]]
  attr(statReturn, "heading") = attr(statObs, "heading")
  
  # return statistics
  statReturn$bootstrap_stats <- stats
  statReturn$observed_stat <- stats[1]
  statReturn$n_sim <- n.sim
  
  class(statReturn) = c("anovaPB", class(statObs))
  return(statReturn)
}

# 跨节点并行计算函数（使用 future）
distributed_anovaPB <- function(objectNull, object, n.sim = 10000, 
                                n_nodes = 10, ncpus_per_node = NULL, ...) {
  
  # 计算每个节点的任务量
  sim_per_node <- n.sim %/% n_nodes
  remainder <- n.sim %% n_nodes
  
  cat("Distributing", n.sim, "simulations across", n_nodes, "nodes\n")
  cat("Each node will run approximately", sim_per_node, "simulations\n")
  
  # 创建节点任务列表
  node_tasks <- 1:n_nodes %>%
    map(~ list(
      node_id = .x,
      n_sim = ifelse(.x <= remainder, sim_per_node + 1, sim_per_node)
    ))
  
  # 全局变量（需要在 worker 中可用）
  objectNull_global <<- objectNull
  object_global <<- object
  
  # 设置并行计划（可以根据需要调整）
  if (is.null(ncpus_per_node)) {
    ncpus_per_node <- max(1, parallel::detectCores() - 1)
  }
  
  # 使用多会话并行（模拟跨节点）
  plan(multisession, workers = min(n_nodes, parallel::detectCores()))
  
  # 跨节点并行执行
  with_progress({
    p <- progressor(steps = n_nodes)
    
    node_results <- node_tasks %>%
      future_map(~ {
        # 更新进度
        p(sprintf("Node %d/%d completed", .x$node_id, n_nodes))
        
        # 在每个"节点"上执行计算
        result <- anovaPB_modified_with_stats(
          objectNull = objectNull_global,
          object = object_global,
          n.sim = .x$n_sim,
          ncpus = ncpus_per_node,
          ...
        )
        
        # 返回精简的结果用于合并
        list(
          node_id = .x$node_id,
          n_sim = .x$n_sim,
          bootstrap_stats = result$bootstrap_stats,
          observed_stat = result$observed_stat,
          p_value = result$"P-value"[rowRef]  # 假设 rowRef 可访问
        )
      }, .progress = FALSE)
  })
  
  # 恢复原始计划
  plan(sequential)
  
  return(node_results)
}

# 合并分布式结果
merge_distributed_results <- function(node_results) {
  if(length(node_results) == 0) return(NULL)
  
  # 合并所有统计量
  all_bootstrap_stats <- c()
  observed_stats <- c()
  total_simulations <- 0
  
  # 第一个节点的结果作为模板
  template_result <- node_results[[1]]
  
  for(i in 1:length(node_results)) {
    result <- node_results[[i]]
    # 排除第一个观测值，只取模拟值
    all_bootstrap_stats <- c(all_bootstrap_stats, result$bootstrap_stats[-1])
    observed_stats <- c(observed_stats, result$observed_stat)
    total_simulations <- total_simulations + result$n_sim
  }
  
  # 观测值应该都相同，取第一个
  final_observed <- observed_stats[1]
  
  # 重构完整的统计量向量：观测值 + 所有模拟值
  merged_stats <- c(final_observed, all_bootstrap_stats)
  
  # 重新计算 p 值
  p_value <- mean(merged_stats[-1] > (merged_stats[1] - 1e-08), na.rm = TRUE)
  
  cat("Merged results from", length(node_results), "nodes\n")
  cat("Total simulations:", total_simulations, "\n")
  cat("Final p-value:", format(p_value, digits = 6), "\n")
  
  # 构建最终结果
  final_result <- list(
    bootstrap_stats = merged_stats,
    observed_stat = final_observed,
    total_simulations = total_simulations,
    p_value = p_value,
    node_results = node_results
  )
  
  return(final_result)
}

# 完整的分布式 ANOVA PB 工作流程
complete_distributed_anovaPB <- function(objectNull, object, n.sim = 10000, 
                                         n_nodes = 10, ncpus_per_node = NULL, ...) {
  
  cat("Starting distributed ANOVA PB with", n.sim, "simulations\n")
  cat("Using", n_nodes, "nodes\n")
  
  # 执行分布式计算
  node_results <- distributed_anovaPB(
    objectNull = objectNull,
    object = object,
    n.sim = n.sim,
    n_nodes = n_nodes,
    ncpus_per_node = ncpus_per_node,
    ...
  )
  
  # 合并结果
  final_result <- merge_distributed_results(node_results)
  
  return(final_result)
}

# 增强的可视化函数
enhanced_bootstrap_plot <- function(distributed_result, bins = 100) {
  if(is.null(distributed_result$bootstrap_stats)) {
    stop("No bootstrap statistics found in result")
  }
  
  stats <- distributed_result$bootstrap_stats
  observed <- distributed_result$observed_stat
  p_value <- distributed_result$p_value
  total_n <- distributed_result$total_simulations
  
  # 创建数据框
  df <- data.frame(statistic = stats[-1])  # 排除观测值
  
  # 绘图
  library(ggplot2)
  p <- ggplot(df, aes(x = statistic)) +
    geom_histogram(bins = bins, fill = "steelblue", alpha = 0.7, color = "white") +
    geom_density(aes(y = ..density.. * total_n * (max(statistic) - min(statistic)) / bins), 
                 color = "darkred", size = 1, alpha = 0.8) +
    geom_vline(xintercept = observed, color = "red", linetype = "dashed", size = 1.5) +
    labs(
      title = "Bootstrap Distribution of Test Statistic",
      subtitle = paste("Observed:", format(observed, digits = 4), 
                       " | Simulations:", format(total_n, big.mark = ","),
                       " | P-value:", format(p_value, digits = 3)),
      x = "Test Statistic",
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12)
    )
  
  # 添加观测值标签
  p <- p + annotate("text", x = observed, y = Inf, 
                    label = paste("Observed =", format(round(observed, 2), nsmall = 2)), 
                    vjust = 2, hjust = 0, color = "red", fontface = "bold", size = 4)
  
  return(p)
}