## =========================
## 0) Packages
## =========================
suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(tidyr); library(forcats)
  library(ggplot2); library(mgcv); library(purrr); library(multcomp)
})

## =========================
## 1) I/O
## =========================
path  <- "D:/datasets/yunfu/interfile_folder/temp"
files <- c(back1 = "back1_deviation_PHQy09.rds",
           gngd  = "gngd_deviation_PHQy09.rds",
           back2 = "back2_deviation_PHQy09.rds")
raw_list <- lapply(file.path(path, files), readRDS)

## =========================
## 2) Helpers (robust cleaning)
## =========================
find_col <- function(df, patterns){
  cand <- names(df)[Reduce(`|`, lapply(patterns, function(p) grepl(p, names(df), ignore.case = TRUE)))]
  if (length(cand) == 0) NA_character_ else cand[1]
}
normalize_sex_factor <- function(x){
  z <- toupper(trimws(as.character(x)))
  z[z %in% c("FEMALE","F","女","WOMAN","WOMEN")] <- "F"
  z[z %in% c("MALE","M","男","MAN","MEN")]       <- "M"
  z[z %in% c("0","2")] <- "F"; z[z %in% c("1")] <- "M"
  z[!(z %in% c("F","M"))] <- NA
  factor(z, levels = c("F","M"))
}
extract_phq_level <- function(x){
  if (is.numeric(x)) lv <- suppressWarnings(as.integer(x)) else {
    s <- as.character(x); lv <- suppressWarnings(as.integer(stringr::str_extract(s, "[0-3]")))
  }
  factor(lv, levels = c(0,1,2,3))
}
clean_df <- function(df, task, file_label){
  dev_col <- find_col(df, c("deviationZ$", "deviation_z$", "dev.*z$","zscore$","deviation$"))
  age_col <- find_col(df, c("^Age_year$", "Age\\s*_?year"))
  phq_col <- find_col(df, c("^PHQ_y09$", "PHQ.*y09", "PHQ.?_?09"))
  sex_col <- find_col(df, c("^Gender$", "^Sex$", "性别"))
  if (is.na(dev_col) || is.na(age_col) || is.na(phq_col)) {
    stop(sprintf("【%s】缺列：deviationZ=%s, Age_year=%s, PHQ=%s",
                 file_label, dev_col, age_col, phq_col))
  }
  keep <- c(age_col, phq_col, dev_col, sex_col); keep <- keep[!is.na(keep)]
  df2  <- df[, keep, drop = FALSE]
  names(df2)[names(df2) == age_col] <- "Age_year"
  names(df2)[names(df2) == phq_col] <- "PHQ_raw"
  names(df2)[names(df2) == dev_col] <- "deviationZ"
  if (!is.na(sex_col)) names(df2)[names(df2) == sex_col] <- "Sex"
  df2 <- df2 %>%
    mutate(Task       = task,
           PHQ_y09    = extract_phq_level(PHQ_raw),
           Age_year   = suppressWarnings(as.numeric(Age_year)),
           deviationZ = suppressWarnings(as.numeric(deviationZ))) %>%
    filter(!is.na(PHQ_y09), !is.na(deviationZ), !is.na(Age_year))
  if ("Sex" %in% names(df2)) {
    sex1 <- normalize_sex_factor(df2$Sex); recog_rate <- mean(!is.na(sex1))
    if (recog_rate >= 0.6) { df2$Sex <- sex1; df2 <- df2 %>% filter(!is.na(Sex)) }
    else df2$Sex <- NULL
  }
  df2
}

df_all <- bind_rows(
  clean_df(raw_list[[1]], "Back1", files[1]),
  clean_df(raw_list[[2]], "GNGd",  files[2]),
  clean_df(raw_list[[3]], "Back2", files[3])
)

## 固定任务顺序：从上到下显示 GNGd → Back1 → Back2
df_all$Task <- factor(df_all$Task, levels = c("GNGd","Back1","Back2"))
tasks <- levels(df_all$Task)

## 可读的任务标签（面板标题）：Go/No-Go, 1-back, 2-back
task_labels <- c(GNGd = "Go/No-Go", Back1 = "1-back", Back2 = "2-back")
facet_labeller <- labeller(Task = function(x) task_labels[x])

## =========================
## 3) 拟合（GAM ）
## =========================
min_n_per_group  <- 2
min_n_combined   <- 4
prefer_ref_value <- "0"

try_fit_gam <- function(dft, task){
  if (nrow(dft) == 0) return(NULL)
  has_sex <- "Sex" %in% names(dft)
  
  # --- 数据准备与分组逻辑 ---
  tab <- dft %>% count(PHQ_y09, name = "n", .drop = FALSE) %>% arrange(desc(n))
  present <- as.character(tab$PHQ_y09[!is.na(tab$PHQ_y09)])
  if (length(present) < 2) return(NULL) # 至少需要两组才能比较
  
  # 确定参考组
  ref <- if (prefer_ref_value %in% present) prefer_ref_value else present[1]
  others <- setdiff(present, ref)
  
  # 检查其他组的样本量
  ok_others <- tab %>% 
    filter(PHQ_y09 %in% others, n >= min_n_per_group) %>% 
    pull(PHQ_y09) %>% 
    as.character()
  
  # 创建用于建模的专属数据框，避免覆盖原始 dft
  dft_model <- dft 
  
  if (length(ok_others) == 0) {
    # [改进] 移除硬编码：只要有其他组存在，就尝试合并
    if (length(others) > 0) {
      # 合并所有非参考组
      dft_model <- dft_model %>% 
        mutate(PHQ_y09 = fct_other(PHQ_y09, keep = ref, other_level = "other_combined"))
      
      # 检查合并后的组样本量是否足够
      n_combined <- dft_model %>% 
        filter(PHQ_y09 == "other_combined") %>% 
        nrow()
      
      if (n_combined < min_n_combined) return(NULL)
    } else {
      # 只有参考组，没有其他组
      return(NULL)
    }
  } else {
    # 保留样本量足够的组
    dft_model <- dft_model %>% 
      filter(PHQ_y09 %in% c(ref, ok_others)) %>% 
      droplevels()
  }
  
  # --- 模型拟合 ---
  # 检查年龄分布
  uniq_age <- n_distinct(dft_model$Age_year)
  if (uniq_age < 3) return(NULL)
  
  # 【修改 1】 k 值确定为 3
  k_use <- 3
  
  # 首选模型公式 (带交互)
  form1 <- if (has_sex)
    deviationZ ~ PHQ_y09 + Sex + s(Age_year, by = PHQ_y09, k = k_use, fx = F) + s(Age_year, k = k_use, fx = F)
  else 
    deviationZ ~ PHQ_y09 + s(Age_year, by = PHQ_y09, k = k_use)
  
  fit <- try(mgcv::gam(form1, data = dft_model, method = "REML"), silent = TRUE)
  
  # 【修改 2】不需要备用模型，如果 form1 无法拟合，输出警告信息即可
  if (inherits(fit, "try-error")) {
    warning(sprintf("【Task: %s】GAM model failed to converge for the specified formula. Returning NULL. Error: %s", task, as.character(fit)))
    return(NULL)
  }
  
  # 返回结果，levels 来自最终用于建模的数据
  list(fit = fit, ref = ref, levels = levels(dft_model$PHQ_y09))
}

meta_list <- map(tasks, ~ try_fit_gam(df_all %>% filter(Task == .x), .x)) %>% set_names(tasks)


## --- 统一生成 t & β 曲线：包含 1vs0,2vs0,3vs0,3vs2,2vs1 ---
tcurve_for_task <- function(task){
  dft <- df_all %>% filter(Task == task)
  meta <- meta_list[[task]]
  
  wanted_pairs <- list(c("1","0"), c("2","0"), c("3","0"), c("3","2"), c("2","1"), c("3", "1"))  
  # 【修改 3】如果没有 GAM 模型，直接输出报错信息，不需要回退方案
  if (is.null(meta)) {
    stop(sprintf("【Task: %s】No valid GAM model was fitted. Cannot generate t-curves.", task))
  }
  
  # --- 以下是原有的 GAM 逻辑 ---
  fit <- meta$fit; lev <- levels(fit$model$PHQ_y09)
  ages <- seq(min(fit$model$Age_year, na.rm=TRUE), max(fit$model$Age_year, na.rm=TRUE), length.out = 200)
  pairs_use <- Filter(function(p) all(p %in% lev), wanted_pairs)
  
  tcurve_contrast_marginal_sex <- function(fit, ages, a, b){
    has_sex <- "Sex" %in% names(fit$model)
    make_nd <- function(phq, sex_level = NULL) {
      nd <- data.frame(Age_year = ages, PHQ_y09 = factor(phq, levels = levels(fit$model$PHQ_y09)))
      resp <- names(fit$model)[1]
      others <- setdiff(names(fit$model), c(resp,"Age_year","PHQ_y09", if (has_sex) "Sex"))
      for (v in others) {
        x <- fit$model[[v]]
        if (is.factor(x)) nd[[v]] <- factor(levels(x)[1], levels = levels(x)) else nd[[v]] <- mean(x, na.rm=TRUE)
      }
      if (has_sex) {
        lv_sex <- levels(fit$model$Sex); if (is.null(sex_level)) sex_level <- lv_sex[1]
        nd$Sex <- factor(sex_level, levels = lv_sex)
      }
      nd
    }
    if (!has_sex) {
      Xa <- predict(fit, newdata = make_nd(a), type = "lpmatrix")
      Xb <- predict(fit, newdata = make_nd(b), type = "lpmatrix")
      D  <- Xa - Xb
    } else {
      wtab <- table(fit$model$Sex); wtab <- wtab/sum(wtab); lv_sex <- levels(fit$model$Sex); D <- 0
      for (sx in lv_sex) {
        Xa <- predict(fit, newdata = make_nd(a, sx), type = "lpmatrix")
        Xb <- predict(fit, newdata = make_nd(b, sx), type = "lpmatrix")
        D  <- D + as.numeric(wtab[sx]) * (Xa - Xb)
      }
    }
    beta_vec <- coef(fit); V <- vcov(fit)
    est <- as.numeric(D %*% beta_vec); se <- sqrt(rowSums((D %*% V) * D))
    tibble(Age = ages, tval = est/se, beta = est)
  }
  
  out <- map_dfr(pairs_use, function(p){
    a <- p[[1]]; b <- p[[2]]
    tc <- tcurve_contrast_marginal_sex(fit, ages, a, b)
    if (is.null(tc)) return(NULL)
    tc %>% mutate(Contrast = sprintf("%s: PHQ %s vs %s", task, a, b))
  })
  
  return(out)
}

tcurves <- map_dfr(tasks, tcurve_for_task)
if (nrow(tcurves) == 0) stop("清洗后各任务都缺少可比较的 PHQ 组或 Age 分布，无法生成曲线。")

## =========================
## 4) Tukey–Kramer（single-step）按年龄校正 + 热图
## =========================
# ... 后续代码无需修改，保持原样 ...
# ... The rest of the code for section 4 and plotting remains the same ...

ALPHA <- 0.05

# 基础绘图数据
heat_df <- tcurves %>%
  filter(is.finite(Age), is.finite(tval), is.finite(beta)) %>%
  mutate(
    Task      = sub(":.*$", "", Contrast),
    Contrast2 = sub("^[^:]+:\\s*", "", Contrast)  # "PHQ 1 vs 0" 等
  )

# 用于 multcomp 的辅助：在某个年龄生成某 PHQ 水平的（性别边际化）设计矩阵
.make_nd <- function(fit, ages, phq, sex_level = NULL){
  lev <- levels(fit$model$PHQ_y09)
  has_sex <- "Sex" %in% names(fit$model)
  nd <- data.frame(Age_year = ages, PHQ_y09 = factor(phq, levels = lev))
  resp <- names(fit$model)[1]
  others <- setdiff(names(fit$model), c(resp,"Age_year","PHQ_y09", if (has_sex) "Sex"))
  for (v in others) {
    x <- fit$model[[v]]
    if (is.factor(x)) nd[[v]] <- factor(levels(x)[1], levels = levels(x)) else nd[[v]] <- mean(x, na.rm = TRUE)
  }
  if (has_sex) {
    lv_sex <- levels(fit$model$Sex); if (is.null(sex_level)) sex_level <- lv_sex[1]
    nd$Sex <- factor(sex_level, levels = lv_sex)
  }
  nd
}

# 对单任务、单年龄做 Tukey–Kramer；包含 1vs0,2vs0,3vs0,3vs2,2vs1
tk_tukey_by_age <- function(task_name){
  meta <- meta_list[[task_name]]
  if (is.null(meta)) return(NULL)  # 该 task 没有可用的 GAM
  fit <- meta$fit
  lev <- levels(fit$model$PHQ_y09)
  if (length(lev) < 2) return(NULL)
  
  ages <- seq(min(fit$model$Age_year, na.rm = TRUE),
              max(fit$model$Age_year, na.rm = TRUE), length.out = 200)
  
  has_sex <- "Sex" %in% names(fit$model)
  w_sex <- if (has_sex) { p <- table(fit$model$Sex); p/sum(p) } else NULL
  lv_sex <- if (has_sex) levels(fit$model$Sex) else NULL
  
  out <- vector("list", length(ages))
  
  for (i in seq_along(ages)){
    a <- ages[i]
    
    X_by_phq <- lapply(lev, function(phq){
      if (!has_sex) {
        predict(fit, newdata = .make_nd(fit, a, phq), type = "lpmatrix")
      } else {
        Reduce(`+`, lapply(lv_sex, function(sx){
          as.numeric(w_sex[sx]) *
            predict(fit, newdata = .make_nd(fit, a, phq, sex_level = sx), type = "lpmatrix")
        }))
      }
    })
    names(X_by_phq) <- lev
    
    pairs <- combn(lev, 2, simplify = FALSE)
    L_list <- lapply(pairs, function(ab){
      a1 <- ab[1]; b1 <- ab[2]
      X_by_phq[[b1]] - X_by_phq[[a1]]
    })
    if (length(L_list) == 0) next
    L <- do.call(rbind, L_list)
    rownames(L) <- vapply(pairs, function(ab) paste0("PHQ ", ab[2], " vs ", ab[1]), "")
    
    gl <- try(multcomp::glht(model = fit, linfct = L), silent = TRUE)
    if (inherits(gl, "try-error")) next
    sm <- summary(gl, test = adjusted("single-step"))
    
    out[[i]] <- tibble::tibble(
      Age = a,
      Task = task_name,
      Contrast2 = rownames(L),
      p_adj_tukey = as.numeric(sm$test$pvalues)
    )
  }
  
  dplyr::bind_rows(out)
}


tk_all <- map_dfr(tasks, tk_tukey_by_age)

# 合并 Tukey–Kramer 校正 p 值
heat_df <- heat_df %>% left_join(tk_all, by = c("Task","Age","Contrast2")) %>%
  mutate(sig_tk = !is.na(p_adj_tukey) & p_adj_tukey < ALPHA)

# y 轴对比顺序
lev_contrast <- c("PHQ 3 vs 2", "PHQ 2 vs 1", "PHQ 3 vs 1", "PHQ 3 vs 0", "PHQ 2 vs 0", "PHQ 1 vs 0")
heat_df$Contrast2 <- factor(heat_df$Contrast2,
                            levels = intersect(lev_contrast, unique(heat_df$Contrast2)))

# 连续显著年龄段 → 黑框
rect_df <- heat_df %>%
  arrange(Task, Contrast2, Age) %>%
  group_by(Task, Contrast2) %>%
  mutate(run = cumsum(c(0, diff(sig_tk) != 0))) %>%
  group_by(Task, Contrast2, run) %>%
  filter(first(sig_tk)) %>%
  summarise(xmin = min(Age), xmax = max(Age), .groups = "drop") %>%
  mutate(
    y_id = as.numeric(factor(Contrast2, levels = levels(heat_df$Contrast2))),
    ymin = y_id - 0.49, ymax = y_id + 0.49
  )

heat_df$Task <- factor(heat_df$Task, levels = c("GNGd", "Back1", "Back2"))
rect_df$Task <- factor(rect_df$Task, levels = c("GNGd", "Back1", "Back2"))



## =========================
## 绘图
## =========================
p_beta_grad <- ggplot(heat_df, aes(x = Age, y = Contrast2, fill = beta)) +
  geom_raster(interpolate = FALSE) +
  geom_rect(data = rect_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = NA, color = "black", linewidth = 0.5) +
  facet_wrap(~ Task, ncol = 1, scales = "free_y", labeller = facet_labeller) +
  scale_fill_gradient2(
    low = "#2c7fb8", mid = "white", high = "#d73027",
    midpoint = 0,
    limits = c(-0.5, 0.1),
    oob = scales::squish,
    name = expression(beta)
  ) +
  labs(
    title = "Beta heatmap (Tukey–Kramer corrected)",
    x = "Age (years)", y = "Contrast"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5)
  )

print(p_beta_grad)