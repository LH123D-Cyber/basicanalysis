#' 差异分析方法建议
#'
#' 检验目标变量在不同分组间是否满足正态性和方差齐性，
#' 并据此推荐合适的差异检验与事后检验方法。
#'
#' @param data 数据框。
#' @param Objective 目标变量列名（字符型，只能选择一个）。
#' @param group 分组列名（字符型，只能选择一个）。
#' @param alpha 显著性水平，默认 0.05。
#'
#' @return 一个 list，包含 `Normality`、`Homogeneity`、`Decision` 等元素。
#' @export
DDT <- function(data,
                Objective,
                group,
                alpha = 0.05) {

  if (length(Objective) != 1) stop("Objective 只能选择一个变量")
  if (length(group) != 1) stop("group 只能指定一个分组列")
  if (!Objective %in% colnames(data))
    stop("Objective 变量不存在于数据中: ", Objective)
  if (!group %in% colnames(data))
    stop("group 列不存在于数据中: ", group)

  x <- data[[Objective]]
  g <- data[[group]]

  if (!is.numeric(x)) stop("Objective 必须是数值型变量")
  if (!is.factor(g)) g <- factor(g)

  # 去除缺失值
  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]; g <- droplevels(g[ok])

  groups <- levels(g)
  grp_data <- split(x, g)
  k <- length(groups)

  if (k < 2) stop("分组列至少需要 2 个水平")

  result <- list(
    variable = Objective,
    group    = group,
    alpha    = alpha,
    n_groups = k,
    groups   = groups,
    n_each   = sapply(grp_data, length)
  )

  # ============================================================
  # 1. 各组正态性（Shapiro-Wilk）
  # ============================================================
  norm_test <- data.frame(
    group  = groups,
    n      = as.integer(sapply(grp_data, length)),
    W      = NA_real_,
    p      = NA_real_,
    normal = NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(groups)) {
    xi <- grp_data[[i]]
    if (length(xi) >= 3 && length(xi) <= 5000) {
      sw <- tryCatch(shapiro.test(xi), error = function(e) NULL)
      if (!is.null(sw)) {
        norm_test$W[i] <- round(as.numeric(sw$statistic), 4)
        norm_test$p[i] <- round(sw$p.value, 4)
        norm_test$normal[i] <- ifelse(sw$p.value > alpha, "Yes", "No")
      } else {
        norm_test$normal[i] <- "Error"
      }
    } else {
      norm_test$normal[i] <- "n<3"
    }
  }
  result$Normality <- norm_test
  all_normal <- all(norm_test$normal == "Yes")

  # ============================================================
  # 2. 方差齐性检验
  # ============================================================
  # 2a. Levene 检验（基于中位数，稳健）
  if (requireNamespace("car", quietly = TRUE)) {
    lev <- tryCatch(
      car::leveneTest(x ~ g, data = data.frame(x = x, g = g)),
      error = function(e) NULL
    )
    if (!is.null(lev)) {
      lev_W <- lev$`F value`[1]
      lev_p <- lev$`Pr(>F)`[1]
    } else {
      lev_W <- NA_real_; lev_p <- NA_real_
    }
  } else {
    # 若未安装 car 包，手动实现基于中位数的 Levene 检验
    med <- tapply(x, g, median, na.rm = TRUE)
    z <- abs(x - med[as.character(g)])
    lev_fit <- aov(z ~ g)
    lev_tab <- summary(lev_fit)[[1]]
    lev_W <- lev_tab$`F value`[1]
    lev_p <- lev_tab$`Pr(>F)`[1]
  }

  # 2b. Bartlett 检验（对正态性敏感）
  brt <- tryCatch(
    bartlett.test(x ~ g, data = data.frame(x = x, g = g)),
    error = function(e) NULL
  )
  brt_K2 <- if (!is.null(brt)) as.numeric(brt$statistic) else NA_real_
  brt_p <- if (!is.null(brt)) brt$p.value else NA_real_

  result$Homogeneity <- list(
    Levene = list(W = round(lev_W, 4),
                  p = round(lev_p, 4),
                  homogeneous = ifelse(!is.na(lev_p) && lev_p > alpha,
                                       "Yes", "No")),
    Bartlett = list(K2 = round(brt_K2, 4),
                    p = round(brt_p, 4),
                    homogeneous = ifelse(!is.na(brt_p) && brt_p > alpha,
                                         "Yes", "No"))
  )
  var_homogeneous <- !is.na(lev_p) && lev_p > alpha

  # ============================================================
  # 3. 决策建议
  # ============================================================
  if (k == 2) {
    # 两组：无需事后检验
    if (all_normal && var_homogeneous) {
      method <- "独立样本 t 检验"
      posthoc <- "无（两组无需事后检验）"
      reason <- "两组数据均满足正态分布，且方差齐性"
    } else if (all_normal && !var_homogeneous) {
      method <- "Welch t 检验"
      posthoc <- "无（两组无需事后检验）"
      reason <- "两组数据均满足正态分布，但方差不齐"
    } else {
      method <- "Wilcoxon 秩和检验"
      posthoc <- "无（两组无需事后检验）"
      reason <- "数据不满足正态分布（Wilcoxon 对分布和方差不敏感）"
    }
  } else {
    # 三组及以上
    if (all_normal && var_homogeneous) {
      method <- "One-way ANOVA（单因素方差分析）"
      posthoc <- "Tukey HSD"
      reason <- "各组数据满足正态分布，且方差齐性"
    } else if (all_normal && !var_homogeneous) {
      method <- "Welch ANOVA（Welch 单因素方差分析）"
      posthoc <- "Games-Howell"
      reason <- "各组数据满足正态分布，但方差不齐"
    } else if (!all_normal && var_homogeneous) {
      method <- "Kruskal-Wallis 秩和检验"
      posthoc <- "Dunn 检验（BH 校正）"
      reason <- "数据不满足正态分布，但方差齐性"
    } else {
      method <- "Kruskal-Wallis 秩和检验"
      posthoc <- "Dunn 检验（BH 校正）"
      reason <- "数据不满足正态分布，且方差不齐（Kruskal-Wallis 对两者均不敏感）"
    }
  }

  # 小样本提示
  if (min(norm_test$n) < 5) {
    reason <- paste0(reason, "；部分组样本量小于 5，正态性检验结果可能不稳定")
  }

  summary_text <- paste0(
    Objective, "：", reason,
    "。建议使用 ", method, " 进行差异分析",
    if (k >= 3) paste0("，并使用 ", posthoc, " 进行事后检验。") else "。"
  )

  result$Decision <- list(
    normality_met   = all_normal,
    homogeneity_met = var_homogeneous,
    reason          = reason,
    test            = method,
    posthoc         = posthoc,
    summary         = summary_text
  )

  # ============================================================
  # 4. 打印输出
  # ============================================================
  cat("══════════ DDT 分析结果 ══════════\n")
  cat("变量:", Objective, " | 分组:", group,
      " | 组数:", k, " | alpha =", alpha, "\n\n")

  cat("【1. 各组正态性 (Shapiro-Wilk)】\n")
  print(norm_test, row.names = FALSE)
  cat("\n")

  cat("【2. 方差齐性】\n")
  cat(sprintf("Levene 检验 (中位数中心化): W = %s, p = %s → %s\n",
              result$Homogeneity$Levene$W,
              result$Homogeneity$Levene$p,
              result$Homogeneity$Levene$homogeneous))
  cat(sprintf("Bartlett 检验            : K² = %s, p = %s → %s\n\n",
              result$Homogeneity$Bartlett$K2,
              result$Homogeneity$Bartlett$p,
              result$Homogeneity$Bartlett$homogeneous))

  cat("【3. 决策建议】\n")
  cat(result$Decision$summary, "\n")

  invisible(result)
}
