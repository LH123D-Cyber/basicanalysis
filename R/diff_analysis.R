#' 差异分析及事后检验
#'
#' 根据指定的差异检验方法和事后检验方法，对目标变量进行组间比较，
#' 并生成包含均值、标准差、标准误和差异字母（CLD）的绘图数据框。
#'
#' @param data 数据框。
#' @param Objective 目标变量列名（字符型，只能选择一个）。
#' @param group 分组列名（字符型，只能选择一个）。
#' @param Diff 差异检验方法，可选 `"ANOVA"`、`"Welch ANOVA"`、`"Kruskal-Wallis"`、`"t"`、`"Welch t"`、`"Wilcoxon"`。
#' @param Posthoc 事后检验方法，如 `"Tukey HSD"`、`"Games-Howell"`、`"Dunn"`；不指定时根据 `Diff` 自动选择。
#' @param Adjust 事后检验 p 值校正方法，默认 `"BH"`。
#' @param summary 逻辑值，是否打印汇总信息，默认 `TRUE`。
#' @param alpha 显著性水平，默认 0.05。
#'
#' @return 一个 list，`summary_data` 为绘图数据框（含 group、mean、sd、se、letter），
#'   其余元素包括 `overall`、`posthoc`、`p_mat` 等。
#' @export
DiffAnalysis <- function(data, Objective, group, Diff, Posthoc = NULL,
                         Adjust = "BH", summary = TRUE, alpha = 0.05) {

  if (length(Objective) != 1) stop("Objective 只能选择一个变量")
  if (length(group) != 1) stop("group 只能指定一个分组列")
  if (!Objective %in% colnames(data)) stop("变量不存在: ", Objective)
  if (!group %in% colnames(data)) stop("分组列不存在: ", group)

  x <- data[[Objective]]
  g <- as.factor(data[[group]])
  if (!is.numeric(x)) stop("Objective 必须是数值型")

  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]; g <- droplevels(g[ok])

  groups <- levels(g); k <- length(groups)
  if (k < 2) stop("分组需要至少 2 个水平")

  diff_key <- tolower(trimws(Diff))

  # ============================================================
  # 1. 总体差异分析
  # ============================================================
  if (diff_key %in% c("anova", "one-way anova", "oneway")) {
    fit <- aov(x ~ g); sm <- summary(fit)[[1]]
    overall <- list(method = "One-way ANOVA",
                    stat = sm$`F value`[1], p = sm$`Pr(>F)`[1],
                    df1 = sm$Df[1], df2 = sm$Df[2])

  } else if (diff_key %in% c("welch anova", "welch", "welch's anova")) {
    wt <- oneway.test(x ~ g, var.equal = FALSE)
    overall <- list(method = "Welch ANOVA",
                    stat = unname(wt$statistic), p = wt$p.value,
                    df1 = unname(wt$parameter[1]), df2 = unname(wt$parameter[2]))

  } else if (diff_key %in% c("kruskal-wallis", "kruskal", "kw", "kruskal wallis")) {
    kt <- kruskal.test(x ~ g)
    overall <- list(method = "Kruskal-Wallis",
                    stat = unname(kt$statistic), p = kt$p.value,
                    df = unname(kt$parameter))

  } else if (diff_key %in% c("t", "t.test", "t-test", "student's t")) {
    if (k != 2) stop("t 检验需要恰好 2 个分组")
    gd <- split(x, g)
    tt <- t.test(gd[[1]], gd[[2]], var.equal = TRUE)
    overall <- list(method = "独立样本 t 检验",
                    stat = unname(tt$statistic), p = tt$p.value)

  } else if (diff_key %in% c("welch t", "welch t-test", "welch t.test")) {
    if (k != 2) stop("Welch t 检验需要恰好 2 个分组")
    gd <- split(x, g)
    tt <- t.test(gd[[1]], gd[[2]], var.equal = FALSE)
    overall <- list(method = "Welch t 检验",
                    stat = unname(tt$statistic), p = tt$p.value)

  } else if (diff_key %in% c("wilcoxon", "wilcox", "mann-whitney",
                             "mann whitney", "mann-whitney u")) {
    if (k != 2) stop("Wilcoxon 检验需要恰好 2 个分组")
    gd <- split(x, g)
    wt <- wilcox.test(gd[[1]], gd[[2]])
    overall <- list(method = "Wilcoxon 秩和检验",
                    stat = unname(wt$statistic), p = wt$p.value)

  } else {
    stop("不支持的 Diff 方法: ", Diff,
         "\n可用: ANOVA / Welch ANOVA / Kruskal-Wallis / t / Welch t / Wilcoxon")
  }

  # ============================================================
  # 2. 事后检验（k >= 3）
  # ============================================================
  p_mat <- NULL
  posthoc_info <- NULL

  if (k >= 3) {
    ph_key <- if (is.null(Posthoc)) "" else tolower(trimws(Posthoc))
    if (ph_key == "") {
      if (diff_key %in% c("anova", "one-way anova", "oneway")) ph_key <- "tukey hsd"
      else if (diff_key %in% c("welch anova", "welch", "welch's anova")) ph_key <- "games-howell"
      else ph_key <- "dunn"
    }
    ph_res <- run_posthoc(x, g, groups, diff_key, ph_key, Adjust)
    p_mat <- ph_res$p_mat
    posthoc_info <- ph_res$method

  } else {
    p_mat <- matrix(1, 2, 2, dimnames = list(groups, groups))
    p_mat[1, 2] <- p_mat[2, 1] <- overall$p
  }

  # ============================================================
  # 3. 描述性统计与差异字母
  # ============================================================
  gd <- split(x, g)
  desc <- data.frame(
    group = groups,
    n     = as.integer(sapply(gd, length)),
    mean  = sapply(gd, mean),
    sd    = sapply(gd, sd),
    se    = sapply(gd, sd) / sqrt(sapply(gd, length)),
    stringsAsFactors = FALSE
  )
  desc$letter <- assign_letters(desc$mean, p_mat, alpha, desc$group)

  # ============================================================
  # 4. 输出
  # ============================================================
  result <- list(
    summary_data = desc,
    overall      = overall,
    posthoc      = posthoc_info,
    p_mat        = p_mat,
    alpha        = alpha,
    variable     = Objective,
    group_var    = group
  )

  if (summary) {
    cat("══════════ DiffAnalysis 结果 ══════════\n")
    cat("变量:", Objective, "| 分组:", group, "| 组数:", k, "\n\n")
    cat(sprintf("【差异分析】%s\n  统计量 = %.4f, p = %.4f\n\n",
                overall$method, overall$stat, overall$p))
    if (!is.null(posthoc_info)) {
      cat("【事后检验】", posthoc_info, "\n")
      cat("两两 p 值矩阵:\n"); print(round(p_mat, 4)); cat("\n")
    }
    cat("【绘图数据框】\n")
    print(desc[, c("group", "mean", "sd", "se", "letter")], row.names = FALSE)
  }

  invisible(result)
}


# ============================================================
# 执行事后检验
# ============================================================
run_posthoc <- function(x, g, groups, diff_key, ph_key, Adjust) {
  k <- length(groups)
  p_mat <- matrix(1, k, k, dimnames = list(groups, groups))
  method_name <- ""

  if (ph_key %in% c("tukey hsd", "tukey", "tukeyhsd", "tukey's hsd")) {
    if (diff_key %in% c("anova", "one-way anova", "oneway")) {
      th <- TukeyHSD(aov(x ~ g))[[1]]
      method_name <- "Tukey HSD"
      for (i in seq_len(nrow(th))) {
        pn <- rownames(th)[i]
        for (g1 in groups) for (g2 in groups) {
          if (g1 != g2 && pn == paste0(g1, "-", g2)) {
            p_mat[g1, g2] <- p_mat[g2, g1] <- th[i, "p adj"]
            break
          }
        }
      }
    } else {
      warning("Tukey HSD 仅适用于 ANOVA，已回退到成对 t 检验")
      ph <- pairwise.t.test(x, g, p.adjust.method = Adjust)
      method_name <- paste0("成对 t 检验 (", Adjust, ")")
      fill_pmat_from_pairwise(ph, groups, p_mat)
    }

  } else if (ph_key %in% c("games-howell", "games howell", "gh")) {
    if (requireNamespace("rstatix", quietly = TRUE)) {
      df_g <- data.frame(x = x, g = g)
      gh <- rstatix::games_howell_test(df_g, x ~ g)
      method_name <- "Games-Howell"
      for (i in seq_len(nrow(gh))) {
        g1 <- as.character(gh$group1[i]); g2 <- as.character(gh$group2[i])
        p_mat[g1, g2] <- p_mat[g2, g1] <- gh$p.adj[i]
      }
    } else {
      warning("未安装 rstatix，已回退到成对 Welch t 检验")
      ph <- pairwise.t.test(x, g, pool.sd = FALSE, p.adjust.method = Adjust)
      method_name <- paste0("成对 Welch t (", Adjust, ")")
      fill_pmat_from_pairwise(ph, groups, p_mat)
    }

  } else if (ph_key %in% c("dunn", "dunn's test", "dunn test")) {
    if (requireNamespace("dunn.test", quietly = TRUE)) {
      dt <- dunn.test::dunn.test(x, g, method = Adjust, table = FALSE, kw = FALSE)
      method_name <- paste0("Dunn 检验 (", Adjust, ")")
      for (i in seq_along(dt$comparisons)) {
        parts <- strsplit(dt$comparisons[i], " - ")[[1]]
        if (length(parts) == 2) {
          g1 <- trimws(parts[1]); g2 <- trimws(parts[2])
          if (g1 %in% groups && g2 %in% groups)
            p_mat[g1, g2] <- p_mat[g2, g1] <- dt$altp[i]
        }
      }
    } else if (requireNamespace("FSA", quietly = TRUE)) {
      dt <- FSA::dunnTest(x ~ g, data = data.frame(x = x, g = g), method = Adjust)
      method_name <- paste0("Dunn 检验 (", Adjust, ")")
      for (i in seq_len(nrow(dt$res))) {
        parts <- strsplit(rownames(dt$res)[i], " - ")[[1]]
        if (length(parts) == 2) {
          g1 <- trimws(parts[1]); g2 <- trimws(parts[2])
          if (g1 %in% groups && g2 %in% groups)
            p_mat[g1, g2] <- p_mat[g2, g1] <- dt$res$P.adj[i]
        }
      }
    } else {
      warning("未安装 dunn.test/FSA，已回退到成对 Wilcoxon 检验")
      ph <- pairwise.wilcox.test(x, g, p.adjust.method = Adjust, exact = FALSE)
      method_name <- paste0("成对 Wilcoxon (", Adjust, ")")
      fill_pmat_from_pairwise(ph, groups, p_mat)
    }

  } else {
    # 其他名称（LSD/Bonferroni 等）→ 通用成对检验
    if (diff_key %in% c("kruskal-wallis", "kruskal", "kw")) {
      ph <- pairwise.wilcox.test(x, g, p.adjust.method = Adjust, exact = FALSE)
      method_name <- paste0("成对 Wilcoxon (", Adjust, ")")
    } else {
      ph <- pairwise.t.test(x, g, p.adjust.method = Adjust)
      method_name <- paste0("成对 t 检验 (", Adjust, ")")
    }
    fill_pmat_from_pairwise(ph, groups, p_mat)
  }

  list(p_mat = p_mat, method = method_name)
}

# 从 pairwise.t.test / pairwise.wilcox.test 的结果填充对称 p 值矩阵
fill_pmat_from_pairwise <- function(ph, groups, p_mat) {
  pv <- ph$p.value
  if (is.null(dim(pv))) return(invisible(NULL))
  for (i in seq_len(nrow(pv))) for (j in seq_len(ncol(pv))) {
    r <- rownames(pv)[i]; c <- colnames(pv)[j]
    if (!is.na(pv[i, j]) && r %in% groups && c %in% groups) {
      p_mat[r, c] <- pv[i, j]
      p_mat[c, r] <- pv[i, j]
    }
  }
  invisible(NULL)
}

# ============================================================
# 差异字母（CLD）：基于 p 值矩阵的贪心插入算法
# ============================================================
assign_letters <- function(means, p_mat, alpha = 0.05, labels = NULL) {
  n <- length(means)
  if (n == 1) return("a")
  if (is.null(labels)) labels <- rownames(p_mat)
  if (is.null(labels)) labels <- as.character(seq_len(n))

  ord <- order(means, decreasing = TRUE)
  same <- p_mat > alpha; diag(same) <- TRUE
  same_s <- same[ord, ord]

  letter_pool <- c(letters, LETTERS)
  cld <- rep("", n); curr <- 1L

  for (i in seq_len(n)) {
    found <- FALSE
    for (L in letter_pool) {
      uses_L <- which(sapply(strsplit(cld, ""), function(s) L %in% s))
      if (length(uses_L) == 0) next
      if (all(same_s[i, uses_L])) { cld[i] <- L; found <- TRUE; break }
    }
    if (!found) { cld[i] <- letter_pool[curr]; curr <- curr + 1L }
  }

  # 修复“应共享字母但未共享”的情况
  changed <- TRUE; iter <- 0
  while (changed && iter < 100) {
    changed <- FALSE; iter <- iter + 1
    for (i in seq_len(n - 1)) for (j in (i + 1):n) {
      if (!same_s[i, j]) next
      li <- strsplit(cld[i], "")[[1]]
      lj <- strsplit(cld[j], "")[[1]]
      if (length(intersect(li, lj)) == 0) {
        for (L in sort(li)) {
          uses_L <- which(sapply(strsplit(cld, ""), function(s) L %in% s))
          if (all(same_s[j, uses_L])) {
            cld[j] <- paste0(cld[j], L); changed <- TRUE; break
          }
        }
      }
    }
  }

  out <- character(n); out[ord] <- cld
  names(out) <- labels
  out
}
