# ============================================================
# basicanalysis 完整单文件版（code.R）
# 用途：备份 / 单文件 source 测试 / GitHub 存档
# 注意：若包已使用 R/ 下的拆分文件，请勿把本文件放入 R/
# ============================================================

# ============================================================
# 整体正态检验
# ============================================================
OND <- function(col,
                QQ = FALSE,
                His = FALSE,
                ShapiroWilk = FALSE,
                Convert = FALSE,
                plot = FALSE,
                data = df,
                outdir = NULL) {

  if (!all(col %in% colnames(data))) {
    stop("以下列名不存在于数据中: ",
         paste(setdiff(col, colnames(data)), collapse = ", "))
  }

  sub <- data[, col, drop = FALSE]
  result <- list()
  result$col  <- col
  result$data_source <- deparse(substitute(data))

  get_mfrow <- function(n) {
    if (n <= 2) return(c(1, n))
    return(c(ceiling(n / 2), 2))
  }

  yeojohnson <- function(x, lambda) {
    out <- numeric(length(x))
    for (i in seq_along(x)) {
      xi <- x[i]
      if (xi >= 0) {
        out[i] <- if (lambda == 0) log(xi + 1) else (xi^lambda - 1) / lambda
      } else {
        if (lambda == 2) out[i] <- -log(-xi + 1)
        else out[i] <- -(((-xi + 1)^(2 - lambda) - 1) / (2 - lambda))
      }
    }
    out
  }

  boxcox_tf <- function(x, lambda) {
    if (abs(lambda) < 1e-10) log(x) else (x^lambda - 1) / lambda
  }

  make_plot_obj <- function(draw_fun, mfrow_par, name) {
    if (!is.null(outdir)) {
      if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
      f <- file.path(outdir, paste0(name, ".png"))
      png(f, width = 1200, height = 800, res = 120)
      par(mfrow = mfrow_par)
      draw_fun()
      dev.off()
      return(list(file = f, mfrow = mfrow_par, is_file = TRUE))
    }
    grDevices::dev.new()
    par(mfrow = mfrow_par)
    draw_fun()
    list(mfrow = mfrow_par, is_file = FALSE, draw_fun = draw_fun)
  }

  if (ShapiroWilk) {
    sw_list <- list()
    for (v in col) {
      sw_list[[v]] <- tryCatch(shapiro.test(sub[[v]]),
                               error = function(e) paste("Error:", e$message))
    }
    result$ShapiroWilk <- sw_list
  }

  if (QQ) {
    result$QQ1 <- make_plot_obj(
      draw_fun = function() {
        for (v in col) {
          qqnorm(sub[[v]], main = paste("Q-Q:", v), pch = 19, col = "blue")
          qqline(sub[[v]], col = "red", lwd = 2)
        }
      },
      mfrow_par = get_mfrow(length(col)),
      name = "QQ1"
    )
  }

  if (His) {
    result$His1 <- make_plot_obj(
      draw_fun = function() {
        for (v in col) {
          hist(sub[[v]], main = paste("Histogram:", v), xlab = v,
               col = "lightblue", border = "white")
        }
      },
      mfrow_par = get_mfrow(length(col)),
      name = "His1"
    )
  }

  if (Convert) {
    convert_results <- data.frame()
    converted_data   <- list()
    full_data        <- data

    for (v in col) {
      x <- sub[[v]]
      candidates <- c("Yeo-Johnson", "x^2")
      if (all(x > 0))  candidates <- c(candidates, "log", "Box-Cox")
      if (all(x >= 0)) candidates <- c(candidates, "sqrt")
      if (all(x != 0)) candidates <- c(candidates, "1/x")

      best <- list(p = -Inf, method = NA_character_, data = NULL,
                   lambda = NA_real_, W = NA_real_)

      for (m in candidates) {
        if (m %in% c("Box-Cox", "Yeo-Johnson")) {
          for (lam in seq(-2, 2, by = 0.1)) {
            xt <- if (m == "Box-Cox") boxcox_tf(x, lam) else yeojohnson(x, lam)
            if (any(!is.finite(xt))) next
            sw <- tryCatch(shapiro.test(xt), error = function(e) NULL)
            if (is.null(sw) || !is.finite(sw$p.value)) next
            if (sw$p.value > best$p) {
              best <- list(p = sw$p.value, method = m, data = xt,
                           lambda = lam, W = as.numeric(sw$statistic))
            }
          }
        } else {
          xt <- switch(m, "log" = log(x), "sqrt" = sqrt(x),
                       "1/x" = 1/x, "x^2" = x^2)
          if (any(!is.finite(xt))) next
          sw <- tryCatch(shapiro.test(xt), error = function(e) NULL)
          if (is.null(sw) || !is.finite(sw$p.value)) next
          if (sw$p.value > best$p) {
            best <- list(p = sw$p.value, method = m, data = xt,
                         lambda = NA_real_, W = as.numeric(sw$statistic))
          }
        }
      }

      convert_results <- rbind(convert_results, data.frame(
        variable     = v,
        method       = best$method,
        lambda       = round(best$lambda, 2),
        W            = round(best$W, 4),
        p.value      = round(best$p, 4),
        normal_after = ifelse(best$p > 0.05, "Yes", "No"),
        stringsAsFactors = FALSE
      ))

      if (best$p > 0.05) {
        converted_data[[v]] <- best$data
        full_data[[v]]     <- best$data
      }
    }

    result$Convert <- convert_results
    result$data    <- full_data
    attr(result$data, "converted_vars") <- names(converted_data)

    if (plot && length(converted_data) > 0) {
      n <- length(converted_data)
      result$QQ2 <- make_plot_obj(
        draw_fun = function() {
          for (v in names(converted_data)) {
            qqnorm(converted_data[[v]],
                   main = paste("Q-Q (converted):", v),
                   pch = 19, col = "darkgreen")
            qqline(converted_data[[v]], col = "red", lwd = 2)
          }
        },
        mfrow_par = get_mfrow(n), name = "QQ2"
      )
      result$His2 <- make_plot_obj(
        draw_fun = function() {
          for (v in names(converted_data)) {
            hist(converted_data[[v]],
                 main = paste("Histogram (converted):", v),
                 xlab = v, col = "lightgreen", border = "white")
          }
        },
        mfrow_par = get_mfrow(n), name = "His2"
      )
    }
  }

  return(result)
}


# ============================================================
# 显示辅助函数：显示或重绘 OND 生成的图形
# ============================================================
OND_show <- function(result, name, data = NULL) {
  obj <- result[[name]]
  if (is.null(obj)) {
    message(name, " 不存在（可能未生成）")
    return(invisible(NULL))
  }

  if (!is.null(data)) {
    col <- result$col
    if (!all(col %in% colnames(data))) {
      stop("新数据中缺少这些列: ",
           paste(setdiff(col, colnames(data)), collapse = ", "))
    }
    sub <- data[, col, drop = FALSE]
    grDevices::dev.new()
    par(mfrow = obj$mfrow)

    if (name %in% c("QQ1", "QQ2")) {
      for (v in col) {
        qqnorm(sub[[v]], main = paste("Q-Q:", v), pch = 19, col = "blue")
        qqline(sub[[v]], col = "red", lwd = 2)
      }
    } else if (name %in% c("His1", "His2")) {
      for (v in col) {
        hist(sub[[v]], main = paste("Histogram:", v), xlab = v,
             col = "lightblue", border = "white")
      }
    }
    return(invisible(NULL))
  }

  if (isTRUE(obj$is_file)) {
    os <- .Platform$OS.type
    if (os == "windows") shell.exec(obj$file)
    else if (Sys.info()[["sysname"]] == "Darwin")
      system(paste("open", shQuote(obj$file)))
    else
      system(paste("xdg-open", shQuote(obj$file)))
  } else {
    grDevices::dev.new()
    par(mfrow = obj$mfrow)
    obj$draw_fun()
  }
  invisible(NULL)
}


# ============================================================
# 统计方法判定
# ============================================================
DDT <- function(data,
                Objective,
                group,
                alpha = 0.05) {

  if (length(Objective) != 1) stop("Objective 只能选择一个变量")
  if (length(group)     != 1) stop("group 只能指定一个分组列")
  if (!Objective %in% colnames(data))
    stop("Objective 变量不存在于数据中: ", Objective)
  if (!group %in% colnames(data))
    stop("group 列不存在于数据中: ", group)

  x <- data[[Objective]]
  g <- data[[group]]

  if (!is.numeric(x)) stop("Objective 必须是数值型变量")
  if (!is.factor(g))  g <- factor(g)

  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]; g <- droplevels(g[ok])

  groups   <- levels(g)
  grp_data <- split(x, g)
  k        <- length(groups)

  if (k < 2) stop("分组列至少需要 2 个水平")

  result <- list(
    variable = Objective,
    group    = group,
    alpha    = alpha,
    n_groups = k,
    groups   = groups,
    n_each   = sapply(grp_data, length)
  )

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
        norm_test$W[i]      <- round(as.numeric(sw$statistic), 4)
        norm_test$p[i]      <- round(sw$p.value, 4)
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
    med <- tapply(x, g, median, na.rm = TRUE)
    z   <- abs(x - med[as.character(g)])
    lev_fit <- aov(z ~ g)
    lev_tab <- summary(lev_fit)[[1]]
    lev_W <- lev_tab$`F value`[1]
    lev_p <- lev_tab$`Pr(>F)`[1]
  }

  brt <- tryCatch(
    bartlett.test(x ~ g, data = data.frame(x = x, g = g)),
    error = function(e) NULL
  )
  brt_K2 <- if (!is.null(brt)) as.numeric(brt$statistic) else NA_real_
  brt_p  <- if (!is.null(brt)) brt$p.value                 else NA_real_

  result$Homogeneity <- list(
    Levene   = list(W = round(lev_W, 4),
                    p = round(lev_p, 4),
                    homogeneous = ifelse(!is.na(lev_p) && lev_p > alpha,
                                         "Yes", "No")),
    Bartlett = list(K2 = round(brt_K2, 4),
                    p  = round(brt_p, 4),
                    homogeneous = ifelse(!is.na(brt_p) && brt_p > alpha,
                                         "Yes", "No"))
  )
  var_homogeneous <- !is.na(lev_p) && lev_p > alpha

  if (k == 2) {
    if (all_normal && var_homogeneous) {
      method  <- "独立样本 t 检验"
      posthoc <- "无（两组无需事后检验）"
      reason  <- "两组数据均满足正态分布，且方差齐性"
    } else if (all_normal && !var_homogeneous) {
      method  <- "Welch t 检验"
      posthoc <- "无（两组无需事后检验）"
      reason  <- "两组数据均满足正态分布，但方差不齐"
    } else {
      method  <- "Wilcoxon 秩和检验"
      posthoc <- "无（两组无需事后检验）"
      reason  <- "数据不满足正态分布（Wilcoxon 对分布和方差不敏感）"
    }
  } else {
    if (all_normal && var_homogeneous) {
      method  <- "One-way ANOVA（单因素方差分析）"
      posthoc <- "Tukey HSD"
      reason  <- "各组数据满足正态分布，且方差齐性"
    } else if (all_normal && !var_homogeneous) {
      method  <- "Welch ANOVA（Welch 单因素方差分析）"
      posthoc <- "Games-Howell"
      reason  <- "各组数据满足正态分布，但方差齐性"
    } else if (!all_normal && var_homogeneous) {
      method  <- "Kruskal-Wallis 秩和检验"
      posthoc <- "Dunn 检验（BH 校正）"
      reason  <- "数据不满足正态分布，但方差齐性"
    } else {
      method  <- "Kruskal-Wallis 秩和检验"
      posthoc <- "Dunn 检验（BH 校正）"
      reason  <- "数据不满足正态分布，且方差不齐（Kruskal-Wallis 对两者均不敏感）"
    }
  }

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


# ============================================================
# 差异分析及事后检验
# ============================================================
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

  # 1. 总体差异分析
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

  # 2. 事后检验（k >= 3）
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

  # 3. 描述性统计与差异字母
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

  # 4. 输出
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
# 执行事后检验（内部函数）
# ============================================================
run_posthoc <- function(x, g, groups, diff_key, ph_key, Adjust) {
  Adjust <- tolower(Adjust)

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
      # list = TRUE 才会返回列表；调整后 p 值字段名为 P.adjusted
      dt <- dunn.test::dunn.test(x, g, method = Adjust,
                                 table = FALSE, kw = FALSE, list = TRUE)
      method_name <- paste0("Dunn 检验 (", Adjust, ")")
      for (i in seq_along(dt$comparisons)) {
        parts <- strsplit(dt$comparisons[i], " - ")[[1]]
        if (length(parts) == 2) {
          g1 <- trimws(parts[1]); g2 <- trimws(parts[2])
          if (g1 %in% groups && g2 %in% groups)
            p_mat[g1, g2] <- p_mat[g2, g1] <- dt$P.adjusted[i]
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


# ============================================================
# 从 pairwise 检验结果填充对称 p 值矩阵（内部函数）
# ============================================================
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


# ============================================================
# 辅助 1：整理散点数据（支持独立的颜色分组）
# ============================================================
.prepare_scatter <- function(Original, pointgroup, pointcolor_by, raw_data,
                             mean, target_levels) {
  if (!is.null(Original)) {
    if (is.list(Original)) {
      if (!is.null(pointcolor_by))
        warning("Original 为 list 时 pointcolor_by 被忽略；",
                "若需散点独立颜色分组，请用数值向量 + pointgroup + pointcolor_by")
      if (is.null(names(Original)) || any(names(Original) == ""))
        stop("Original 为 list 时，每个元素必须有名字（组名 = 数值向量）")
      df <- do.call(rbind, lapply(names(Original), function(g) {
        data.frame(group = g, .value = Original[[g]], stringsAsFactors = FALSE)
      }))
    } else if (is.numeric(Original)) {
      if (is.null(pointgroup))
        stop("Original 为数值向量时必须提供 pointgroup")
      if (length(pointgroup) != length(Original))
        stop("Original 与 pointgroup 长度不一致")
      df <- data.frame(group = as.character(pointgroup),
                       .value = Original, stringsAsFactors = FALSE)
      if (!is.null(pointcolor_by)) {
        if (length(pointcolor_by) != length(Original))
          stop("Original 与 pointcolor_by 长度不一致")
        df$.point_color <- as.character(pointcolor_by)
      }
    } else {
      stop("Original 需为命名 list 或数值向量")
    }
  } else if (!is.null(raw_data)) {
    if (!is.data.frame(raw_data)) stop("raw_data 必须是数据框")
    if (!"group" %in% colnames(raw_data)) stop("raw_data 缺少 group 列")
    if (!mean %in% colnames(raw_data))
      stop("raw_data 缺少 '", mean, "' 列")
    df <- data.frame(group = as.character(raw_data$group),
                     .value = raw_data[[mean]], stringsAsFactors = FALSE)
    if (!is.null(pointcolor_by)) {
      if (!is.character(pointcolor_by) || length(pointcolor_by) != 1)
        stop("raw_data 配合 pointcolor_by 时，pointcolor_by 应为 raw_data 中的列名")
      if (!pointcolor_by %in% colnames(raw_data))
        stop("raw_data 中不存在列: ", pointcolor_by)
      df$.point_color <- as.character(raw_data[[pointcolor_by]])
    }
  } else {
    return(NULL)
  }

  src_levels <- unique(df$group)
  common     <- intersect(src_levels, target_levels)
  if (length(common) == 0) {
    stop("散点位置分组与 data$group 的水平完全不一致！\n",
         "  · data$group 水平 : ", paste(target_levels, collapse = ", "), "\n",
         "  · 传入分组水平    : ", paste(src_levels,    collapse = ", "), "\n",
         "→ 请检查 pointgroup 是否与 data$group 使用同一套分组")
  }
  if (length(common) < length(src_levels)) {
    warning("散点分组中有 ", length(src_levels) - length(common),
            " 个水平不在 data$group 中，已忽略: ",
            paste(setdiff(src_levels, target_levels), collapse = ", "))
  }
  df$group <- factor(df$group, levels = target_levels)
  df <- df[!is.na(df$group), , drop = FALSE]
  df
}


# ============================================================
# 辅助 2：添加散点图层（优先按 .point_color 上色）
# ============================================================
.add_scatter <- function(p, scatter_df, pointcolor, color, group_levels,
                         pointcolor_name = NULL) {
  if (is.null(scatter_df) || nrow(scatter_df) == 0) return(p)

  has_point_color <- ".point_color" %in% colnames(scatter_df)

  # 场景 A：散点颜色按 pointcolor_by 分组
  if (has_point_color) {
    pc_levels <- unique(scatter_df$.point_color)

    if (!is.null(pointcolor) && length(pointcolor) > 1) {
      if (!is.null(names(pointcolor))) {
        missing <- setdiff(pc_levels, names(pointcolor))
        if (length(missing) > 0) {
          stop("pointcolor 缺少以下散点颜色分组的颜色: ",
               paste(missing, collapse = ", "), "\n",
               "  · 散点颜色分组水平 : ", paste(pc_levels, collapse = ", "), "\n",
               "  · names(pointcolor) : ", paste(names(pointcolor), collapse = ", "), "\n",
               "→ 请核对 pointcolor 的命名是否与 pointcolor_by 的水平一致")
        }
        values <- pointcolor[pc_levels]
      } else {
        if (length(pointcolor) != length(pc_levels))
          stop("pointcolor 为无名向量时，长度必须等于散点颜色分组水平数 (",
               length(pc_levels), ")")
        values <- setNames(pointcolor, pc_levels)
      }
    } else if (!is.null(pointcolor) && length(pointcolor) == 1) {
      values <- setNames(rep(pointcolor, length(pc_levels)), pc_levels)
    } else {
      values <- NULL
    }

    scatter_df$.point_color <- factor(scatter_df$.point_color, levels = pc_levels)
    legend_name <- if (is.null(pointcolor_name)) "颜色" else pointcolor_name

    if (is.null(values)) {
      p <- p + geom_jitter(data = scatter_df,
                           aes(x = group, y = .value, color = .point_color),
                           width = 0.15, size = 1.5, alpha = 0.5,
                           inherit.aes = FALSE) +
        labs(color = legend_name)
    } else {
      p <- p + geom_jitter(data = scatter_df,
                           aes(x = group, y = .value, color = .point_color),
                           width = 0.15, size = 1.5, alpha = 0.5,
                           inherit.aes = FALSE) +
        scale_color_manual(values = values, name = legend_name)
    }
    return(p)
  }

  # 场景 B：无 pointcolor_by，沿用旧逻辑
  use_group_color <- FALSE
  color_values    <- NULL

  if (!is.null(pointcolor)) {
    if (length(pointcolor) > 1) {
      if (!is.null(names(pointcolor))) {
        missing_g <- setdiff(group_levels, names(pointcolor))
        if (length(missing_g) > 0)
          stop("pointcolor 缺少以下分组的颜色: ",
               paste(missing_g, collapse = ", "))
        color_values <- pointcolor[group_levels]
      } else {
        if (length(pointcolor) != length(group_levels))
          stop("pointcolor 为无名向量时，长度必须等于分组数")
        color_values <- setNames(pointcolor, group_levels)
      }
      use_group_color <- TRUE
    } else {
      color_values <- pointcolor
    }
  } else if (length(color) > 1) {
    if (!is.null(names(color))) color_values <- color[group_levels]
    else color_values <- setNames(color, group_levels)
    use_group_color <- TRUE
  } else {
    color_values <- color
  }

  if (use_group_color) {
    p + geom_jitter(data = scatter_df,
                    aes(x = group, y = .value, color = group),
                    width = 0.15, size = 1.5, alpha = 0.5,
                    inherit.aes = FALSE) +
      scale_color_manual(values = color_values, guide = "none")
  } else {
    p + geom_jitter(data = scatter_df,
                    aes(x = group, y = .value),
                    width = 0.15, size = 1.5, alpha = 0.5,
                    color = color_values, inherit.aes = FALSE)
  }
}


# ============================================================
# Barplot：普通差异柱状图
# ============================================================
Barplot <- function(data, mean = "mean", Whiskers = "se",
                    point = FALSE, x = NULL, color = "steelblue",
                    alpha = 0.7, label = TRUE, loca = NULL,
                    raw_data = NULL,
                    Original = NULL, pointgroup = NULL,
                    pointcolor = NULL, pointcolor_by = NULL,
                    pointcolor_name = NULL,
                    bargroup = "group") {

  requireNamespace("ggplot2", quietly = TRUE) ||
    stop("请先安装 ggplot2: install.packages('ggplot2')")
  library(ggplot2)

  if (!is.data.frame(data)) stop("data 必须是数据框")
  if (!mean %in% colnames(data)) stop("mean 列不存在: ", mean)
  if (!"group" %in% colnames(data)) stop("数据中必须有 'group' 列")
  if (!Whiskers %in% c("sd", "se")) stop("Whiskers 只能是 'sd' 或 'se'")
  if (!Whiskers %in% colnames(data)) stop("数据缺少 ", Whiskers, " 列")
  if (label && !"letter" %in% colnames(data)) stop("label=TRUE 但缺少 letter 列")
  if (!is.character(bargroup) || length(bargroup) != 1)
    stop("bargroup 必须是单个列名字符串")
  if (!bargroup %in% colnames(data))
    stop("data 中不存在柱状图颜色分组列: ", bargroup)

  if (!is.null(x)) {
    data <- data[data$group %in% x, , drop = FALSE]
    data$group <- factor(data$group, levels = x)
  } else {
    data$group <- factor(data$group, levels = data$group)
  }
  data <- data[order(data$group), , drop = FALSE]
  group_levels <- levels(data$group)

  data$ymin <- data[[mean]] - data[[Whiskers]]
  data$ymax <- data[[mean]] + data[[Whiskers]]
  if (is.null(loca)) loca <- diff(range(c(data$ymin, data$ymax))) * 0.05

  if (bargroup != "group") {
    chk <- tapply(as.character(data[[bargroup]]), as.character(data$group),
                  function(v) length(unique(v)))
    if (any(chk > 1))
      warning("部分 x 分组在 '", bargroup,
              "' 列中有多个取值，柱子填充色可能不一致")
  }

  scatter_df <- if (point) .prepare_scatter(Original, pointgroup, pointcolor_by,
                                            raw_data, mean, group_levels) else NULL

  multi_color <- length(color) > 1
  if (multi_color) {
    bar_levels <- if (bargroup == "group") group_levels
    else unique(as.character(data[[bargroup]]))

    if (!is.null(names(color))) {
      miss <- setdiff(bar_levels, names(color))
      if (length(miss) > 0)
        stop("color 缺少以下 '", bargroup, "' 水平的颜色: ",
             paste(miss, collapse = ", "))
      fill_values <- color[bar_levels]
    } else {
      if (length(color) != length(bar_levels))
        stop("color 为无名向量时长度必须等于 '", bargroup, "' 的水平数 (",
             length(bar_levels), ")")
      fill_values <- setNames(color, bar_levels)
    }
  }

  if (multi_color) {
    p <- ggplot(data, aes(x = group, y = .data[[mean]],
                          fill = .data[[bargroup]])) +
      geom_col(alpha = alpha, color = "black", width = 0.65) +
      scale_fill_manual(values = fill_values, guide = "none")
  } else {
    p <- ggplot(data, aes(x = group, y = .data[[mean]])) +
      geom_col(fill = color, alpha = alpha, color = "black", width = 0.65)
  }

  p <- p +
    geom_errorbar(aes(ymin = ymin, ymax = ymax),
                  color = "black", width = 0.2, linewidth = 0.8) +
    labs(x = NULL, y = mean) +
    theme_bw(base_size = 13) +
    theme(axis.text = element_text(color = "black", size = 12),
          panel.grid = element_blank())

  if (label) {
    p <- p + geom_text(aes(label = letter, y = ymax + loca),
                       vjust = 0, size = 5)
  }

  if (point) {
    if (!is.null(scatter_df)) {
      p <- .add_scatter(p, scatter_df, pointcolor, color, group_levels, pointcolor_name)
    } else {
      p <- p + geom_point(aes(y = .data[[mean]]), size = 2, color = "white")
    }
  }

  return(p)
}


# ============================================================
# CircleBar：支持 circle=TRUE 的极坐标版本
# ============================================================
CircleBar <- function(data, mean = "mean", Whiskers = "se",
                      point = FALSE, x = NULL, color = "steelblue",
                      alpha = 0.7, label = TRUE, loca = NULL,
                      raw_data = NULL,
                      Original = NULL, pointgroup = NULL,
                      pointcolor = NULL, pointcolor_by = NULL,
                      pointcolor_name = NULL,
                      bargroup = "group",
                      circle = FALSE, Inter = 0) {

  requireNamespace("ggplot2", quietly = TRUE) ||
    stop("请先安装 ggplot2: install.packages('ggplot2')")
  library(ggplot2)

  if (!is.data.frame(data)) stop("data 必须是数据框")
  if (!mean %in% colnames(data)) stop("mean 列不存在: ", mean)
  if (!"group" %in% colnames(data)) stop("数据中必须有 'group' 列")
  if (!Whiskers %in% c("sd", "se")) stop("Whiskers 只能是 'sd' 或 'se'")
  if (!Whiskers %in% colnames(data)) stop("数据缺少 ", Whiskers, " 列")
  if (label && !"letter" %in% colnames(data)) stop("label=TRUE 但缺少 letter 列")
  if (circle && Inter < 0) stop("Inter 必须 >= 0")
  if (!is.character(bargroup) || length(bargroup) != 1)
    stop("bargroup 必须是单个列名字符串")
  if (!bargroup %in% colnames(data))
    stop("data 中不存在柱状图颜色分组列: ", bargroup)

  if (!is.null(x)) {
    data <- data[data$group %in% x, , drop = FALSE]
    data$group <- factor(data$group, levels = x)
  } else {
    data$group <- factor(data$group, levels = data$group)
  }
  data <- data[order(data$group), , drop = FALSE]
  group_levels <- levels(data$group)

  data$ymin <- data[[mean]] - data[[Whiskers]]
  data$ymax <- data[[mean]] + data[[Whiskers]]
  if (is.null(loca)) loca <- diff(range(c(data$ymin, data$ymax))) * 0.05

  if (bargroup != "group") {
    chk <- tapply(as.character(data[[bargroup]]), as.character(data$group),
                  function(v) length(unique(v)))
    if (any(chk > 1))
      warning("部分 x 分组在 '", bargroup,
              "' 列中有多个取值，柱子填充色可能不一致")
  }

  scatter_df <- if (point) .prepare_scatter(Original, pointgroup, pointcolor_by,
                                            raw_data, mean, group_levels) else NULL

  multi_color <- length(color) > 1
  if (multi_color) {
    bar_levels <- if (bargroup == "group") group_levels
    else unique(as.character(data[[bargroup]]))
    if (!is.null(names(color))) {
      miss <- setdiff(bar_levels, names(color))
      if (length(miss) > 0)
        stop("color 缺少以下 '", bargroup, "' 水平的颜色: ",
             paste(miss, collapse = ", "))
      fill_values <- color[bar_levels]
    } else {
      if (length(color) != length(bar_levels))
        stop("color 为无名向量时长度必须等于 '", bargroup, "' 的水平数 (",
             length(bar_levels), ")")
      fill_values <- setNames(color, bar_levels)
    }
  }

  # 普通柱状图（circle = FALSE）
  if (!circle) {
    if (multi_color) {
      p <- ggplot(data, aes(x = group, y = .data[[mean]],
                            fill = .data[[bargroup]])) +
        geom_col(alpha = alpha, color = "black", width = 0.65) +
        scale_fill_manual(values = fill_values, guide = "none")
    } else {
      p <- ggplot(data, aes(x = group, y = .data[[mean]])) +
        geom_col(fill = color, alpha = alpha, color = "black", width = 0.65)
    }
    p <- p +
      geom_errorbar(aes(ymin = ymin, ymax = ymax),
                    color = "black", width = 0.2, linewidth = 0.8) +
      labs(x = NULL, y = mean) +
      theme_bw(base_size = 13) +
      theme(axis.text = element_text(color = "black", size = 12),
            panel.grid = element_blank())

    if (label) {
      p <- p + geom_text(aes(label = letter, y = ymax + loca),
                         vjust = 0, size = 5)
    }
    if (point) {
      if (!is.null(scatter_df)) {
        p <- .add_scatter(p, scatter_df, pointcolor, color, group_levels, pointcolor_name)
      } else {
        p <- p + geom_point(aes(y = .data[[mean]]), size = 2, color = "white")
      }
    }
    return(p)
  }

  # 圆形柱状图（circle = TRUE）
  y_top_max <- max(data$ymax + loca, na.rm = TRUE)
  y_range   <- diff(range(c(data$ymin, data$ymax)))
  y_upper   <- y_top_max + y_range * 0.15

  if (multi_color) {
    p <- ggplot(data, aes(x = group, y = .data[[mean]],
                          fill = .data[[bargroup]])) +
      geom_col(alpha = alpha, color = "black", width = 0.65) +
      scale_fill_manual(values = fill_values, guide = "none")
  } else {
    p <- ggplot(data, aes(x = group, y = .data[[mean]])) +
      geom_col(fill = color, alpha = alpha, color = "black", width = 0.65)
  }

  p <- p +
    geom_errorbar(aes(ymin = ymin, ymax = ymax),
                  color = "black", width = 0.2, linewidth = 0.8) +
    coord_polar(start = 0) +
    ylim(-Inter, y_upper) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid  = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks  = element_blank(),
      axis.title  = element_blank(),
      axis.text.x = element_text(size = 12, color = "black")
    )

  if (label) {
    p <- p + geom_text(aes(label = letter, y = ymax + loca),
                       vjust = 0, size = 5)
  }

  if (point) {
    if (!is.null(scatter_df)) {
      p <- .add_scatter(p, scatter_df, pointcolor, color, group_levels, pointcolor_name)
    } else {
      p <- p + geom_point(aes(y = .data[[mean]]), size = 2, color = "white")
    }
  }

  return(p)
}


# ============================================================
# plotclass：不做差异分析时的 0/1 分类热图
# ============================================================
plotclass <- function(data, col, x, y,
                      color1 = "#cfe2f3",
                      color2 = "#e06666",
                      gap    = 1,
                      na_color = "grey90",
                      xlab = NULL, ylab = NULL,
                      legend_name = NULL) {

  requireNamespace("ggplot2", quietly = TRUE) ||
    stop("请先安装 ggplot2: install.packages('ggplot2')")
  library(ggplot2)

  if (!is.data.frame(data)) stop("data 必须是数据框")
  for (nm in c(col, x, y)) {
    if (!nm %in% colnames(data)) stop("数据中不存在列: ", nm)
  }

  df <- data.frame(
    .x    = factor(data[[x]]),
    .y    = factor(data[[y]]),
    .val  = data[[col]],
    stringsAsFactors = FALSE
  )

  bad <- setdiff(unique(df$.val[!is.na(df$.val)]), c(0, 1))
  if (length(bad) > 0)
    warning("列 '", col, "' 含有非 0/1 的值: ",
            paste(bad, collapse = ", "), "（将按 NA 处理）")
  df$.val[!is.na(df$.val) & !(df$.val %in% c(0, 1))] <- NA

  df$.val_chr <- ifelse(is.na(df$.val), NA,
                        ifelse(df$.val == 1, "1", "0"))
  df$.val_chr <- factor(df$.val_chr, levels = c("0", "1"))

  if (is.null(xlab)) xlab <- x
  if (is.null(ylab)) ylab <- y
  if (is.null(legend_name)) legend_name <- col

  p <- ggplot(df, aes(x = .x, y = .y, fill = .val_chr)) +
    geom_tile(color = "white", linewidth = gap) +
    scale_fill_manual(
      values = c("0" = color1, "1" = color2),
      na.value = na_color,
      name   = legend_name,
      drop   = FALSE
    ) +
    labs(x = xlab, y = ylab) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid       = element_blank(),
      axis.text.x      = element_text(angle = 45, hjust = 1,
                                      color = "black", size = 11),
      axis.text.y      = element_text(color = "black", size = 11),
      axis.ticks       = element_blank(),
      legend.position  = "right",
      plot.background  = element_rect(fill = "white", color = NA)
    ) +
    coord_fixed()

  p
}


# ============================================================
# plotRadialTree：三层径向树图
# ============================================================
plotRadialTree <- function(data, Objective, group, subgroup,
                           color = NULL,
                           edgealpha = 1,
                           nodealpha = 0.6,
                           nodesize = c(1, 30),
                           root_name = NULL,
                           root_color = "#888888",
                           agg_fun = sum,
                           edge_width = 0.2,
                           text_size = 2.5,
                           text_expand = 1.1,
                           size_scale = "sqrt",
                           return = "plot") {

  for (pkg in c("tidygraph", "ggraph", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("请先安装 ", pkg, ": install.packages('", pkg, "')")
  }
  library(ggplot2)

  if (!is.data.frame(data)) stop("data 必须是数据框")
  for (nm in list(Objective, group, subgroup)) {
    if (length(nm) != 1 || !is.character(nm))
      stop("Objective / group / subgroup 必须是单个列名字符串")
  }
  need <- c(Objective, group, subgroup)
  miss <- setdiff(need, colnames(data))
  if (length(miss) > 0)
    stop("数据中缺少列: ", paste(miss, collapse = ", "))

  if (is.null(root_name)) root_name <- Objective

  d <- data.frame(
    g1  = as.character(data[[group]]),
    g2  = as.character(data[[subgroup]]),
    val = data[[Objective]],
    stringsAsFactors = FALSE
  )
  d <- d[!is.na(d$g1) & !is.na(d$val), , drop = FALSE]
  if (nrow(d) == 0) stop("过滤 NA 后无有效数据")

  total_val <- agg_fun(d$val, na.rm = TRUE)
  g1_levels <- unique(d$g1)
  g1_sum    <- tapply(d$val, d$g1, agg_fun, na.rm = TRUE)

  sub_agg <- aggregate(d$val, by = list(g1 = d$g1, g2 = d$g2),
                       FUN = agg_fun, na.rm = TRUE)
  names(sub_agg)[3] <- "value"
  sub_agg$g1 <- factor(sub_agg$g1, levels = g1_levels)
  sub_agg <- sub_agg[order(sub_agg$g1, sub_agg$g2), , drop = FALSE]

  nodes_root <- data.frame(
    node_name   = "ROOT",
    node_label  = root_name,
    node_branch = "ROOT",
    abu         = as.numeric(total_val),
    group       = "ROOT",
    node_level  = "root",
    stringsAsFactors = FALSE
  )

  nodes_g1 <- data.frame(
    node_name   = paste0("G1_", g1_levels),
    node_label  = g1_levels,
    node_branch = g1_levels,
    abu         = as.numeric(g1_sum[g1_levels]),
    group       = g1_levels,
    node_level  = "group",
    stringsAsFactors = FALSE
  )

  nodes_g2 <- data.frame(
    node_name   = paste0("G2_", sub_agg$g1, "_", sub_agg$g2),
    node_label  = sub_agg$g2,
    node_branch = as.character(sub_agg$g1),
    abu         = sub_agg$value,
    group       = as.character(sub_agg$g1),
    node_level  = "subgroup",
    stringsAsFactors = FALSE
  )

  nodes <- rbind(nodes_root, nodes_g1, nodes_g2)
  rownames(nodes) <- NULL

  edges_root <- data.frame(
    from = "ROOT",
    to   = nodes_g1$node_name,
    stringsAsFactors = FALSE
  )
  edges_g1 <- data.frame(
    from = paste0("G1_", sub_agg$g1),
    to   = paste0("G2_", sub_agg$g1, "_", sub_agg$g2),
    stringsAsFactors = FALSE
  )
  edges <- rbind(edges_root, edges_g1)

  n_g1 <- length(g1_levels)
  if (is.null(color)) {
    color <- scales::hue_pal()(n_g1)
    names(color) <- g1_levels
  } else {
    if (is.null(names(color))) {
      if (length(color) < n_g1)
        stop("color 长度不足：需要至少 ", n_g1, " 个颜色")
      color <- color[seq_len(n_g1)]
      names(color) <- g1_levels
    } else {
      miss_col <- setdiff(g1_levels, names(color))
      if (length(miss_col) > 0)
        stop("color 缺少以下分组的颜色: ", paste(miss_col, collapse = ", "))
      color <- color[g1_levels]
    }
  }
  color_map <- c(color, "ROOT" = root_color)

  nodes$.size_val <- switch(size_scale,
                            sqrt = sqrt(nodes$abu),
                            log  = log1p(nodes$abu),
                            nodes$abu
  )

  graph <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = TRUE)

  p <- ggraph::ggraph(graph, layout = 'dendrogram', circular = TRUE) +
    ggraph::geom_edge_diagonal(aes(color = node1.group),
                               alpha = edgealpha,
                               linewidth = edge_width,
                               show.legend = FALSE) +
    ggraph::geom_node_point(aes(size = .size_val, color = group),
                            alpha = nodealpha) +
    ggraph::geom_node_text(aes(x = x * text_expand,
                               y = y * text_expand,
                               label = node_label,
                               angle = ggraph::node_angle(x, y),
                               filter = leaf,
                               color = group),
                           size = text_size,
                           hjust = 'outward',
                           show.legend = FALSE) +
    ggraph::geom_node_text(aes(label = node_label,
                               filter = !leaf,
                               color = group),
                           fontface = "bold",
                           size = text_size + 0.5,
                           show.legend = FALSE) +
    scale_size(range = nodesize, guide = "none") +
    scale_color_manual(values = color_map, guide = "none") +
    ggraph::scale_edge_color_manual(values = color_map, guide = "none") +
    scale_x_continuous(limits = c(-1.35, 1.35)) +
    scale_y_continuous(limits = c(-1.25, 1.25)) +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "none",
          plot.background = element_rect(fill = "white", color = "white"))

  if (identical(return, "plot")) return(p)
  if (identical(return, "data")) return(list(nodes = nodes, edges = edges, graph = graph))
  if (identical(return, "both")) return(list(plot = p, nodes = nodes, edges = edges, graph = graph))
  p
}
