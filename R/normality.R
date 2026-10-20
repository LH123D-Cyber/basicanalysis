#' 整体正态性检验与数据转换
#'
#' 对多个变量进行 Shapiro-Wilk 正态性检验、Q-Q 图、直方图，
#' 并自动尝试多种转换方法以寻找使数据满足正态分布的最佳方案。
#'
#' @param col 需要检验的列名（字符向量）。
#' @param QQ 逻辑值，是否绘制原始数据的 Q-Q 图，默认 `FALSE`。
#' @param His 逻辑值，是否绘制原始数据的直方图，默认 `FALSE`。
#' @param ShapiroWilk 逻辑值，是否执行 Shapiro-Wilk 检验，默认 `FALSE`。
#' @param Convert 逻辑值，是否对每个变量自动选择最佳正态化转换，默认 `FALSE`。
#' @param plot 逻辑值，是否对转换后满足正态的数据绘制 Q-Q 图和直方图，默认 `FALSE`。
#' @param data 数据框，默认使用全局环境中的 `df`。
#' @param outdir 可选，若提供目录路径则将图形保存为 PNG，否则在屏幕设备中绘制。
#'
#' @return 一个 list，包含 `ShapiroWilk`、`QQ1`、`His1`、`Convert`、`data`、`QQ2`、`His2` 等元素。
#' @export
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
  result$col <- col
  result$data_source <- deparse(substitute(data))

  # 根据变量个数自动计算多面板布局的行列数
  get_mfrow <- function(n) {
    if (n <= 2) return(c(1, n))
    return(c(ceiling(n / 2), 2))
  }

  # Yeo-Johnson 变换，支持负值
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

  # Box-Cox 变换，要求输入全部为正数
  boxcox_tf <- function(x, lambda) {
    if (abs(lambda) < 1e-10) log(x) else (x^lambda - 1) / lambda
  }

  # 图形对象构建：支持保存为 PNG 或返回可重放的绘制函数
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

  # 原始数据的 Shapiro-Wilk 检验
  if (ShapiroWilk) {
    sw_list <- list()
    for (v in col) {
      sw_list[[v]] <- tryCatch(shapiro.test(sub[[v]]),
                               error = function(e) paste("Error:", e$message))
    }
    result$ShapiroWilk <- sw_list
  }

  # 原始数据 Q-Q 图
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

  # 原始数据直方图
  if (His) {
    result$His1 <- make_plot_obj(
      draw_fun = function() {
        for (v in col) {
          hist(sub[[v]], main = paste("直方图:", v), xlab = v,
               col = "lightblue", border = "white")
        }
      },
      mfrow_par = get_mfrow(length(col)),
      name = "His1"
    )
  }

  # 自动正态化转换：遍历候选方法，以 Shapiro p 值最大者为最优
  if (Convert) {
    convert_results <- data.frame()
    converted_data <- list()
    full_data <- data

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
          # 网格搜索最优 lambda
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
          xt <- switch(m,
                       "log" = log(x), "sqrt" = sqrt(x),
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

      # 转换成功（p > 0.05）则替换完整数据集中的对应列
      if (best$p > 0.05) {
        converted_data[[v]] <- best$data
        full_data[[v]] <- best$data
      }
    }

    result$Convert <- convert_results
    result$data <- full_data
    attr(result$data, "converted_vars") <- names(converted_data)

    # 转换后满足正态的数据绘制 Q-Q 图和直方图
    if (plot && length(converted_data) > 0) {
      n <- length(converted_data)
      result$QQ2 <- make_plot_obj(
        draw_fun = function() {
          for (v in names(converted_data)) {
            qqnorm(converted_data[[v]],
                   main = paste("Q-Q（转换后）:", v),
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
                 main = paste("直方图（转换后）:", v),
                 xlab = v, col = "lightgreen", border = "white")
          }
        },
        mfrow_par = get_mfrow(n), name = "His2"
      )
    }
  }

  return(result)
}


#' 显示 OND 生成的图形
#'
#' 对 `OND()` 返回的图形对象进行重放。若提供 `data`，则用新数据重绘原始分布。
#'
#' @param result `OND()` 的返回值。
#' @param name 图形名称，如 `"QQ1"`、`"His1"`、`"QQ2"`、`"His2"`。
#' @param data 可选，若提供则用该数据框重绘（仅重绘原始分布，不重新转换）。
#'
#' @return 无返回值，仅产生图形。
#' @export
OND_show <- function(result, name, data = NULL) {
  obj <- result[[name]]
  if (is.null(obj)) {
    message(name, " 不存在（可能未生成）")
    return(invisible(NULL))
  }

  # 提供 data 时，用新数据重绘原始分布
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
        hist(sub[[v]], main = paste("直方图:", v), xlab = v,
             col = "lightblue", border = "white")
      }
    }
    return(invisible(NULL))
  }

  # 无 data 时，按录制方式重放
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
