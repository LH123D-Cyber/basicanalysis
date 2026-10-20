#' 环形差异柱状图
#'
#' 与 `Barplot()` 类似，但支持 `circle = TRUE` 将柱状图转换为极坐标下的环形柱状图。
#'
#' @inheritParams Barplot
#' @param circle 逻辑值，`TRUE` 绘制圆形柱状图，`FALSE` 为普通柱状图。
#' @param Inter 圆形柱状图圆心留白半径，需 >= 0。
#'
#' @return 一个 `ggplot` 对象。
#' @export
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

  # ---- 检查 ----
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

  # ---- 排分组顺序 ----
  if (!is.null(x)) {
    data <- data[data$group %in% x, , drop = FALSE]
    data$group <- factor(data$group, levels = x)
  } else {
    data$group <- factor(data$group, levels = data$group)
  }
  data <- data[order(data$group), , drop = FALSE]
  group_levels <- levels(data$group)
  n_grp <- length(group_levels)

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

  # ---- 散点数据 ----
  scatter_df <- if (point) .prepare_scatter(Original, pointgroup, pointcolor_by,
                                            raw_data, mean, group_levels) else NULL

  # ---- 柱子颜色 ----
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
      fill_values <- stats::setNames(color, bar_levels)
    }
  }

  # ============================================================
  # 普通柱状图
  # ============================================================
  if (!circle) {
    if (multi_color) {
      p <- ggplot2::ggplot(data, ggplot2::aes(x = group, y = .data[[mean]],
                                              fill = .data[[bargroup]])) +
        ggplot2::geom_col(alpha = alpha, color = "black", width = 0.65) +
        ggplot2::scale_fill_manual(values = fill_values, guide = "none")
    } else {
      p <- ggplot2::ggplot(data, ggplot2::aes(x = group, y = .data[[mean]])) +
        ggplot2::geom_col(fill = color, alpha = alpha, color = "black", width = 0.65)
    }
    p <- p +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = ymin, ymax = ymax),
                             color = "black", width = 0.2, linewidth = 0.8) +
      ggplot2::labs(x = NULL, y = mean) +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::theme(axis.text = ggplot2::element_text(color = "black", size = 12),
                     panel.grid = ggplot2::element_blank())

    if (label) {
      p <- p + ggplot2::geom_text(ggplot2::aes(label = letter, y = ymax + loca),
                                  vjust = 0, size = 5)
    }
    if (point) {
      if (!is.null(scatter_df)) {
        p <- .add_scatter(p, scatter_df, pointcolor, color, group_levels, pointcolor_name)
      } else {
        p <- p + ggplot2::geom_point(ggplot2::aes(y = .data[[mean]]),
                                     size = 2, color = "white")
      }
    }
    return(p)
  }

  # ============================================================
  # 圆形柱状图
  # ============================================================
  y_top_max <- max(data$ymax + loca, na.rm = TRUE)
  y_range <- diff(range(c(data$ymin, data$ymax)))
  y_upper <- y_top_max + y_range * 0.15

  if (multi_color) {
    p <- ggplot2::ggplot(data, ggplot2::aes(x = group, y = .data[[mean]],
                                            fill = .data[[bargroup]])) +
      ggplot2::geom_col(alpha = alpha, color = "black", width = 0.65) +
      ggplot2::scale_fill_manual(values = fill_values, guide = "none")
  } else {
    p <- ggplot2::ggplot(data, ggplot2::aes(x = group, y = .data[[mean]])) +
      ggplot2::geom_col(fill = color, alpha = alpha, color = "black", width = 0.65)
  }

  p <- p +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ymin, ymax = ymax),
                           color = "black", width = 0.2, linewidth = 0.8) +
    ggplot2::coord_polar(start = 0) +
    ggplot2::ylim(-Inter, y_upper) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid  = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks  = ggplot2::element_blank(),
      axis.title  = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 12, color = "black")
    )

  if (label) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = letter, y = ymax + loca),
                                vjust = 0, size = 5)
  }

  if (point) {
    if (!is.null(scatter_df)) {
      p <- .add_scatter(p, scatter_df, pointcolor, color, group_levels, pointcolor_name)
    } else {
      p <- p + ggplot2::geom_point(ggplot2::aes(y = .data[[mean]]),
                                   size = 2, color = "white")
    }
  }

  return(p)
}
