#' 二元分类数据热图
#'
#' 将二元（0/1）分类变量绘制为热图，每个单元格代表一个 (x, y) 组合的值，
#' 0 和 1 使用不同颜色，单元格之间保留白色间隙。
#'
#' @param data 数据框。
#' @param col 二元分类列名（值应为 0 或 1）。
#' @param x x 轴列名（通常为时间）。
#' @param y y 轴列名（通常为区县或样本名）。
#' @param color1 值 0 的颜色。
#' @param color2 值 1 的颜色。
#' @param gap 单元格白边粗细。
#' @param na_color 缺失值颜色。
#' @param xlab x 轴标题，默认取 `x`。
#' @param ylab y 轴标题，默认取 `y`。
#' @param legend_name 图例标题，默认取 `col`。
#'
#' @return 一个 `ggplot` 对象。
#' @export
plotclass <- function(data, col, x, y,
                      color1 = "#cfe2f3",
                      color2 = "#e06666",
                      gap    = 1,
                      na_color = "grey90",
                      xlab = NULL, ylab = NULL,
                      legend_name = NULL) {

  requireNamespace("ggplot2", quietly = TRUE) ||
    stop("请先安装 ggplot2: install.packages('ggplot2')")

  # ---- 参数检查 ----
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

  # 只保留 0/1（以及 NA），其它值给出提示并按 NA 处理
  bad <- setdiff(unique(df$.val[!is.na(df$.val)]), c(0, 1))
  if (length(bad) > 0)
    warning("列 '", col, "' 含有非 0/1 的值: ",
            paste(bad, collapse = ", "), "（将按 NA 处理）")
  df$.val[!is.na(df$.val) & !(df$.val %in% c(0, 1))] <- NA

  # 转成字符串标签以便图例显示
  df$.val_chr <- ifelse(is.na(df$.val), NA,
                        ifelse(df$.val == 1, "1", "0"))
  df$.val_chr <- factor(df$.val_chr, levels = c("0", "1"))

  if (is.null(xlab)) xlab <- x
  if (is.null(ylab)) ylab <- y
  if (is.null(legend_name)) legend_name <- col

  # ---- 绘图 ----
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .x, y = .y, fill = .val_chr)) +
    ggplot2::geom_tile(color = "white", linewidth = gap) +
    ggplot2::scale_fill_manual(
      values = c("0" = color1, "1" = color2),
      na.value = na_color,
      name   = legend_name,
      drop   = FALSE
    ) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1,
                                               color = "black", size = 11),
      axis.text.y      = ggplot2::element_text(color = "black", size = 11),
      axis.ticks       = ggplot2::element_blank(),
      legend.position  = "right",
      plot.background  = ggplot2::element_rect(fill = "white", color = NA)
    ) +
    ggplot2::coord_fixed()

  p
}
