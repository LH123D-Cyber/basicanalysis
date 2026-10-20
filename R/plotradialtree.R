#' 三层环形树状图
#'
#' 自动从原始数据整理出三层节点结构（根节点 = 总和，第二层 = 分组汇总，
#' 第三层 = 分组 × 子分组的具体值），并用 `ggraph` 绘制圆形树状图。
#'
#' @param data 数据框。
#' @param Objective 指标列名，根节点为其总和。
#' @param group 第二层分组列名。
#' @param subgroup 第三层分组列名。
#' @param color 第二层各组的颜色向量；`NULL` 时自动生成。
#' @param edgealpha 边透明度。
#' @param nodealpha 节点透明度。
#' @param nodesize 节点大小范围。
#' @param root_name 根节点标签，默认取 `Objective`。
#' @param root_color 根节点和根边颜色。
#' @param agg_fun 聚合函数，默认 `sum`。
#' @param edge_width 边线宽。
#' @param text_size 标签字号。
#' @param text_expand 叶子标签向外的距离倍数。
#' @param size_scale 节点大小映射方式，`"sqrt"`、`"log"` 或 `"linear"`。
#' @param return 返回类型，`"plot"`、`"data"` 或 `"both"`。
#'
#' @return 根据 `return` 返回 `ggplot` 对象或节点/边数据。
#' @export
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

  # ---- 依赖 ----
  for (pkg in c("tidygraph", "ggraph", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("请先安装 ", pkg, ": install.packages('", pkg, "')")
  }

  # ---- 参数检查 ----
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

  # ---- 准备数据 ----
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
  g1_sum <- tapply(d$val, d$g1, agg_fun, na.rm = TRUE)

  # 第三层汇总
  sub_agg <- stats::aggregate(d$val, by = list(g1 = d$g1, g2 = d$g2),
                              FUN = agg_fun, na.rm = TRUE)
  names(sub_agg)[3] <- "value"
  sub_agg$g1 <- factor(sub_agg$g1, levels = g1_levels)
  sub_agg <- sub_agg[order(sub_agg$g1, sub_agg$g2), , drop = FALSE]

  # ---- 节点表 ----
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

  # ---- 边表 ----
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

  # ---- 颜色 ----
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

  # ---- 大小缩放 ----
  nodes$.size_val <- switch(size_scale,
                            sqrt = sqrt(nodes$abu),
                            log  = log1p(nodes$abu),
                            nodes$abu
  )

  graph <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = TRUE)

  # ---- 绘图 ----
  p <- ggraph::ggraph(graph, layout = 'dendrogram', circular = TRUE) +
    ggraph::geom_edge_diagonal(ggplot2::aes(color = node1.group),
                               alpha = edgealpha,
                               linewidth = edge_width,
                               show.legend = FALSE) +
    ggraph::geom_node_point(ggplot2::aes(size = .size_val, color = group),
                            alpha = nodealpha) +
    ggraph::geom_node_text(ggplot2::aes(x = x * text_expand,
                                        y = y * text_expand,
                                        label = node_label,
                                        angle = ggraph::node_angle(x, y),
                                        filter = leaf,
                                        color = group),
                           size = text_size,
                           hjust = 'outward',
                           show.legend = FALSE) +
    ggraph::geom_node_text(ggplot2::aes(label = node_label,
                                        filter = !leaf,
                                        color = group),
                           fontface = "bold",
                           size = text_size + 0.5,
                           show.legend = FALSE) +
    ggplot2::scale_size(range = nodesize, guide = "none") +
    ggplot2::scale_color_manual(values = color_map, guide = "none") +
    ggraph::scale_edge_color_manual(values = color_map, guide = "none") +
    ggplot2::scale_x_continuous(limits = c(-1.35, 1.35)) +
    ggplot2::scale_y_continuous(limits = c(-1.25, 1.25)) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none",
                   plot.background = ggplot2::element_rect(
                     fill = "white", color = "white"))

  # ---- 返回 ----
  if (identical(return, "plot")) return(p)
  if (identical(return, "data")) return(list(nodes = nodes, edges = edges, graph = graph))
  if (identical(return, "both")) return(list(plot = p, nodes = nodes, edges = edges, graph = graph))
  p
}
