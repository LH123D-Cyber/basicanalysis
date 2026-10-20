# ============================================================
# 辅助函数 1：整理散点数据（支持独立的颜色分组）
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

  # 位置分组一致性诊断
  src_levels <- unique(df$group)
  common <- intersect(src_levels, target_levels)
  if (length(common) == 0) {
    stop("散点位置分组与 data$group 的水平完全不一致！\n",
         "  · data$group 水平 : ", paste(target_levels, collapse = ", "), "\n",
         "  · 传入分组水平    : ", paste(src_levels, collapse = ", "), "\n",
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
# 辅助函数 2：添加散点图层（优先按 .point_color 上色）
# ============================================================
.add_scatter <- function(p, scatter_df, pointcolor, color, group_levels,
                         pointcolor_name = NULL) {
  if (is.null(scatter_df) || nrow(scatter_df) == 0) return(p)

  has_point_color <- ".point_color" %in% colnames(scatter_df)

  # ---- 场景 A：散点颜色按 pointcolor_by 分组 ----
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
      p <- p + ggplot2::geom_jitter(data = scatter_df,
                                    ggplot2::aes(x = group, y = .value, color = .point_color),
                                    width = 0.15, size = 1.5, alpha = 0.5,
                                    inherit.aes = FALSE) +
        ggplot2::labs(color = legend_name)
    } else {
      p <- p + ggplot2::geom_jitter(data = scatter_df,
                                    ggplot2::aes(x = group, y = .value, color = .point_color),
                                    width = 0.15, size = 1.5, alpha = 0.5,
                                    inherit.aes = FALSE) +
        ggplot2::scale_color_manual(values = values, name = legend_name)
    }
    return(p)
  }

  # ---- 场景 B：无 pointcolor_by，沿用按组着色的逻辑 ----
  use_group_color <- FALSE
  color_values <- NULL

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
    p + ggplot2::geom_jitter(data = scatter_df,
                             ggplot2::aes(x = group, y = .value, color = group),
                             width = 0.15, size = 1.5, alpha = 0.5,
                             inherit.aes = FALSE) +
      ggplot2::scale_color_manual(values = color_values, guide = "none")
  } else {
    p + ggplot2::geom_jitter(data = scatter_df,
                             ggplot2::aes(x = group, y = .value),
                             width = 0.15, size = 1.5, alpha = 0.5,
                             color = color_values, inherit.aes = FALSE)
  }
}
