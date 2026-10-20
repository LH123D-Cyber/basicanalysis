# 声明 ggplot2 非标准求值中使用的变量，避免 R CMD check 的 NOTE
utils::globalVariables(c(
  ".data", ".value", ".point_color", ".x", ".y", ".val", ".val_chr",
  ".size_val", ".bar_fill", "node1.group", "leaf", "x", "y"
))
