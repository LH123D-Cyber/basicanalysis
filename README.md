# basicanalysis

<!-- badges: start -->
[![R-CMD-check](https://github.com/LH123D-Cyber/basicanalysis/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/LH123D-Cyber/basicanalysis/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`basicanalysis` 是一个面向环境监测数据的统计分析与可视化工具包，提供从**数据正态性检验**、**差异分析方法决策**、**组间差异分析与事后检验**，到**多种出版级图形**（柱状图、环形柱状图、热图、三层环形树状图）的一站式流程。

## 功能模块

| 模块 | 函数 | 说明 |
|---|---|---|
| 正态性检验与转换 | `OND()`、`OND_show()` | Shapiro-Wilk 检验、Q-Q 图、直方图，自动搜索最优正态化转换（Yeo-Johnson、Box-Cox、log 等） |
| 差异分析方法决策 | `DDT()` | 根据各组正态性和方差齐性，自动推荐差异检验与事后检验方法 |
| 差异分析与事后比较 | `DiffAnalysis()` | 支持 ANOVA / Welch ANOVA / Kruskal-Wallis / t / Welch t / Wilcoxon，事后检验支持 Tukey HSD / Games-Howell / Dunn 等，自动生成差异字母（CLD） |
| 普通差异柱状图 | `Barplot()` | 均值 + 误差棒 + 差异字母，可叠加原始观测散点，支持散点独立配色 |
| 环形差异柱状图 | `CircleBar()` | 在 `Barplot()` 基础上支持极坐标环形显示 |
| 二元分类热图 | `plotclass()` | 0/1 分类变量的热图，单元格白色间隙 |
| 三层环形树状图 | `plotRadialTree()` | 自动整理"总和 → 分组汇总 → 分组×子分组"三层结构并绘制环形树状图 |

## 安装

从 GitHub 安装开发版：

```r
# 方法一：使用 pak（推荐）
# install.packages("pak")
pak::pak("LH123D-Cyber/basicanalysis")

# 方法二：使用 devtools
# install.packages("devtools")
devtools::install_github("LH123D-Cyber/basicanalysis")
```

部分功能依赖可选包，按需安装：

```r
install.packages(c("car", "rstatix", "dunn.test", "FSA"))
```

- `car`：Levene 方差齐性检验（未安装时使用内置简化版）
- `rstatix`：Games-Howell 事后检验（未安装时回退到成对 Welch t 检验）
- `dunn.test` 或 `FSA`：Dunn 事后检验（未安装时回退到成对 Wilcoxon 检验）

## 快速开始

```r
library(basicanalysis)
```

假设数据框 `dt` 包含以下列：`class`（区县）、`group`（年份）、`DO`（溶解氧）、`RE`（0/1 二分类）等。

### 1. 正态性检验与自动转换

```r
res <- OND(col = c("pH", "DO", "EC"),
           QQ = TRUE, His = TRUE,
           ShapiroWilk = TRUE,
           Convert = TRUE, plot = TRUE,
           data = dt)

# 原始数据的 Shapiro-Wilk 检验结果
res$ShapiroWilk

# 各变量的转换方法与转换后检验结果
res$Convert
#   variable      method lambda      W p.value normal_after
# 1       pH Yeo-Johnson   0.60 0.9823  0.7124         Yes
# 2       DO         log     NA 0.9754  0.5310         Yes
# 3       EC    Box-Cox   -0.40 0.9611  0.2903         Yes

# 转换后的完整数据框（成功的列已替换）
head(res$data)

# 图形重放
OND_show(res, "QQ1")
OND_show(res, "His1")
OND_show(res, "QQ2")
OND_show(res, "His2")
```

### 2. 差异分析方法决策

```r
# 判断 DO 在不同区县之间应使用哪种差异检验
ddt <- DDT(data = dt, Objective = "DO", group = "class")
```

输出会自动检验每组正态性和方差齐性，并给出推荐方法，例如：

```
【3. 决策建议】
DO：数据不满足正态分布，且方差不齐（Kruskal-Wallis 对两者均不敏感）。
建议使用 Kruskal-Wallis 秩和检验 进行差异分析，
并使用 Dunn 检验（BH 校正） 进行事后检验。
```

### 3. 差异分析与事后比较

```r
diff_res <- DiffAnalysis(data = dt,
                         Objective = "DO",
                         group = "class",
                         Diff = "Kruskal-Wallis",
                         Posthoc = "Dunn",
                         Adjust = "BH")

# 绘图数据框（含均值、sd、se、差异字母）
diff_res$summary_data
#        group  n     mean       sd       se letter
# 1  Changchun  6 10.56500 1.024389 0.418204      a
# 2      Jilin  6 11.79500 1.893830 0.773115      a
# 3     Siping  6 10.28000 3.134263 1.279506      a
# 4   Liaoyuan  6 10.29000 1.539708 0.628583      a
# 5    Tonghua  6 10.78000 3.239126 1.322327      a
# 6    Baishan  6 10.30000 1.977624 0.807363      a
# 7   Songyuan  6 13.55000 7.699870 3.143462      a
# 8   Baicheng  6 10.35000 2.159976 0.881813      a

# 两两比较 p 值矩阵
diff_res$p_mat
```

### 4. 绘制差异柱状图

```r
library(ggplot2)

p <- Barplot(data       = diff_res$summary_data,
             mean       = "mean",
             Whiskers   = "se",
             point      = TRUE,
             x          = c("Changchun", "Jilin", "Siping", "Liaoyuan",
                            "Tonghua", "Baishan", "Songyuan", "Baicheng"),
             color      = c("#889f4e", "#f8d042", "#5a97d0", "#c1d7ef",
                            "#5dce9c", "#b58ade", "#4fd1c5", "#f08a5d"),
             alpha      = 0.7,
             label      = TRUE,
             loca       = 0.5,
             Original   = dt$DO,         # 原始观测
             pointgroup = dt$class,      # 散点位置：区县
             pointcolor_by = dt$group,   # 散点颜色：年份（与位置解耦）
             pointcolor = c("2015" = "#36663e", "2016" = "#6a9b3a",
                            "2017" = "#b7c94f", "2018" = "#f2d656",
                            "2019" = "#d97c3a", "2020" = "#a52c2c"),
             pointcolor_name = "年份")
p
ggsave("Barplot_DO.png", p, width = 9, height = 5.5, dpi = 300)
```

### 5. 环形差异柱状图

```r
p_circle <- CircleBar(data       = diff_res$summary_data,
                      mean       = "mean",
                      Whiskers   = "se",
                      point      = TRUE,
                      x          = c("Changchun", "Jilin", "Siping", "Liaoyuan",
                                     "Tonghua", "Baishan", "Songyuan", "Baicheng"),
                      color      = c("#889f4e", "#f8d042", "#5a97d0", "#c1d7ef",
                                     "#5dce9c", "#b58ade", "#4fd1c5", "#f08a5d"),
                      alpha      = 0.7,
                      label      = TRUE,
                      loca       = 0.5,
                      Original   = dt$DO,
                      pointgroup = dt$class,
                      pointcolor_by = dt$group,
                      pointcolor = c("2015" = "#36663e", "2016" = "#6a9b3a",
                                     "2017" = "#b7c94f", "2018" = "#f2d656",
                                     "2019" = "#d97c3a", "2020" = "#a52c2c"),
                      pointcolor_name = "年份",
                      circle     = TRUE,
                      Inter      = 5)
p_circle
ggsave("CircleBar_DO.png", p_circle, width = 8, height = 8, dpi = 300)
```

### 6. 二元分类热图

```r
# 行 = 区县，列 = 年份，值 = RE（0/1）
p_heat <- plotclass(data   = dt,
                    col    = "RE",
                    x      = "group",
                    y      = "class",
                    color1 = "#cfe2f3",   # 0 的颜色
                    color2 = "#e06666",   # 1 的颜色
                    gap    = 1,
                    xlab   = "年份",
                    ylab   = "区县",
                    legend_name = "RE")
p_heat
ggsave("plotclass_RE.png", p_heat, width = 8, height = 5, dpi = 300)
```

### 7. 三层环形树状图

```r
# 根节点 = DO 总和
# 第二层 = 各区县 DO 汇总
# 第三层 = 各区县 × 各年份的具体值
p_tree <- plotRadialTree(
  data      = dt,
  Objective = "DO",
  group     = "class",
  subgroup  = "group",
  color     = c("Changchun" = "#889f4e", "Jilin"     = "#f8d042",
                "Siping"    = "#5a97d0", "Liaoyuan"  = "#c1d7ef",
                "Tonghua"   = "#5dce9c", "Baishan"   = "#b58ade",
                "Songyuan"  = "#4fd1c5", "Baicheng"  = "#f08a5d"),
  edgealpha = 1,
  nodealpha = 0.7,
  nodesize  = c(1, 30),
  root_name = "DO"
)
p_tree
ggsave("plotRadialTree_DO.png", p_tree, width = 10, height = 10, dpi = 300)
```

`plotRadialTree()` 还支持返回节点/边数据：

```r
res <- plotRadialTree(dt, "DO", "class", "group", return = "both")
res$plot      # 图
res$nodes     # 1（根）+ 8（区县）+ 48（区县×年份）= 57 行
res$edges     # 56 条边
res$graph     # tbl_graph 对象，可继续用 ggraph 定制
```

## 完整工作流

```r
library(basicanalysis)

# 1. 先检验正态性（若需要可自动转换）
norm_res <- OND(col = c("DO", "EC", "TN"),
                ShapiroWilk = TRUE, Convert = TRUE,
                data = dt)
norm_res$Convert

# 2. 用决策函数判断该用哪种差异检验
ddt <- DDT(data = dt, Objective = "DO", group = "class")

# 3. 按推荐方法执行差异分析
diff_res <- DiffAnalysis(dt, "DO", "class",
                         Diff = ddt$Decision$test,
                         Posthoc = ddt$Decision$posthoc)

# 4. 绘制柱状图或环形柱状图
Barplot(diff_res$summary_data,
        Whiskers = "se", point = TRUE,
        Original = dt$DO, pointgroup = dt$class)
```

## 作者

- **Hao Li**（维护者）— Plague and Brucellosis Prevention and Control Base,
  Chinese Center for Disease Control and Prevention & Chinese Academy of Preventive Medicine
- **Cheng Ju**、**Xiaolei Zhou**、**Chi Song**、**Xiaochen Zhang**
  — 同上

## 许可证

MIT © China University of Geosciences, Wuhan;
Plague and Brucellosis Prevention and Control Base,
Chinese Center for Disease Control and Prevention & Chinese Academy of Preventive Medicine
