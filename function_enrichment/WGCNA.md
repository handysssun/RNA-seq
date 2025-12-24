---
title: "Untitled Analysis"
author: "handy"
date: "`r Sys.Date()`"
output:
  rmdformats::readthedown:
    self_contained: true
    thumbnails: true
    lightbox: true
    gallery: true
    css: "D:/Rlibrary/rmdformats/style.css"                 
  prettydoc::html_pretty:   
    theme: hpstr            #prettydoc包的hpstr风格
    toc: yes
  html_document:
    theme: flatly          # 现代主题（可选：cosmo，flatly, yeti, readable）
    highlight: tango       # 代码高亮
    toc: true              # 自动生成目录
    toc_float: true        # 浮动目录（侧边栏）
    number_sections: true  # 章节编号
    code_folding: show     # 代码折叠（"show" 或 "hide"）
  word_document: default   # 同时支持 Word
  # pdf_document: default    # 同时支持 PDF 输出（需 LaTeX）
editor_options:
  markdown:
    wrap: 72
---

```{r setup, include=FALSE,echo=FALSE}
# =============================================================
# 🛠️ 全局设置 —— 所有 chunk 的默认行为（不显示自身）
# =============================================================
knitr::opts_chunk$set(
  echo = TRUE,              # 显示代码（教学/可复现性）
  eval = TRUE,              # 执行代码（设 FALSE 可跳过整块）
  include = TRUE,           # 包含输出（设 FALSE 静默执行）
  message = FALSE,          # 隐藏包加载消息
  warning = TRUE,          # 显示警告
  error = FALSE,            # 出错即停（调试时可临时设 TRUE 继续渲染）
  fig.align = "center",     # 图形居中
  fig.width = 7,            # 默认图形宽（英寸）
  fig.height = 5,           # 默认图形高
  dpi = 120,                # 图形分辨率（PDF 建议 300）
  cache = TRUE,            # 缓存结果加速（大数据可设 TRUE；注意：修改代码后需 clean_cache）
  autodep = TRUE            # 自动追踪 chunk 依赖（配合 cache 使用）
)

# =============================================================
# 📦 加载常用包 —— 建议显式加载，避免隐式依赖
# =============================================================
library(tidyverse)    # 数据处理/绘图核心
library(dplyr)
library(openxlsx)
library(stringr)
library(knitr)        # 表格美化
library(kableExtra)   # 高级表格
library(patchwork)    # 图形组合
# library(limma)      # 微阵列（按需取消注释）
# library(DESeq2)     # RNA-seq（按需取消注释）

# =============================================================
# ⚙️ 其他全局设置
# =============================================================
options(scipen = 999)  # 禁用科学计数法（如 1e-04 → 0.0001）
options(digits = 4)    # 数值显示4位小数
theme_set(theme_minimal(base_size = 12))  # 全局 ggplot 主题
```

# Overall design

# Analysis

WGCNA (加权基因共表达网络分析) 是一种用于从基因表达数据中寻找共表达基因模块并进行生物学分析的方法。在 R 中进行 WGCNA 分析的全流程通常包括以下步骤：

# 安装和加载必要的 R 包

```{R}
# install.packages("WGCNA")
library(WGCNA)
```

# 数据准备

WGCNA 需要一个基因表达数据矩阵，通常是一个基因（行）与样本（列）对应的矩阵。

## 数据导入

```{R}
load("GSE46408_RAW/expr.Rdata")
expr <- t(expr)
# 检查数据
expr[1:5,1:5]

```

# 数据预处理

在进行 WGCNA 分析之前，通常需要对数据进行一定的预处理，如去除缺失值、对数据进行标准化等。

## 检查缺失值

```{R}
gsg <- goodSamplesGenes(expr, verbose = 3)
if (!gsg$allOK) {
  # 去除不合格的基因或样本
  expr <- expr[gsg$goodSamples, gsg$goodGenes]
}
# 检查数据
expr[1:5,1:5]
```

## 对数据进行标准化（可选）

有时我们需要对数据进行标准化，以确保每个基因的表达值的均值为 0，标准差为 1。

```{R}
expr <- t(scale(t(expr)))
```

# 构建基因共表达网络

## 选择合适的软阈值 (power)

WGCNA 的第一步是选择一个适当的软阈值，以确保网络是无尺度的。这是通过计算不同阈值下网络的拓扑特性来完成的。
选择一个 `R^2` 接近 0.9 的阈值作为合适的 `power`。

```{R}
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(expr, powerVector = powers, verbose = 5)
# 绘制软阈值与网络拓扑的关系图
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     type = "n", xlab = "Soft threshold power", 
     ylab = "Scale Free Topology Model Fit, signed R^2")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, col="red")
```

## 构建共表达网络

```{R}
power <- 14  # 根据前面的选择R^2最接近0.9时，软阈值选择14
net <- blockwiseModules(expr, power = power, TOMType = "unsigned", 
                        minModuleSize = 30, reassignThreshold = 0, 
                        mergeCutHeight = 0.25, numericLabels = TRUE, 
                        pamRespectsDendro = FALSE, verbose = 3)
```

# 模块分析

## 获取模块信息

```{R}
moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)

```


## 绘制树状图

```{r}
for (b in 1:length(net$dendrograms)) {
  plotDendroAndColors(
    net$dendrograms[[b]],
    moduleColors[ net$blockGenes[[b]] ],
    "Module colors",
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05,
    main = paste("Gene dendrogram and module colors - block", b)
  )
}

```


## 提取模块基因

```{R}
moduleGenes <- split(colnames(expr), moduleColors)
# moduleGenes 是一个 list，每个元素是一个模块的基因向量
genes_blue <- moduleGenes[["blue"]] #如果你只想要某个模块（比如 blue）

```

## 寻找目标基因在哪个模块

```{r}
"DHX9" %in% moduleGenes          # TRUE/FALSE 是否存在

which(moduleGenes == "DHX9")     # 返回位置(可能多个)

match("DHX9", moduleGenes)       # 返回第一个位置(找不到是 NA)

```

# 与临床特征相关性分析

可以将模块与临床特征进行相关性分析，如性别、年龄、病理数据等。

## 准备临床特征数据

```{R}
clinicalData <- read.csv("clinical_data.csv")
# 保证临床数据行名 = 样本名
clinicalData <- clinicalData[match(rownames(MEs), rownames(clinicalData)), , drop = FALSE]
stopifnot(all(rownames(clinicalData) == rownames(MEs)))

# 快速把非数值列转成因子再转数值（仅适合二分类且你确认顺序）
clinicalData_num <- clinicalData
for (cn in colnames(clinicalData_num)) {
  if (!is.numeric(clinicalData_num[[cn]])) {
    clinicalData_num[[cn]] <- as.numeric(as.factor(clinicalData_num[[cn]]))
  }
}

```

## 计算模块与临床特征的相关性

```{R}
moduleTraitCor <- cor(MEs, clinicalData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
```

## 绘制模块与临床特征的相关性图

```{R}
textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")")
dim(textMatrix) <- dim(moduleTraitCor)

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(clinicalData_num),
  yLabels = colnames(MEs),
  ySymbols = colnames(MEs),
  textMatrix = textMatrix,
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  setStdMargins = FALSE,
  cex.text = 0.6,
  zlim = c(-1, 1),
  main = "Module-trait relationships"
)

```

# 模块功能注释

使用 GO 或 KEGG 富集分析来注释模块的功能。你可以使用像 `clusterProfiler` 或 `GOstats` 等 R 包进行富集分析。

## GO 富集分析

```{R}
library(clusterProfiler)
library(org.Hs.eg.db)

# 以 blue 模块做 GO 富集
module <- "blue"
genes <- moduleGenes[[module]]   # 这个 genes 应该是 SYMBOL（如果你列名就是SYMBOL）

ego <- enrichGO(
  gene = genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

dotplot(ego)
```

# 模块保存和导出

最后，可以将分析结果保存或导出，供进一步的分析和展示。

```{R}
# 导出成“基因-模块”两列表
moduleTable <- data.frame(
  gene = colnames(datExpr),
  module = moduleColors,
  stringsAsFactors = FALSE
)
write.csv(moduleTable, "module_gene_table.csv", row.names = FALSE)
# 每个模块单独一个文件
for (m in names(moduleGenes)) {
  write.table(moduleGenes[[m]],
              file = paste0("module_", m, "_genes.txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}

```

# 结果解读

根据模块与临床特征的相关性分析，确定与表型显著相关的模块，并结合富集分析结果来解释这些模块的生物学意义。


---

这个过程涵盖了 WGCNA 分析的基本流程，但也可以根据具体的数据和需求进行适当调整。每一步的选择和参数调整都需要根据你的数据进行优化。

