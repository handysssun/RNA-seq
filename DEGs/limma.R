# 由于limma包能够处理连续型数据，因此它成为了处理标准化信号值进行差异分析的首选。
# limma主要用于微阵列数据分析，但同样适用于其他类型的数据，只要数据满足正态分布或近似正态分布。
rm(list = ls())
setwd("E:/R_Proj/20240604_GSE54584_PC3_siGGCT/")

# 差异表达分析
library(limma)
library(Biobase)

# 三个处理组时分析如下
# 读取表达矩阵
# 后续分析需要矩阵格式，此处把基因名用作行名
expr <- read.csv("expr_sc.csv", header = TRUE,row.names = 1)
# 表达矩阵的列名作为样本信息的行名
metadata  <- read.csv("sample_info.csv", header = TRUE,row.names = 1)

# 出现assayData和phenoData名字不一致时进行
"""
# 对样本名称进行排序
expr <- expr[, order(expr_sample_names)]
metadata <- metadata[order(metadata_sample_names), ]

# 确认排序后的样本名称一致
expr_sample_names <- colnames(expr)
metadata_sample_names <- metadata$sample

print(expr_sample_names)
print(metadata_sample_names)
"""

# 创建设计矩阵
design <- model.matrix(~0 + metadata $group)
colnames(design) <- c("control","treat1","treat2")
print(design)

# 创建eset
# 创建 AnnotatedDataFrame 对象
pheno_data <- new("AnnotatedDataFrame", data = metadata)
# 创建 ExpressionSet 对象
eset <- ExpressionSet(assayData = as.matrix(expr), phenoData = pheno_data)

# 数据标准化（可选）
eset <- normalizeBetweenArrays(eset)

# 拟合线性模型
fit <- lmFit(eset, design)
# 创建比较矩阵
contrast_matrix <- makeContrasts(
  Treatment1_vs_Control = treat1 - control,
  Treatment2_vs_Control = treat2 - control,
  Treatment2_vs_Treatment1 = treat2 - treat1,
  levels = design
)
fit2 <- contrasts.fit(fit, contrast_matrix)
# 贝叶斯检验
fit2 <- eBayes(fit2)


# 若出现方差为零的情况可以对数据进行过滤（可选）
'''
non_zero_var_genes <- apply(eset, 1, var) != 0
filtered_data <- eset[non_zero_var_genes, ]
# 拟合线性模型
fit <- lmFit(filtered_data, design)
# 创建比较矩阵
contrast_matrix <- makeContrasts(
  Treatment1_vs_Control = treat1 - control,
  Treatment2_vs_Control = treat2 - control,
  Treatment2_vs_Treatment1 = treat2 - treat1,
  levels = design
)
fit2 <- contrasts.fit(fit, contrast_matrix)
# 贝叶斯检验
fit2 <- eBayes(fit2)
'''

# 提取差异表达基因
# 提取 Treatment1 vs Control 的结果
results_Treatment1_vs_Control <- topTable(fit2, coef = "Treatment1_vs_Control", adjust = "BH", number = Inf)
write.csv(results_Treatment1_vs_Control, file = "results_Treatment1_vs_Control.csv", quote = FALSE, row.names = TRUE)

# 提取 Treatment2 vs Control 的结果
results_Treatment2_vs_Control <- topTable(fit2, coef = "Treatment2_vs_Control", adjust = "BH", number = Inf)
write.csv(results_Treatment2_vs_Control, file = "results_Treatment2_vs_Control.csv", quote = FALSE, row.names = TRUE)

# 提取 Treatment2 vs Treatment1 的结果
results_Treatment2_vs_Treatment1 <- topTable(fit2, coef = "Treatment2_vs_Treatment1", adjust = "BH", number = Inf)
write.csv(results_Treatment2_vs_Treatment1, file = "results_Treatment2_vs_Treatment1.csv", quote = FALSE, row.names = TRUE)


# 两个处理组时分析如下
# 读取表达矩阵
# 后续分析需要矩阵格式，此处把基因名用作行名
expr <- read.csv("expr_sc.csv", header = TRUE,row.names = 1)
# 表达矩阵的列名作为样本信息的行名
metadata  <- read.csv("sample_info.csv", header = TRUE,row.names = 1)

# metadata
metadata <- data.frame(rep(c("Control","Treatment"),each = 3))
colnames(metadata) <- "group"
rownames(metadata) <- colnames(c1)
metadata$group <- factor(metadata$group)

# 创建设计矩阵,group列需要为因子类型
design <- model.matrix(~0 + metadata$group)
metadata$group <- factor(metadata$group, levels = c("Control", "Treatment"))
colnames(design) <- levels(metadata$group)
print(design)

# 创建eSet,需要Biobase包，assayData需要是矩阵格式的
pheno_data <- new("AnnotatedDataFrame", data = metadata)
eset <- ExpressionSet(assayData = as.matrix(expr), phenoData = pheno_data)

# 数据标准化（可选）
eset <- normalizeBetweenArrays(eset)

# 拟合线性模型
fit <- lmFit(exprs(eset), design)

# 创建比较矩阵
contrast_matrix <- makeContrasts(Treatment - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
# 贝叶斯检验
fit2 <- eBayes(fit2)

# 提取差异分析结果
results <- topTable(fit2, adjust = "BH", number = Inf)
# print(head(results))
write.csv(results, file = "differential_expression_results.csv", quote = FALSE, row.names = TRUE)




