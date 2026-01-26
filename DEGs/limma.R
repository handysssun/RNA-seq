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

expr <- df1

group1 <- "WT_C"
group2 <- "KO_C"

# metadata
metadata <- data.frame(rep(c(group1,group2),each = 3))
colnames(metadata) <- "group"
rownames(metadata) <- colnames(expr)
metadata$group <- factor(metadata$group)

# 创建设计矩阵,group列需要为因子类型
design <- model.matrix(~0 + metadata$group)
metadata$group <- factor(metadata$group, levels = c(group1, group2))
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
contrast_formula <- paste0(group2, " - ", group1)
contrast_matrix <- makeContrasts(contrast_formula, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
# 贝叶斯检验
fit2 <- eBayes(fit2)

# 提取差异分析结果
results <- topTable(fit2, adjust = "BH", number = Inf)

# 合并表达矩阵和symbol
results_df <- merge(expr,results,by.x = 0,by.y = 0)

# print(head(results))
write.csv(results_df, file = paste0(group1,"_VS_",group2,"_DE_results.csv"),
          quote = FALSE, row.names = TRUE)

# 去除地表达基因----------------------------------------------------------------------------------
library(limma)

# 1. 数据预处理
# 假设你的 FPKM 矩阵叫 fpkm_matrix
# 过滤掉在超过一半的样本中表达量都极低 (比如 < 1) 的基因
keep <- rowSums(fpkm_matrix > 1) >= (ncol(fpkm_matrix) / 2)
fpkm_filtered <- fpkm_matrix[keep, ]

# 2. Log2 转换 (这是 limma 处理连续数据的核心需求)
# 加 1 是为了防止对 0 取 log 导致产生无穷大
exprs_data <- log2(fpkm_filtered + 1)

# 3. 创建设计矩阵 (假设同前)
group_list <- factor(c(rep("Control", 3), rep("Treatment", 3)))
design <- model.matrix(~0 + group_list)
colnames(design) <- levels(group_list)

# 4. 差异分析
fit <- lmFit(exprs_data, design)
contrast.matrix <- makeContrasts(Treatment - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 5. 获取结果
res <- topTable(fit2, coef = 1, n = Inf, adjust.method = "fdr")



# 去除地表达基因和voom转换----------------------------------------------------------------------------------

# 为什么需要 filterByExpr()?
#在多重假设检验（FDR 校正）中，检验的次数越多，校正后的 P 值就越难显著。

#作用：剔除那些“背景噪音”基因（比如在所有样本中 Count 都是 0 或 1 的基因）。

#逻辑：它不仅看平均值，还会考虑你的分组。例如，如果一个基因只在处理组表达，而在对照组不表达，它会聪明地保留这个基因，因为它极可能是潜在的差异基因。

#为什么需要 voom 转换?
#limma 最初是为芯片数据设计的，芯片信号通常呈正态分布。而 RNA-seq 的原始 Count 数据符合负二项分布，且存在“均值-方差依赖”关系（表达量越低，噪音比例越高）。

#作用：它计算每个观察值的精度权重。

#可视化意义：当你运行 voom(..., plot=TRUE) 时，你会看到一条向下的曲线。

#横轴：平均表达量。

#纵轴：标准差。

#voom 的目的就是让这条曲线变平，消除这种依赖性，让数据能够像芯片数据一样被线性模型处理。
library(edgeR)
library(limma)

# 1. 创建 DGEList 对象
dge <- DGEList(counts = counts)

# 2. 使用 filterByExpr 过滤低表达基因
# 它会自动根据样本分组情况，剔除在大多数样本中都不表达的基因
keep <- filterByExpr(dge, design)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# 3. 计算标准化因子 (TMM 标准化)
dge <- calcNormFactors(dge)

# 4. 使用 voom 进行转换
# 这步将 count 转换为 log2-CPM，并计算权重以处理均值-方差关系
v <- voom(dge, design, plot=TRUE)

# 5. 随后即可进入 limma 拟合
fit <- lmFit(v, design)
# ... 接后续的 eBayes 步骤

