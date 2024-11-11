# DESeq2差异分析

# 安装DESeq2包
# BiocManager::install("DESeq2")
library(DESeq2)
# 示例数据---------------
rc <- matrix(runif(30), nrow=5)
rownames(rc) <- paste0("gene",1:5)
colnames(rc) <- c(paste0('ctrl',1:3),paste0('treat',1:3))
head(rc)

metadata <- data.frame(group = c(rep("ctrl",3),rep("treat",3)))
rownames(metadata) <- colnames(rc)
metadata$group <- factor(metadata$group, levels = c("ctrl", "treat"))# 分组需要转换为因子型
metadata

# 差异矩阵-------------
dds <- DESeqDataSetFromMatrix(countData = rc, colData = metadata, design= ~ group)
nrow(dds)

# 筛选分析-----------------
dds_filter <- dds[ rowSums(counts(dds))>1, ]# 筛选有表达值的
dds_out <- DESeq(dds_filter)
dds_res <- results(dds_out)

# 差异分析结果汇总-----------------
summary(dds_res)
table(dds_res$padj<0.05)

# 按p排序------------------------
result_sort_p <- dds_res[order(dds_res$padj),]

# 筛选差异基因-----------------------
degs <- subset(result_sort_p, padj<0.05 & (log2FoldChange > 1 | log2FoldChange < -1))

# 整理差异分析结果-----------------------
dres_diff_data <- merge(as.data.frame(dds_res),as.data.frame(counts(dds_out,normalize=TRUE)),by="row.names",sort=FALSE)
write.csv(dres_diff_data,file = "result_data/DESeq2_result.csv",row.names = FALSE,quote = FALSE)

# 主成分分析-----------------------
rld <- rlog(dds, blind = FALSE)
plotPCA(rld,intgroup=c("group"))

# 绘制热图
library("genefilter")
library("pheatmap")
topVarGene <- head(order(rowVars(assay(rld)),decreasing = TRUE),50)
mat  <- assay(rld)[ topVarGene, ]
pheatmap(mat, annotation_col=col_data)
pheatmap(mat,cluster_row=T,scale="row", annotation_col=col_data) 
res0.5 <- results(dds, contrast = c("condition","Basal","LP"),alpha=0.05)
#??һ?ֻ?ͼ??ʽ
mat  <- mat - rowMeans(mat)
pheatmap(mat, annotation_col=col_data)