# 差异表达分析
rm(list = ls())

# 导入必要的R包
library(edgeR)

# 示例数据---------------
rc <- matrix(runif(30), nrow=5)
rownames(rc) <- paste0("gene",1:5)
colnames(rc) <- c(paste0('ctrl',1:3),paste0('treat',1:3))
head(rc)

metadata <- data.frame(group = c(rep("ctrl",3),rep("treat",3)))
rownames(metadata) <- colnames(rc)
metadata$group <- factor(metadata$group, levels = c("ctrl", "treat"))# 分组需要转换为因子型
metadata

# 创建DGEList对象------------------
y <- DGEList(counts = rc, group = metadata$group)

# 过滤低表达基因
keep <- rowSums(cpm(y) > 1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]

# 计算归一化因子
y <- calcNormFactors(y)

# 创建设计矩阵
design <- model.matrix(~ group, data = y$samples)

# 估计离散度
y <- estimateDisp(y, design)

# 拟合广义线性模型
fit <- glmFit(y, design)

# 进行似然比检验
lrt <- glmLRT(fit)

# 检查差异表达基因
top_tags <- topTags(lrt, n = Inf)

# 将差异表达基因结果保存为文件
write.csv(top_tags, file = "edgeR_result.csv",quote = FALSE,row.names = TRUE)

# 提取差异表达基因的详细信息
de_genes <- topTags(lrt, n = Inf)$table

# 转换-----------------
library(clusterProfiler)
exchange <- bitr(rownames(de_genes), fromType = 'ENSEMBL', 
                 toType = 'SYMBOL', 
                 OrgDb = org.Hs.eg.db)
de_genes$ENSEMBL <- rownames(de_genes)

merged_de_genes <- merge(de_genes, exchange, 
                         by = "ENSEMBL")

write.csv(merged_de_genes, 
          file = "differential_expression_symbol.csv", 
          quote = FALSE)

# 筛选差异表达基因
de_genes <- decideTestsDGE(lrt, p.value = 0.05)
summary(de_genes)

# 获取显著性差异表达基因列表
significant_genes <- topTags(lrt, n = Inf)

# 绘制MA图
plotMD(lrt)

# 绘制火山图
volcanoplot(lrt)

# 基因注释（假设使用org.Hs.eg.db作为人类基因注释）
library(org.Hs.eg.db)
annotated_genes <- select(org.Hs.eg.db, keys = rownames(significant_genes), columns = c("SYMBOL", "ENTREZID", "GENENAME"), keytype = "ENSEMBL")

# 功能富集分析
library(clusterProfiler)
ego <- enrichGO(gene = annotated_genes$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01)


