# clusterProfiler包富集分析
# 加载包
library(openxlsx)#读取.xlsx文件
library(dplyr)
library(stringr)#基因ID转换
library(clusterProfiler)#GO,KEGG,GSEA
library(org.Hs.eg.db)
library(org.Mm.eg.db)

# 指定物种库
GO_database <- 'org.Mm.eg.db' #GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
KEGG_database <- 'mmu' #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html
# "hsa" "mmu"
# "org.Hs.eg.db" "org.Mm.eg.db"

# 指定阈值
p_cutoff = 0.05
q_cutoff = 0.05

# 加载和整理基因列表
degs <- read.xlsx("GO/logFC1_DEGs.xlsx")
gene <- bitr(degs$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID', OrgDb = GO_database)
gene_list <- gene$ENTREZID

# GO富集分析
ego <- enrichGO(gene_list,#GO富集分析
                OrgDb = GO_database,
                keyType = "ENTREZID",#设定读取的gene ID类型
                ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                pAdjustMethod = "BH",
                pvalueCutoff = p_cutoff,#设定p值阈值
                qvalueCutoff = q_cutoff,#设定q值阈值
                readable = T)
write.csv(as.data.frame(ego),"clsuterprofiler_GO.csv",row.names = FALSE)


# KEGG富集分析
ekegg <- enrichKEGG(gene_list,
                    organism = KEGG_database,
                    pvalueCutoff = p_cutoff,
                    qvalueCutoff = q_cutoff)
write.csv(as.data.frame(ekegg), "clsuterprofiler_KEGG.csv", row.names = FALSE)





