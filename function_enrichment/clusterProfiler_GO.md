---
title: "R Notebook"
output: html_notebook
---

# clusterProfiler包富集和可视化

```{r }
# 加载包
library(openxlsx)#读取.xlsx文件
library(dplyr)
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例
library(org.Hs.eg.db)
library(fgsea)
library(gprofiler2)
library(GOstats)
library(ggvolcano)
library(ggrepel)
library(RColorBrewer)
library(tidyverse)
```

指定富集分析的物种库
```{R}
GO_database <- 'org.Hs.eg.db' #GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
KEGG_database <- 'hsa' #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html
```

载入富集分析的基因并转换
```{r}
degs <- read.xlsx("GO/logFC1_DEGs.xlsx")
keytypes(org.Hs.eg.db)# ID类型

library(clusterProfiler)
library(org.Hs.eg.db)

gene <- bitr(degs$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID', OrgDb = GO_database)
gene_list<- gene$ENTREZID

```

进行富集(一步全出)
```{r}
clusterprofiler_GO_enrichment <- function(input_file,
                                          p_cutoff = 0.05,
                                          q_cutoff = 0.05,
                                          GO_database = "org.Hs.eg.db",
                                          KEGG_database = "hsa"){
  library(clusterProfiler)
  GO <- enrichGO( input_file,#GO富集分析
                  OrgDb = GO_database,
                  keyType = "ENTREZID",#设定读取的gene ID类型
                  ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                  pAdjustMethod = "BH",
                  pvalueCutoff = p_cutoff,#设定p值阈值
                  qvalueCutoff = q_cutoff,#设定q值阈值
                  readable = T)
  return(GO)
  write.csv(as.data.frame(GO),"clsuterprofiler_GO.csv",row.names = FALSE)
  KEGG <- enrichKEGG(gene_list,
                     organism = KEGG_database,
                     pvalueCutoff = p_cutoff,
                     qvalueCutoff = q_cutoff)
  return(KEGG)
  write.csv(as.data.frame(KEGG), "clsuterprofiler_KEGG.csv", row.names = FALSE)
}
clusterprofiler_GO_enrichment(input_file = gene_list)
```

分别富集
```{r}
library(clusterProfiler)
ego <- enrichGO( gene_list,#GO富集分析
			  OrgDb = GO_database,
			  keyType = "ENTREZID",#设定读取的gene ID类型
			  ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
			  pAdjustMethod = "BH",
			  pvalueCutoff = 0.05,#设定p值阈值
			  qvalueCutoff = 0.05,#设定q值阈值
			  readable = T)
			  
ekegg <- enrichKEGG(gene = gene$ENTREZID, #基因列表文件中的基因名称
					keyType = 'kegg', #KEGG 富集
					organism = "hsa",#例如，hsa 代表人类，其它物种更改这行即可
					pAdjustMethod = "BH",
					pvalueCutoff = 0.05, #指定 p 值阈值（可指定 1 以输出全部）
					qvalueCutoff = 0.05) #指定 q 值阈值（可指定 1 以输出全部）
	
```

cluster Profiler可视化
简单绘制条形图和气泡图
```{r}
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# 绘制条形图
barplot(ego, showCategory=10, title="GO Enrichment Analysis")

# 绘制气泡图
dotplot(ego, showCategory=10, title="GO Enrichment Analysis")
```

分别绘制不同ont的图
```{r}
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# 按 GO 分类分成不同的数据框
ego_bp <- subset(ego, ONTOLOGY == "BP")
ego_cc <- subset(ego, ONTOLOGY == "CC")
ego_mf <- subset(ego, ONTOLOGY == "MF")
go_bp <- go_bp %>% arrange(pvalue)
go_cc <- go_cc %>% arrange(pvalue)
go_mf <- go_mf %>% arrange(pvalue)

# 分别绘制 BP、CC 和 MF 分类的条形图和气泡图
barplot(ego_bp, showCategory=10, title="GO Enrichment Analysis: Biological Process (BP)")
dotplot(ego_bp, showCategory=10, title="GO Enrichment Analysis: Biological Process (BP)")

barplot(ego_cc, showCategory=10, title="GO Enrichment Analysis: Cellular Component (CC)")
dotplot(ego_cc, showCategory=10, title="GO Enrichment Analysis: Cellular Component (CC)")

barplot(ego_mf, showCategory=10, title="GO Enrichment Analysis: Molecular Function (MF)")
dotplot(ego_mf, showCategory=10, title="GO Enrichment Analysis: Molecular Function (MF)")

# 分面绘制
dotplot(ego, showCategory= c(catebp,catecc,catemf),split = "ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale="free")+
  ggtitle("Dotplot for GO Enrichment")+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40))  # 设置标签的最大宽度为 15
```

按照筛选条目进行绘制dotplot
```{r}
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
catebp <- c("regulation of cell cycle phase transition",
            "DNA replication",
            "cell cycle checkpoint signaling",
            "intrinsic apoptotic signaling pathway",
            "regulation of protein ubiquitination",
            "cellular response to external stimulus",
            "mitochondrial transport",
            "cell growth",
            "cellular response to reactive oxygen species",
            "membrane assembly")
dotplot(ego, showCategory=catebp) + ggtitle("Dotplot for GO:BP Enrichment")
ggsave("单独的功能富集/HCC/BP.pdf", width = 6,height = 10)

catecc <- c("focal adhesion",
            "cell-substrate junction",
            "chromosomal region",
            "nuclear envelope",
            "mitochondrial matrix",
            "nuclear chromosome",
            "ubiquitin ligase complex",
            "outer membrane",
            "organelle outer membrane",
            "mitochondrial inner membrane")
dotplot(ego, showCategory= catecc) + ggtitle("Dotplot for GO:CC Enrichment")
ggsave("单独的功能富集/HCC/CC.pdf", width = 6,height = 10)

catemf <- c("cadherin binding",
            "DNA-binding transcription factor binding",
            "ubiquitin-like protein ligase binding",
            "phosphatase activator activity",
            "ATP-dependent activity",
            "protein N-terminus binding",
            "catalytic activity",
            "type II transforming growth factor beta receptor binding",
            "MAP kinase kinase activity",
            "ATPase regulator activity")
catemf <- ego@result$Description[c(1312,1313,1317,1328,1334,1345,1346,1360,1406,1420)]
dotplot(ego, showCategory= catemf) + ggtitle("Dotplot for GO:MF Enrichment")
ggsave("单独的功能富集/HCC/MF.pdf", width = 6,height = 10)

catekegg <- c("Cell cycle",
              "Cellular senescence",
              "Ubiquitin mediated proteolysis",
              "Protein processing in endoplasmic reticulum",
              "p53 signaling pathway",
              "Chemical carcinogenesis - reactive oxygen species",
              "DNA replication",
              "Hippo signaling pathway",
              "FoxO signaling pathway",
              "Proteoglycans in cancer")
dotplot(ekegg, showCategory=catekegg) + ggtitle("Dotplot for KEGG Enrichment")
ggsave("单独的功能富集/HCC/KEGG.pdf", width = 6,height = 10)
```

分面绘制dotplot
```{r}
library(stringr)
dotplot(ego, showCategory= c(catebp,catecc,catemf),split = "ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale="free")+
  ggtitle("Dotplot for GO Enrichment")+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40))  # 设置标签的最大宽度为 15

ggsave("BPCCMF.pdf",width = 7,height = 12)
```

cnet和emap图
```{r}
library(clusterProfiler)
library(enrichplot)
# 绘制基因-通路网络图
cnetplot(ego, showCategory=5, title="GO Enrichment Analysis Network")

# 绘制富集结果的关联图
emapplot(ego, showCategory=10, title="GO Enrichment Analysis Enrichment Map")

```

筛选条目绘制cnet
```{r}
cate <- c("mitotic cell cycle checkpoint signaling",
          "response to oxygen levels",
          "regulation of mitotic cell cycle phase transition",
          "cellular response to oxygen levels",
          "regulation of cell cycle phase transition")


cnetplot(GO, color.params = list(foldChange = geneList$foldChange), circular = TRUE)

p <- cnetplot(ego, showCategory=cate, title="GO Enrichment Analysis Network",circular = TRUE, colorEdge = TRUE)+   
  theme(plot.title = element_text(size = 36, face = "bold"),  # 标题字体
        axis.title = element_text(size = 15),                 # 坐标轴标题字体
        axis.text = element_text(size = 12),                  # 坐标轴刻度字体
        legend.title = element_text(size = 15),               # 图例标题字体
        legend.text = element_text(size = 12)                 # 图例文字字体
  )
# 调整标签重叠参数
p + geom_text_repel(max.overlaps = 2000)
p
ggsave("cnetplot.pdf",p,width = 32,height = 28,limitsize = FALSE)
```

upsetplot
```{r}
library(ggupset)
enrich_data <- summary(ego)
filtered_data <- enrich_data[enrich_data$Description %in% catebp,]
ego_filtered <- ego
ego_filtered@result <- filtered_data
upsetplot(ego_filtered)
enrich_data <- as.data.frame(ekegg)
filtered_data <- enrich_data[enrich_data$Description %in% catekegg,]
ekegg_filtered <- ekegg
ekegg_filtered@result <- filtered_data
upsetplot(ekegg_filtered)
upsetplot(ego, nset = 5)
```

筛选条目绘制upsetplot
```{r}
enrich_data <- summary(ego)
filtered_data <- enrich_data[enrich_data$Description %in% catebp,]
ego_filtered <- ego
ego_filtered@result <- filtered_data
upsetplot(ego_filtered)
```

其他
```{r}
# 
pathview(gene.data = gene_data, pathway.id = "hsa04110", species = "hsa")#将基因表达信息映射到 KEGG 通路图
#gene.data: 包含基因表达变化的向量。
#pathway.id: KEGG 通路的 ID，例如 "hsa04110"。
#species: 物种代码，hsa 表示人类。
#out.suffix: 输出文件的后缀。
#kegg.native: 是否输出原生 KEGG 图形。TRUE 输出原始 KEGG 图，FALSE 则输出用户自定义颜色的图。
# 
cnetplot(ego, foldChange = )
# 
heatplot(ego, foldChange = pc_enrich_list$logFC,showCategory = 5)#KEGG 通路和基因的热图
# 
ridgeplot(ekegg,showCategory = 10)#展示基因在多个通路中的分布,需要gseaResult
# 
gseaplot(ekegg, geneSetID = "hsa04110")# 需要gseaResult
```

保存结果恢复成enrichResult
```{r}
#如果是clusterProfiler的enrichGO(),gseGO(),enricher(),gseGO()等函数的结果ego保存成的文件，已关闭R环境。
#可导入文件，新建enrichResult对象ego，再进行下一步可视化。
data<-read.table("go_enrich.csv",sep="\t",header=T,quote="")

geneID_all <- unlist(apply(as.matrix(data$geneID),1,function(x) unlist(strsplit(x,'/')))) #得到富集到的所用geneID

ego<-new("enrichResult", result=data, gene=geneID_all, pvalueCutoff=0.01,pAdjustMethod="BH",qvalueCutoff=0.05,ontology="BP",keytype="GID",universe='Unknown',geneSets=list(),organism="Unknown",readable=FALSE) #把data内容赋值给ego的result，geneID_all赋值给gene，每个富集到的GO对应的gene集应该赋值给geneSets(数据是字典(键值对是GOID和geneIDs)组成的列表，这里直接给了空的)，ontology与enrichGO分析的ont参数一致，这里的pvalueCutoff=0.01,pAdjustMethod="BH",qvalueCutoff=0.05根据富集分析参数的设置，或者随意设置或者不设置也不会影响可视化。
```

