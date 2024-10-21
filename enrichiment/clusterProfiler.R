# 20240725富集分析
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
library(ggVolcano)
library(ggrepel)
library(RColorBrewer)
library(tidyverse)

degs <- read.xlsx("GO/logFC1_DEGs.xlsx")
# 指定富集分析的物种库
GO_database <- 'org.Hs.eg.db' #GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
KEGG_database <- 'hsa' #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html

# gene ID转换
# keytypes(org.Hs.eg.db)
# [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
# [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
# [11] "GENETYPE"     "GO"           "GOALL"        "IPI"          "MAP"         
# [16] "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"        
# [21] "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
# [26] "UNIPROT"

# 分析前先准备gene_list
library(clusterProfiler)
library(org.Hs.eg.db)
gene <- bitr(degs$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID', OrgDb = GO_database)
gene_list<- gene$ENTREZID
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
  write.csv(as.data.frame(GO),"clsuterprofiler_GO.csv",row.names = FALSE)
  KEGG <- enrichKEGG(gene_list,
                     organism = KEGG_database,
                     pvalueCutoff = p_cutoff,
                     qvalueCutoff = q_cutoff)
  write.csv(as.data.frame(KEGG), "clsuterprofiler_KEGG.csv", row.names = FALSE)
}
clusterprofiler_GO_enrichment(input_file = gene_list)


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
					pvalueCutoff = 0.1, #指定 p 值阈值（可指定 1 以输出全部）
					qvalueCutoff = 0.2) #指定 q 值阈值（可指定 1 以输出全部）
			  
#enrichplot可视化

# 保存的clusterProfiler富集分析结果做可视化
#如果是clusterProfiler的enrichGO(),gseGO(),enricher(),gseGO()等函数的结果ego保存成的文件，已关闭R环境。
#可导入文件，新建enrichResult对象ego，再进行下一步可视化。
data<-read.table("go_enrich.csv",sep="\t",header=T,quote="")

geneID_all <- unlist(apply(as.matrix(data$geneID),1,function(x) unlist(strsplit(x,'/')))) #得到富集到的所用geneID

ego<-new("enrichResult", result=data, gene=geneID_all, pvalueCutoff=0.01,pAdjustMethod="BH",qvalueCutoff=0.05,ontology="BP",keytype="GID",universe='Unknown',geneSets=list(),organism="Unknown",readable=FALSE) #把data内容赋值给ego的result，geneID_all赋值给gene，每个富集到的GO对应的gene集应该赋值给geneSets(数据是字典(键值对是GOID和geneIDs)组成的列表，这里直接给了空的)，ontology与enrichGO分析的ont参数一致，这里的pvalueCutoff=0.01,pAdjustMethod="BH",qvalueCutoff=0.05根据富集分析参数的设置，或者随意设置或者不设置也不会影响可视化。

# 绘制条形图
barplot(ego, showCategory=10, title="GO Enrichment Analysis")

# 绘制气泡图
dotplot(ego, showCategory=10, title="GO Enrichment Analysis")

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

dotplot(ego, showCategory= c(catebp,catecc,catemf),split = "ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale="free")+
  ggtitle("Dotplot for GO Enrichment")+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40))  # 设置标签的最大宽度为 15

ggsave("BPCCMF.pdf",width = 7,height = 12)
# 绘制基因-通路网络图
cnetplot(ego, showCategory=5, title="GO Enrichment Analysis Network")

# 绘制富集结果的关联图
emapplot(ego, showCategory=10, title="GO Enrichment Analysis Enrichment Map")


# 创建可视化函数-点图
plot_go <- function(df, title,OUTPUT_FILE) {
  ggplot(df, aes(y = reorder(Description, -pvalue), x = Count)) +
    geom_point(aes(size = Count, color = pvalue)) +
    scale_color_gradient(low = "red", high = "blue") +
    labs(title = title,
         y = "GO Term",
         x = "Count") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),       # 标题居中
      axis.title.x = element_text(hjust = 0.5),     # X轴标题居中
      axis.title.y = element_text(hjust = 0.5),     # Y轴标题居中
      axis.text.x = element_text(angle = 45, hjust = 1)  # X轴标签倾斜
    )+
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40) )
  ggsave(OUTPUT_FILE,height = 7,width = 7)
}

# 分别对 BP、CC、MF 进行可视化
plot_bp <- plot_go(go_bp[1:20,], "GO Enrichment Dotplot: Biological Process (BP)","BP.pdf")

plot_cc <- plot_go(go_cc[1:20,], "GO Enrichment Dotplot: Cellular Component (CC)","CC.pdf")

plot_mf <- plot_go(go_mf[1:20,], "GO Enrichment Dotplot: Molecular Function (MF)","MF.pdf")

# 挑选的条目进行绘图
# https://cloud.tencent.com/developer/article/2319648
cate <- c("mononuclear cell differentiation",
          "leukocyte mediated immunity",
          "immune response-activating signal transduction",
          "activation of immune response",
          "positive regulation of cytokine production")


cnetplot(GO, color.params = list(foldChange = geneList$foldChange), circular = TRUE)

p <- cnetplot(GO, showCategory=cate, title="GO Enrichment Analysis Network",circular = TRUE, colorEdge = TRUE)+   
  theme(plot.title = element_text(size = 20, face = "bold"),  # 标题字体
        axis.title = element_text(size = 15),                 # 坐标轴标题字体
        axis.text = element_text(size = 12),                  # 坐标轴刻度字体
        legend.title = element_text(size = 15),               # 图例标题字体
        legend.text = element_text(size = 12)                 # 图例文字字体
  )
# 调整标签重叠参数
p + geom_text_repel(max.overlaps = 500)
p
ggsave("cnetplot.pdf",p,width = 34,height = 30,limitsize = FALSE)

catekegg <- KEGG@result$Description
catekegg <- catekegg[c(5,1,11,15,4,24,34,30,36,38)]

dotplot(KEGG,showCategory=catekegg)  + ggtitle("dotplot for KEGG")
ggsave("keggdot.pdf",width = 8,height = 12)






# clusterProfiler包GOKEGG可视化
library(clusterProfiler)
library(enrichplot)
dotplot(ego, showCategory=20) + ggtitle("Dotplot for KEGG Enrichment")
barplot(kegg_result, showCategory=20, title="Barplot for KEGG Enrichment")
# 
ekegg <- pairwise_termsim(ekegg)    #相关性  
emapplot(ekegg, showCategory = 20, layout = "kk", cex_category = 1.2, cex_line = 1)
emapplot(ekegg,cex_category=1.5,layout= "kk")
# 
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
