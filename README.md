# RNA-seq

## RNA-seq数据的分析流程

### 1. counts、fpkm、rpkm、tmp的区别与分析方法

[counts转换为FPKM](https://github.com/handysssun/RNA-seq/blob/main/pre-proceed/counts_to_FPKM.R)

### 2. limma、edgR、DESeq2差异表达分析的区别

2.1 差异分析

[limma](https://github.com/handysssun/RNA-seq/blob/main/DEGs/limma.R) [limma两组](https://github.com/handysssun/RNA-seq/blob/main/DEGs/limma_two_groups.Rmd)

[DESeq2](https://github.com/handysssun/RNA-seq/blob/main/DEGs/DESeq2.R)

[edgeR](https://github.com/handysssun/RNA-seq/blob/main/DEGs/edgeR.R)

2.2 差异可视化

[volcano_plot](https://github.com/handysssun/RNA-seq/blob/main/DEGs/Volcano_plot.R)

[pHeatmap](https://github.com/handysssun/RNA-seq/blob/main/DEGs/pHeatmap.R)

[histogram](https://github.com/handysssun/RNA-seq/blob/main/DEGs/histogram.R)

3. 功能富集分析

[使用clusterProfiler包进行GO和KEGG富集分析](https://github.com/handysssun/RNA-seq/blob/main/function_enrichment/clusterProfiler.R)

[clusterProfiler的富集结果可以进行可视化](https://github.com/handysssun/RNA-seq/blob/main/function_enrichment/clusterProfiler_GO.Rmd)

[也可以使用DAVID在线富集后使用R进行可视化](https://github.com/handysssun/RNA-seq/blob/main/function_enrichment/DAVID_annotion_visualization.R)

4. 富集可视化

5. GSEA

6. 共表达及网络

7. 相关性分析

[对两个基因进行相关性分析并绘制散点+回归线图](https://github.com/handysssun/RNA-seq/blob/main/cor/cor_plot.R)
