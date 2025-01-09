# 加载包、读取数据
library(ggvolcano)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(tidyverse)

degdata <- dres_diff_data
colnames(degdata)

# 设定差异基因阈值
degdata$status <- 'stable'
degdata$status[degdata$log2FoldChange>1 & degdata$padj<0.05] <- 'up'
degdata$status[degdata$log2FoldChange< -1 & degdata$padj<0.05] <- 'down'

# 挑选出差异基因子集，用于后续添加label
labeldata <- subset(degdata, abs(log2FoldChange)>1 & padj<0.05)
id <- order(-log10(labeldata$padj),decreasing = T)
labeldata <- labeldata[id,]

# (1)基础版-------------------------------------------
ggplot(degdata,aes(log2FoldChange,-log10(padj)))+
  # 绘制散点
  geom_point(aes(color=status),
             size=1.5)+
  # 绘制垂直线
  geom_vline(xintercept = c(log2(1/2),log2(2)), linetype = "dashed")+
  # 绘制水平线
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  theme_bw()+
  # 设置点的颜色
  scale_color_manual(values = c('red','grey','blue'))+
  # 设置标签注释
  # geom_text_repel(data = labeldata[1:10,], 
  #                 aes(label = label,color=status),
  #                 size=2.5)+
  # 横轴标题
  xlab('Log2FC')+
  # 纵轴标题
  ylab('-Log10(adjPvalue)')+ 
  ggtitle("Volcano Plot of Basel vs. Early",)
ggsave('result_plot/Basel_vs._Early.pdf',width = 7,height = 5.5)

# (2)进阶版-------------------------------------------
ggplot(data = degdata,
       mapping = aes(
         x=logFC,
         y=-log10(adj.P.Val)))+
  # 绘制散点
  geom_point(aes(color=status,
                 size= -log10(adj.P.Val)),
             alpha=1)+
  # 绘制垂直线
  geom_vline(xintercept = c(log2(1/2),log2(2)), linetype = "dashed")+
  # 绘制水平线
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  theme_bw()+
  # 设置点的颜色
  scale_color_manual(values = c('#ffcc00','grey','#0066ff'))+
  # 设置点的大小范围
  scale_size_continuous(range = c(0.3,3)) + 
  # 设置标签注释
  geom_text_repel(data = labeldata[1:10,],
                  mapping = aes(label = label,color=status),
                  size=2)+
  # 横轴标题
  xlab('Log2FC')+
  # 纵轴标题
  ylab('-Log10(adjPvalue)')
ggsave('plot2.pdf',width = 7,height = 5.5)

# (3)彩色渐变版-------------------------------------------
ggplot(data = degdata,
       mapping = aes(
         x=logFC,
         y=-log10(adj.P.Val)))+
  # 绘制散点
  geom_point(aes(color=-log10(adj.P.Val),
                 size= -log10(adj.P.Val)),
             alpha=1)+
  # 绘制垂直线
  geom_vline(xintercept = c(log2(1/2),log2(2)), 
             linetype = "dashed",
             color='#363636')+
  # 绘制水平线
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             color='#363636')+
  theme_bw()+
  # 设置点的颜色
  scale_color_gradientn(colours = brewer.pal(11,'RdYlBu') %>% rev())+
  # 设置点的大小范围
  scale_size_continuous(range = c(0.3,3)) + 
  # 设置标签注释
  geom_text_repel(data = labeldata[1:10,],
                  mapping = aes(label = label,
                                color=-log10(adj.P.Val)),
                  size=2.5)+
  # 横轴标题
  xlab('Log2FC')+
  # 纵轴标题
  ylab('-Log10(adjPvalue)')+
  # 去除一个图例，设置另一个图例的标题
  guides(size=FALSE,
         color=guide_colorbar(title = '-Log10(adjPvalue)')) # guide_colorbar是设置连续型变量图例的
ggsave('plot3.pdf',width = 7,height = 5.5)
