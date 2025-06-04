install.packages("pheatmap") #安装pheatmap包
library(pheatmap) #加载pheatmap包
#?pheatmap #查看pheatmap包里面的详细介绍
#?pheatmap::pheatmap #查看pheatmap包里pheatmap函数的具体参数

#简单绘制相关性热图
#数据框格式可直接使用cor函数求相关系数
matrix <- cor(data) 
#基于相关系数作图
pheatmap(matrix)
#在热图的单位格中显示数值
pheatmap(matrix,display_numbers=T)



#直接绘制
pheatmap(test)

#设置col、row方向的树高
pheatmap(test,treeheight_row=100,treeheight_col=20)

#取消列聚类，并且更改颜色
pheatmap(test,
         treeheight_row=100,
         treeheight_col=20,
         cluster_cols=FALSE,
         color=colorRampPalette(c("green","black","red"))(1000))  
		 
#自定义图例颜色范围
col_fun = colorRamp2(c(0,1,6),#范围
                     c("blue", "white", "red"))#将蓝色赋给0，白色赋给1，红色赋给6，设定完成后将col_fun赋给color参数

#取消单元格间的边框，调整字体大小，并且保存在桌面文件中
pheatmap(test,
         treeheight_row=100,
         treeheight_col=20,
         cluster_cols=FALSE,
         color=colorRampPalette(c("green","black","red"))(1000),
         border_color=NA,
         fontsize=10,
         fontsize_row=8,
         fontsize_col=16,
         file='C:/Users/xu/Desktop/test.jpg')   


#增加分组信息，使得pheatmap显示行或列的分组信息
#这部分以及之后的内容参考了第四篇参考文献
annotation_col = data.frame(CellType = factor(rep(c("X1", "X2"), 5)), Time = 1:5)  #增加Time，CellType分组信息
rownames(annotation_col) = paste("Test", 1:10, sep = "")   
annotation_row = data.frame(GeneClass = factor(rep(c("P1", "P2", "P3"), c(10, 7, 3))))  #增加GeneClass分组信息
rownames(annotation_row) = paste("Gene", 1:20, sep = "") 
pheatmap(test, annotation_col = annotation_col, annotation_row = annotation_row)
#使用annotation_colors参数设定各个分组的颜色  
ann_colors = list(Time = c("white", "green"),cellType = c(X1= "#1B9E77", X2 = "#D95F02"),GeneClass = c(P1 = "#7570B3", P2 = "#E7298A", P3 = "#66A61E"))   
pheatmap(test, annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = ann_colors) 
# cutree_rows, cutree_cols可以根据行列的聚类数将热图分隔开；
pheatmap(test,cutree_rows=3,cutree_cols=2)

#查看聚类结果
result <- pheatmap(test)
summary(result)

# 行的聚类排列顺序
result$tree_row$order 
#  得到行名的顺序
rownames(test)[result$tree_row$order]
# 查看按行聚类后的热图顺序结果
head(test[result$tree_row$order,])
# 查看按照行和列聚类之后，得到的热图顺序结果
head(test[result$tree_row$order,result$tree_col$order])

# 绘图前进行zscore标准化
df_zscore <- apply(plot_data, 1, function(row){
  (row - mean(row)) / sd(row)
})

plot_data <- t(df_zscore)


bk <- seq(-1.5,1.5, length.out=1000)# 设置颜色范围
p <- pheatmap(plot_data,cluster_cols = FALSE,
              color = colorRampPalette(colors = c("blue","white","red"))(length(bk)),
              breaks = bk,
              show_colnames = TRUE,
              show_rownames = TRUE,
              cluster_rows = TRUE,
              scale = "column",
              cellheight = 3,
              fontsize_row = 3
              
              
)

matrix2 <- round(t(apply(plot_data,1, scale)),2)# 由于之前的绘图过程中对数据进行了标准化，此处手动按照行标准化，取两位小数
colnames(matrix2) <- colnames(plot_data)

exprTable <- as.data.frame(matrix2)#若是镜像列名需要转置matrix,行不需要
row_dist = dist(exprTable,method = "euclidean")
hclust_1 <- hclust(row_dist)# 生成聚类文件
p#最初的热图存于一个变量中
p$tree_row$order#聚类后的行名索引,为数字
old <- rownames(plot_data)[p$tree_row$order]#最初热图按聚类结果排序的行名

#我自己的案例是把前1到15和后16到31换一下顺序
manual_order = old[c(44:101,1:43)]#修改后的镜像顺序
dend = reorder(as.dendrogram(hclust_1), wts=order(match(manual_order, rownames(exprTable))))# 重新生成聚类文件
row_cluster <- as.hclust(dend)

# 将数据进行分区聚类
# 切分
group1 <- plot_data[1:52, ]     # 第一部分
group2 <- plot_data[53:554, ]   # 第二部分
# 分别进行聚类（层次聚类为例）
# 使用欧氏距离 + 完全链接
hc1 <- hclust(dist(group1), method = "complete")
hc2 <- hclust(dist(group2), method = "complete")
# 结合聚类结果排列行顺序
# 获取聚类后的行顺序
order1 <- hc1$order
order2 <- hc2$order

# 合并排序
new_order <- c(rownames(group1)[order1], rownames(group2)[order2])

# 按新顺序重新排列原数据
plot_data_ordered <- plot_data[new_order, ]


pheatmap(plot_data,cluster_cols = FALSE,
         color = colorRampPalette(colors = c("blue","white","red"))(length(bk)),
         breaks = bk,
         show_colnames = TRUE,
         show_rownames = TRUE,
         cluster_rows = row_cluster,
         scale = "column",
         cellheight = 3,
         fontsize_row = 3
         
         
)
