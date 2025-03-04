# counts FPKM RPKM TPM之间的转换

# 加载参考基因组文件，计算基因长度
### 获取基因长度
library(GenomicFeatures) #加载R包
## 1.读取GTF文件
txdb <- makeTxDbFromGFF("E:/R_script/genomic.gtf") #TxDb：R语言中用于储存gtf文件的一种格式
## 2.提取外显子信息
exonic <- exonsBy(txdb, by="gene") #可以提取基因外显子部分，计算counts数据比对的为基因外显子序列
## 3.外显子长度求和
exonic.gene.sizes <- sum(width(reduce(exonic))) #reduce可以删除外显子重复序列部分
## 4.结果整理
mygeneids <- data.frame(gene_id=names(exonic.gene.sizes),width=exonic.gene.sizes)

library(org.Mm.eg.db) #人org.Hs.eg.db
library(org.Hs.eg.db)
library(AnnotationDbi)
# 如果注释是ENSEMBLID，需要转换为SYMBOL
# mygeneids$gene_id <- gsub("\\..*", "",mygeneids$gene_id) #去除版本号
# mygeneids$gene_symbol <- mapIds(org.Mm.eg.db,
#                                 keys=mygeneids$gene_id,
#                                 column="SYMBOL",
#                                 keytype="ENSEMBL",
#                                 multiVals="first") #名称转换由ENSEMBL转换为SYMBOL
# rownames(mygeneids) <- c(1:nrow(mygeneids))
# library(dplyr)
# gene_length <- mygeneids %>% 
#   distinct(gene_symbol,.keep_all = T) %>% 
#   filter(!is.na(gene_symbol)) %>%
#   select(-gene_id)
# gene_length <- cbind(gene_length$gene_symbol,gene_length$width)
# 如果是symbol就直接转换
# gene_length <- mygeneids
colnames(gene_length) <- c("gene_symbol","width")
gene_length <- data.frame(gene_length)

SCI_fpkm_count <- read_excel("data/GSE196928_Complete_geneList.xlsx") #加载数据
SCI_fpkm_count <- SCI_fpkm_count[,c(1,7:9,35:37)] #选择自己需要的数据，本次只选择了三组couts与fpkm数据
colnames(SCI_fpkm_count)[1] <- "gene_symbol"
SCI_fpkm_count_1 <- merge(SCI_fpkm_count,gene_length,by = "gene_symbol") #将数据与基因长度合并
counts <- SCI_fpkm_count_1 %>%
  select(c(1,5:8)) %>%
  column_to_rownames("gene_symbol") #提取counts数据
counts$width <- as.numeric(counts$width)
counts <- as.matrix(counts) #转换为矩阵
View(counts) #看一下counts数据格式

# counts to CPM
cpm <- apply(X =subset(counts, select = c(-width)), 
             MARGIN =2, 
             FUN =function(x){
               x/sum(as.numeric(x)) * 10^6
             })

# counts to RPKM/FPKM
geneLengths <- as.vector(subset(counts, select = c(width))) # 创建基因长度的向量，为下一步标准化做准备  
rpkm <- apply(X = subset(counts, select = c(-width)),
              MARGIN = 2, 
              FUN = function(x) {
                10^9 * x / geneLengths / sum(as.numeric(x))
              }) # 计算rpkm

# counts to TPM
rpk <- apply( subset(counts, select = c(-width)), 2, 
              function(x) x/(geneLengths/1000)) #求基因长度归一化值
tpm <- apply(X = rpk, 
             MARGIN = 2, 
             FUN = function(x){
               x / sum(as.numeric(x)) * 10^6
             })#使用 rpk 值根据样本大小进行标准化

# FPKM/RPKM to TPM
fpkmToTpm <- function(fpkm) {
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}## rpkm转换为tpm


# 另一种方法
# 读取基因长度文件，并对长度文件和表达文件进行排序和整理
length <- read.csv("gene_length.csv",header = T,sep = ",")
head(length)
rownames(length) <- length$X
head(length)
rownames(length) <- substr(rownames(length),1,15)
rownames(data) <- substr(rownames(data),1,15)
data[1:5,1:5]
head(length)
inter <- intersect(rownames(data),rownames(length))
data <- data[inter,]
length <- length[inter,]
data[1:5,1:5]
identical(rownames(length),rownames(data))
if (identical(rownames(length), rownames(data))){
  print("GTF and expression matix now have the same gene and gene in same order")
}

# 转换成FPKM
countToFpkm <- function(counts, effLen){
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

fpkms <- apply(data, 2, countToFpkm, effLen = length$eff_length)
fpkms.m<-data.frame(fpkms)
colnames(fpkms.m)<-colnames(data)
dim(fpkms.m)
#查看前三个基因的TPM值
fpkms.m[1:3,]
write.csv(fpkms.m,"ACC_FPKM.csv")





