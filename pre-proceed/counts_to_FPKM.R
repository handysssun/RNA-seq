# counts FPKM RPKM TPM之间的转换
library(GenomicFeatures)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(tibble)
library(readxl)

# 加载参考基因组文件，计算基因长度
txdb <- makeTxDbFromGFF("resource_data/Homo_sapiens.GRCh38.115.gff3.gz")
exonic <- exonsBy(txdb, by="gene")
exonic.gene.sizes <- sum(width(reduce(exonic)))
mygeneids <- data.frame(gene_id=names(exonic.gene.sizes), gene_length=exonic.gene.sizes)

# 如果注释是ENSEMBLID，需要转换为SYMBOL
mygeneids$gene_id <- gsub("\\..*", "",mygeneids$gene_id) #去除版本号
mygeneids$gene_symbol <- mapIds(org.Hs.eg.db,
                                keys=mygeneids$gene_id,
                                column="SYMBOL",
                                keytype="ENSEMBL",
                                multiVals="first") #名称转换由ENSEMBL转换为SYMBOL
rownames(mygeneids) <- c(1:nrow(mygeneids))
library(dplyr)
gene_length <- mygeneids %>%
  distinct(gene_symbol,.keep_all = T) %>%
  filter(!is.na(gene_symbol)) %>%
  select(-gene_id)
# gene_length <- cbind(gene_length$gene_symbol,gene_length$gene_length)


# 如果是symbol就直接转换
# gene_length <- mygeneids
# colnames(gene_length) <- c("gene_symbol", "gene_length")
# gene_length <- data.frame(gene_length)

# 加载整理readcount数据，第一列为gene_symbol，不能有重复
df <- eset
rc <- data.frame(gene_symbol = rownames(eset),eset)
colnames(rc)[1] <- "gene_symbol"
merge_df <- merge(rc, gene_length, by = "gene_symbol")

counts <- merge_df %>%
  select(1:10) %>%
  column_to_rownames("gene_symbol")

counts$gene_length <- as.numeric(counts$gene_length)
counts <- as.matrix(counts)

# counts to CPM
cpm <- apply(X = subset(counts, select = c(-gene_length)), 
             MARGIN = 2, 
             FUN = function(x) {
               x / sum(as.numeric(x)) * 10^6
             })

# counts to RPKM/FPKM
rpkm <- apply(X = subset(counts, select = c(-gene_length)),
              MARGIN = 2, 
              FUN = function(x) {
                10^9 * x / subset(counts, select = c(gene_length)) / sum(as.numeric(x))
              })

# counts to TPM
rpk <- apply(subset(counts, select = c(-gene_length)), 2, 
             function(x) x / (geneLengths / 1000))
tpm <- apply(X = rpk, 
             MARGIN = 2, 
             FUN = function(x) {
               x / sum(as.numeric(x)) * 10^6
             })

# FPKM/RPKM to TPM
fpkmToTpm <- function(fpkm) {
  if (any(fpkm <= 0)) {
    stop("FPKM values must be greater than zero.")
  }
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

# 转换成FPKM
countToFpkm <- function(counts, effLen){
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

fpkms <- apply(counts[,1:ncol(counts) - 1], 2, countToFpkm, effLen = counts$gene_length)# 传递counts和单独的geneLength列
fpkms.m<-data.frame(fpkms)
# colnames(fpkms.m) <- colnames(eset)
dim(fpkms.m)
#查看前三个基因的FPKM值
fpkms.m[1:3,]

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

fpkms <- apply(counts[,1:ncol(counts) - 1], 2, countToFpkm, effLen = counts$gene_length)# 传递counts和单独的geneLength列
fpkms.m<-data.frame(fpkms)
# colnames(fpkms.m) <- colnames(eset)
dim(fpkms.m)
#查看前三个基因的FPKM值
fpkms.m[1:3,]
write.csv(fpkms.m,"ACC_FPKM.csv")
