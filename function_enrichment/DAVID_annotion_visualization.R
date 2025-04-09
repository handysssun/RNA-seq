# 载入DAVID富集结果
david_anno <- read.table("4_DAVID_annotion_chart_E8B89F6D8AA01744182996757.txt",header = TRUE,sep = "\t")

# 指定顺序向量
order_vec <- c("GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY")

library(stringi)
library(ggplot2)
library(dplyr)

enrich <- david_anno# read and reorder annotion result
enrich_signif=enrich[which(enrich$PValue<0.05),]# significant term
enrich_signif=enrich_signif[,c(1:3,5)]# keep Categroy/Term/PValue
head(enrich_signif)
enrich_signif=data.frame(enrich_signif)
enrich_signif_sorted <- enrich_signif[order(factor(enrich_signif[[1]], levels = order_vec), enrich_signif$PValue), ]# reorder by counts and category


KEGG = enrich_signif_sorted[enrich_signif_sorted$Category == "KEGG_PATHWAY",][1:10,]# Extracting partial results-KEGG
KEGG$Term<-stri_sub(KEGG$Term,10,100)# sub_strings from $2 to $3，stringi包
KEGG <- KEGG[order(KEGG$Count,decreasing = FALSE),]# reorder by count 
KEGG$Term <- factor(KEGG$Term,levels = KEGG$Term)# order Terms with factor type


ggplot(KEGG,aes(x=Count,y=Term))+
  geom_point(aes(color=PValue,size=Count))+
  scale_color_gradient(low='slateblue4',high='firebrick3')+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  theme_light()+
  theme(panel.grid.minor = element_line(color = "grey80"),panel.grid.major = element_line(color = "grey90"))
ggsave("5_kegg_dotplot.pdf",height = 7,width = 5)


BP = enrich_signif_sorted[enrich_signif_sorted$Category == "GOTERM_BP_DIRECT",][1:10,]# Extracting partial results-BP
BP$Term<-stri_sub(BP$Term,12,100)# sub_strings from $2 to $3，stringi包
BP <- BP[order(BP$Count,decreasing = FALSE),]# reorder by count 
BP$Term <- factor(BP$Term,levels = BP$Term)# order Terms with factor type

ggplot(BP,aes(x=Count,y=Term))+
  geom_point(aes(color=PValue,size=Count))+
  scale_color_gradient(low='slateblue4',high='firebrick3')+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  theme_light()+
  theme(panel.grid.minor = element_line(color = "grey80"),
        panel.grid.major = element_line(color = "grey90")
  )
ggsave("5_bp_dotplot.pdf",height = 7,width = 5)

CC = enrich_signif_sorted[enrich_signif_sorted$Category == "GOTERM_CC_DIRECT",][1:10,]# Extracting partial results-CC
CC$Term<-stri_sub(CC$Term,12,100)# sub_strings from $2 to $3，stringi包
CC <- CC[order(CC$Count,decreasing = FALSE),]# reorder by count 
CC$Term <- factor(CC$Term,levels = CC$Term)# order Terms with factor type

ggplot(CC,aes(x=Count,y=Term))+
  geom_point(aes(color=PValue,size=Count))+
  scale_color_gradient(low='slateblue4',high='firebrick3')+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  theme_light()+
  theme(panel.grid.minor = element_line(color = "grey80"),
        panel.grid.major = element_line(color = "grey90")
  )
ggsave("5_cc_dotplot.pdf",height = 7,width = 5)

MF = enrich_signif_sorted[enrich_signif_sorted$Category == "GOTERM_MF_DIRECT",][1:10,]# Extracting partial results-MF
MF$Term<-stri_sub(MF$Term,12,100)# sub_strings from $2 to $3，stringi包
MF <- MF[order(MF$Count,decreasing = FALSE),]# reorder by count 
MF$Term <- factor(MF$Term,levels = MF$Term)# order Terms with factor type

ggplot(MF,aes(x=Count,y=Term))+
  geom_point(aes(color=PValue,size=Count))+
  scale_color_gradient(low='slateblue4',high='firebrick3')+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  theme_light()+
  theme(panel.grid.minor = element_line(color = "grey80"),
        panel.grid.major = element_line(color = "grey90")
  )
ggsave("5_mf_dotplot.pdf",height = 7,width = 5)