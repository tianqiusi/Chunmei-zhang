setwd("D:/R/富集分析")

# 加载富集分析所需要的R包 ------------------------------------------------------------
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")

library(clusterProfiler) #富集分析主要的包
library(org.Hs.eg.db)#查找物种注释信息
library(ggplot2)#分面绘图所需

# 读取数据：基因差异分析的结果 ----------------------------------------------------------
upGene <- read.csv("upGene_1_0.05.csv",row.names = 1)
Genes <- bitr(rownames(upGene),  
                    fromType = "SYMBOL", #输入数据的类型
                    toType = c("ENTREZID"), #要转换的数据类型
                    OrgDb = org.Hs.eg.db) #物种

# GO富集分析与结果存储 -------------------------------------------------------------
GO <- enrichGO(gene = Genes$ENTREZID, #输入基因的"ENTREZID"
               OrgDb = org.Hs.eg.db,#注释信息
               keyType = "ENTREZID",
               ont = "all",     #可选条目BP/CC/MF
               pAdjustMethod = "BH", #p值的校正方式
               pvalueCutoff = 1,   #pvalue的阈值
               qvalueCutoff = 1, #qvalue的阈值
               minGSSize = 5,
               maxGSSize = 5000,
               readable = TRUE)   #是否将entrez id转换为symbol
write.table(data.frame(ID=row.names(GO@result),GO@result),
            file = "GO_enrich_result.txt",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

# 结果可视化--柱状图 -----------------------------------------------------------
pdf(file="GO柱状图.pdf",width = 10,height = 8)
barplot(GO, drop = TRUE, 
        showCategory =6,split="ONTOLOGY") + 
        facet_grid(ONTOLOGY~., scale='free')
dev.off()
head(GO@result)

# 结果可视化--点图 ---------------------------------------------------------------
pdf(file="GO点图.pdf",width = 10,height = 8)
dotplot(GO,showCategory = 6,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')  #方式一：分面
dotplot(
  GO,
  x = "GeneRatio",
  color = "p.adjust",
  showCategory = 15,
  size = NULL,
  split = NULL,
  font.size = 12,
  title = "",
  orderBy = "x",
  label_format = 30)  #方式二：不分面

# KEGG富集分析与结果存储 -----------------------------------------------------------
KEGG <- enrichKEGG(gene = Genes$ENTREZID,
                   organism = "hsa", 
                   keyType = "kegg", #KEGG数据库
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)

barplot(KEGG,
        x = "GeneRatio",
        color = "p.adjust",
        showCategory = 15,
        title = "KEGG_enrichment")#标题
dotplot(KEGG)
View(as.data.frame(KEGG))#以数据框的形式展示KEGG的结果
browseKEGG(KEGG,"hsa04060")#打开KEGG官网查看通路的详细信息

####结果存储
KEGG_results <- DOSE::setReadable(KEGG,
                                  OrgDb = org.Hs.eg.db,
                                  keyType = "ENTREZID") #将ENTREZID转换为genesymbol
write.table(data.frame(ID=rownames(KEGG_results@result),KEGG_results@result),
            file = "KEGG_enrich_result.txt",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)
