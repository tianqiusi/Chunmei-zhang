setwd("D:/R/富集分析结果美化")

# 下载所需要的R包 ----------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# 检查 clusterProfiler 包是否已经安装
if (!require("clusterProfiler", quietly = TRUE)) {
  # 如果未安装，则使用 BiocManager 安装
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("clusterProfiler")
}
# 检查 org.Hs.eg.db 包是否已经安装
if (!require("org.Hs.eg.db", quietly = TRUE)) {
  # 如果未安装，则使用 BiocManager 安装
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("org.Hs.eg.db")
}
# 检查 ggplot2 包是否已经安装
if (!require("ggplot2", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("ggplot2")
}
# 检查 dplyr 包是否已经安装
if (!require("dplyr", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("dplyr")
}
# 检查 ggsci 包是否已经安装
if (!require("ggsci", quietly = TRUE)) {
  # 如果未安装，则使用 install.packages 函数安装
  install.packages("ggsci")
}
# 加载所需要的R包 ----------------------------------------------------------------
library(clusterProfiler) #富集分析主要的包
library(org.Hs.eg.db)#查找物种注释信息
library(ggplot2)#分面绘图所需
library(dplyr)#整理数据
library(ggsci)#配色

# 读取数据：基因差异分析的结果 ----------------------------------------------------------
upGene <- read.csv("upGene_1_0.05.csv",row.names = 1)
Genes <- bitr(rownames(upGene),  
              fromType = "SYMBOL", #输入数据的类型
              toType = c("ENTREZID"), #要转换的数据类型
              OrgDb = org.Hs.eg.db) #物种

# 富集分析 --------------------------------------------------------------------
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
GO_result <- as.data.frame(GO)


# 基础柱状图 -------------------------------------------------------------------
barplot(GO)
barplot(GO, drop = TRUE, 
        showCategory =6,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')


# 或许你想展示某些特定条目 ------------------------------------------------------------
#比如挑选每个分类中前10显著的条目
BP_top10 <- GO_result %>%
  filter(ONTOLOGY == "BP") %>%  # 筛选 ONTOLOGY 为 BP 的子集
  arrange(p.adjust) %>%        # 按照 p.adjust 的升序排列
  head(10)                      # 提取前 10行

CC_top10 <- GO_result %>%
  filter(ONTOLOGY == "CC") %>%  # 筛选 ONTOLOGY 为 CC 的子集
  arrange(p.adjust) %>%        # 按照 p.adjust 的升序排列
  head(10)                      # 提取前 10行

MF_top10 <- GO_result %>%
  filter(ONTOLOGY == "MF") %>%  # 筛选 ONTOLOGY 为 MF 的子集
  arrange(p.adjust) %>%        # 按照 p.adjust 的升序排列
  head(10)                      # 提取前 10行

merge_data <- rbind(BP_top10,CC_top10,MF_top10)
merge_data$ONTOLOGY <- factor(merge_data$ONTOLOGY, levels = c("BP", "CC","MF"))
merge_data$logPvalue <- (-log(merge_data$p.adjust))
# 美化版本 --------------------------------------------------------------------
colnames(merge_data)
merge_data$Description <- factor(merge_data$Description, levels = unique(merge_data$Description[order(merge_data$ONTOLOGY)]))

ggplot(merge_data, aes(x = Description, y = logPvalue, fill = ONTOLOGY)) + #指定数据为merge_data，以Description作为x轴，以logPvalue作为y轴，以ONTOLOGY作为颜色填充。
  geom_bar(stat = "identity", width = 0.8, color = "black") + #添加一个柱状图层，使用identity统计方法（不进行任何变换），每个柱子的宽度为0.8，颜色为黑色
  scale_x_discrete(limits = unique(merge_data$Description[order(merge_data$ONTOLOGY)])) + #设置x轴的标签顺序，根据ONTOLOGY对Description排序，以便按照ONTOLOGY的顺序显示。unique()函数保证每个标签只出现一次。
  coord_flip() + #翻转坐标轴，交换x轴和y轴，使得柱状图变为横向。
  scale_y_continuous(expand = c(0,0))+ #设置y轴的范围为连续型变量（continuous），设置expand参数为c(0,0)，意味着y轴不需要扩展（不需要增加额外的空白）。
  theme(panel.background = element_rect(fill = "white"))+ #设置背景为白色，使用element_rect()函数设置元素的矩形形状
  theme(axis.text = element_text(size = 16, family = "Arial"),
        axis.title = element_text(size=18,family="Arial",colour = "black"),
        legend.text = element_text(size = 16, family = "Arial"),
        legend.title = element_text(size = 18, family = "Arial"),
        legend.position = "top")+ #设置各种文本元素的字体、大小和颜色，包括坐标轴文本、坐标轴标题、图例文本和图例标题。将图例放置在顶部。
  labs(x = "GO term", y = "-log P-value")+
  scale_fill_nejm(alpha = 0.8) #使用NEJM（New England Journal of Medicine）颜色调色板，设置颜色透明度为0.8，对填充颜色进行比例缩放。
