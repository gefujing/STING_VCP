
# 设置环境
rm(list = ls())  
library(clusterProfiler)
library(ggthemes)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)
library(enrichplot)
library(GseaVis)

# 加载数据
deg = read.csv(file = "./3-BMDM/1-degs-bvsa.csv")

# ID转换
s2e = bitr(deg$gene_name, 
           fromType = "SYMBOL",
           toType = "ENTREZID",
           OrgDb = org.Mm.eg.db)

deg = inner_join(deg, s2e, by=c("gene_name" = "SYMBOL"))

# # 输入数据
# gene_diff = deg$ENTREZID[deg$change != "stable"] 

# # KEGG富集
# ekk = enrichKEGG(gene = gene_diff, organism = "hsa")
# ekk = setReadable(ekk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
# ekk_table = ekk@result
# write.csv(x = ekk_table, file = "./kegg_enrich.csv")

# # KEGG富集可视化
# p1 = barplot(ekk)
# p1
# ggsave(plot = p1, filename = "./kegg_enrich.pdf", width = 4.5, height = 3.8)
# 
# # GO富集分析
# ego = enrichGO(gene = gene_diff, OrgDb = org.Hs.eg.db, ont = "ALL", readable = TRUE)
# ego_table = ego@result
# write.csv(x = ego_table, file = "./go_enrich.csv")
# 
# # GO可视化
# dotplot(ego, split = "ONTOLOGY") + 
#   facet_grid(ONTOLOGY ~ ., space = "free_y",scales = "free_y") 

# GSEA富集分析

## 数据构建
deg = deg[order(deg$log2.fc., decreasing = T),]
genelist = deg$log2.fc.
names(genelist) = deg$ENTREZID

## 利用KEGG进行富集
gsea_kk = gseKEGG(
  geneList = genelist, # 根据logFC排序的基因集
  organism = "mmu",    # 人的拉丁名缩写
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

## 导出表格
gsea_kk_table = gsea_kk@result
write.csv(gsea_kk_table, file = "./3-BMDM/2-kegg-dvsall.csv")

## 可视化
geneSetID = c("mmu04623", "mmu05164", "mmu04620")
gseaNb(object = gsea_kk,
       geneSetID = geneSetID,
       subPlot = 2,
       termWidth = 35,
       legend.position = c(0.8,0.8),
       addPval = T,
       pvalX = 0.05, pvalY = 0.05)



## 利用GO进行富集
gsea_go = gseGO(
  geneList = genelist, # 根据logFC排序的基因集
  ont = "ALL", # 可选"BP"、"MF"、"CC"三大类或"ALL"
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)


## 导出表格
gsea_go_table = gsea_go@result
write.csv(gsea_go_table, file = "./3-BMDM/2-go-dvsall.csv")

## 可视化
geneSetID = c("GO:0035458", "GO:0051607", "GO:0002263", "GO:0035456", "GO:0034340")
gseaNb(object = gsea_go,
       geneSetID = geneSetID,
       subPlot = 2,
       termWidth = 35,
       legend.position = c(0.8,0.8),
       addPval = T,
       pvalX = 0.05, pvalY = 0.05,
       curveCol = jjAnno::useMyCol('paired', 5))






