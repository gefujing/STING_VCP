
# 准备环境
rm(list=ls()) 
options(stringsAsFactors = F) 
library(Seurat)
library(ggplot2)
library(stringr)
library(tidyverse)
library(dplyr)
library(stringr)
library(data.table)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)
library(enrichplot)
library(GSVA)
library(msigdbr)

load("./5-episub/episub.Rdata")
load(file = "./1-qc/mycol.Rdata")
sce.markers = read.csv(file = "./5-episub/episub.degs.csv", header = T)
rownames(sce.markers) = sce.markers$X

# 功能聚类
features = c("OAS1", "TMEM173", "IRF3", "IRF5", "CXCL9", "CXCL10", "CCL5", "CCL13", "STAT1", "IFNAR1", "TLR4", "IFI27")
VlnPlot(episub, features = features, group.by = "group", cols = col_vec, ncol = 6, pt.size = 0)
ggsave(filename = "./5-episub/1-sce.markerss-vln.pdf", width = 12, height = 6.7) 


## KEGG
## 获取上下调基因
gene_up = rownames(sce.markers[sce.markers$avg_log2FC > 0,])
gene_down = rownames(sce.markers[sce.markers$avg_log2FC < 0,])
## 把SYMBOL改为ENTREZID
gene_up = as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db, keys = gene_up, columns = 'ENTREZID', keytype = 'SYMBOL')[,2]))
gene_down = as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db, keys = gene_down, columns = 'ENTREZID', keytype = 'SYMBOL')[,2]))
gene_diff = c(gene_up, gene_down)
## KEGG
gene_diff = unique(gene_diff)
kk.diff = enrichKEGG(gene = gene_diff, organism = "hsa", pvalueCutoff = 0.9, qvalueCutoff = 0.9)
kkdiff = kk.diff@result
write.csv(kkdiff, file = "./5-episub/episub.kegg.csv")

pathways = c("Th17 cell differentiation", "T cell receptor signaling pathway", "Th1 and Th2 cell differentiation",
             "Toll-like receptor signaling pathway", "Fc gamma R-mediated phagocytosis", "Chemokine signaling pathway", 
             "IL-17 signaling pathway", "NOD-like receptor signaling pathway", "Leukocyte transendothelial migration",
             "Cytosolic DNA-sensing pathway")
dotplot(kk.diff, showCategory = pathways, color = "pvalue")
ggsave(filename = "./5-episub/2-sce.markers-KEGG.pdf", width = 5.5, height = 4.5)


##GSEA
## 上一步差异分析得到差异基因列表sce.markers后取出，p值和log2FC
nrsce.markers = sce.markers[ ,c("avg_log2FC", "p_val")]
colnames(nrsce.markers) = c("log2FoldChange", "pvalue") ##更改列名

## 把SYMBOL转换为ENTREZID，可能有部分丢失
gene <- bitr(rownames(nrsce.markers), fromType = "SYMBOL", toType =  "ENTREZID", OrgDb = org.Hs.eg.db)

## 基因名、ENTREZID、logFC一一对应起来
gene$logFC = nrsce.markers$log2FoldChange[match(gene$SYMBOL, rownames(nrsce.markers))]

## 构建genelist
geneList = gene$logFC
names(geneList) = gene$ENTREZID 
geneList = sort(geneList, decreasing = T) # 降序，按照logFC的值来排序

## GSEA分析
# 准备基因集
downgene = read.csv(file = "./4-macsub/gsea_genelist.csv", header = F)
downgene = downgene$V1
downgene = str_to_upper(downgene)
downgene = bitr(downgene, fromType = "SYMBOL", toType =  "ENTREZID", OrgDb = org.Hs.eg.db)
downgene = downgene$ENTREZID

gene_set = data.frame(
  gs_name = "cGASgene",
  entrez_gene = downgene
)

# df = sce.markers[downgene,]


# GSEA
gsea_result = GSEA(geneList = geneList, 
                   TERM2GENE = gene_set, 
                   minGSSize = 1,
                   maxGSSize = 1000, 
                   pvalueCutoff = 1,
                   verbose = FALSE)

gsea_result_table = gsea_result@result
write.csv(gsea_result_table, file = "./5-episub/episub.gsea.csv")

# 可视化
gseaplot2(gsea_result, geneSetID = "cGASgene", pvalue_table = T)


geneSetID = c("cGASgene")
gseaNb(object = gsea_result,
       geneSetID = geneSetID,
       subPlot = 2,
       termWidth = 35,
       legend.position = c(0.8,0.8),
       addPval = T,
       pvalX = 1, pvalY = 1,
       curveCol = jjAnno::useMyCol('paired', 5))

