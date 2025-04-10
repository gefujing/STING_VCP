
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
library(ggsignif)
library(GseaVis)

load("./4-macsub/macsub.Rdata")
load(file = "./1-qc/mycol.Rdata")
sce.markers = read.csv(file = "./4-macsub/macsub.degs.csv", header = T)
rownames(sce.markers) = sce.markers$X

# 功能聚类
features = c("OAS1", "TMEM173", "IRF3", "IRF5", "CXCL9", "CXCL10", "CCL5", "CCL13", "STAT1", "IFNAR1", "TLR4", "IFI27")
VlnPlot(macsub, features = features, group.by = "group", cols = col_vec, ncol = 6, pt.size = 0)
ggsave(filename = "./4-macsub/1-sce.markerss-vln.pdf", width = 12, height = 6.7) 

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
write.csv(kkdiff, file = "./4-macsub/macsub.kegg.csv")

pathways = c("Th17 cell differentiation", "T cell receptor signaling pathway", "Th1 and Th2 cell differentiation",
             "Toll-like receptor signaling pathway", "Fc gamma R-mediated phagocytosis", "Chemokine signaling pathway", 
             "IL-17 signaling pathway", "NOD-like receptor signaling pathway", "Leukocyte transendothelial migration",
             "Cytosolic DNA-sensing pathway")
dotplot(kk.diff, showCategory = pathways, color = "pvalue")


pdata = kkdiff[kkdiff$Description %in% pathways, ]
pdata = pdata[order(pdata$Count, decreasing = T),]
pdata$order = factor(rev(as.integer(1:10)),
                     labels = rev(pdata$Description))

ggplot(pdata, aes(x = Count, y = order, fill = pvalue)) + 
  geom_bar(stat='identity') +
  scale_fill_gradient(low = "red",high = "blue" ) +
  labs(x = "Count", 
       y = "Pathways")+
  theme_classic()+
  theme(text = element_text(colour = "black"))
ggsave(filename = "./4-macsub/2-sce.markers-KEGG1.pdf", width = 5.5, height = 3.5)


ggplot(pdata, aes(x = Count, y = order, fill = order)) + 
  geom_bar(stat='identity') +
  scale_fill_manual(values = col_vec) +
  labs(x = "Count", 
       y = "Pathways")+
  theme_classic()+
  theme(text = element_text(colour = "black"))
ggsave(filename = "./4-macsub/2-sce.markers-KEGG2.pdf", width = 7.5, height = 3.5)


## KEGG
gene_diff = unique(gene_diff)
go.diff = enrichGO(gene = gene_diff, OrgDb = "org.Hs.eg.db", ont = "ALL", pvalueCutoff = 0.9, qvalueCutoff = 0.9, readable=T)
godiff = go.diff@result
write.csv(kkdiff, file = "./4-macsub/macsub.go.csv")

pathways = c("immune response-activating signaling pathway", "mononuclear cell differentiation", "positive regulation of cytokine production",
             "chemotaxis", "regulation of T cell activation", "leukocyte mediated immunity", 
             "T cell differentiation", "lymphocyte mediated immunity", "mononuclear cell proliferation",
             "cell chemotaxis")

pdata = godiff[godiff$Description %in% pathways, ]
pdata = pdata[order(pdata$Count, decreasing = T),]
pdata$order = factor(rev(as.integer(1:10)),
                     labels = rev(pdata$Description))

ggplot(pdata, aes(x = Count, y = order, fill = pvalue)) + 
  geom_bar(stat='identity') +
  scale_fill_gradient(low = "red",high = "blue" ) +
  labs(x = "Count", 
       y = "Pathways")+
  theme_classic()+
  theme(text = element_text(colour = "black"))
ggsave(filename = "./4-macsub/2-sce.markers-GO1.pdf", width = 6, height = 3.5)


ggplot(pdata, aes(x = Count, y = order, fill = order)) + 
  geom_bar(stat='identity') +
  scale_fill_manual(values = col_vec) +
  labs(x = "Count", 
       y = "Pathways")+
  theme_classic()+
  theme(text = element_text(colour = "black"))
ggsave(filename = "./4-macsub/2-sce.markers-GO2.pdf", width = 8.5, height = 3.5)


##GSEA
## 上一步差异分析得到差异基因列表sce.markers后取出，p值和log2FC
nrsce.markers = sce.markers[ ,c("avg_log2FC", "p_val")]
colnames(nrsce.markers) = c("log2FoldChange", "pvalue") ##更改列名

## 把SYMBOL转换为ENTREZID，可能有部分丢失
gene = bitr(rownames(nrsce.markers), fromType = "SYMBOL", toType =  "ENTREZID", OrgDb = org.Hs.eg.db)

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
write.csv(gsea_result_table, file = "./4-macsub/macsub.gsea.csv")

# 可视化
gseaplot2(gsea_result, geneSetID = "cGASgene", pvalue_table = T)

gseaNb(gsea_result, )

geneSetID = c("cGASgene")
gseaNb(object = gsea_result,
       geneSetID = geneSetID,
       subPlot = 2,
       termWidth = 35,
       legend.position = c(0.8,0.8),
       addPval = T,
       pvalX = 1, pvalY = 1,
       curveCol = jjAnno::useMyCol('paired', 5))
















