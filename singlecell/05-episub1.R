
# 准备环境
rm(list=ls()) 
options(stringsAsFactors = F) 
library(Seurat)
library(ggplot2)
library(stringr)
library(clustree)
library(cowplot)
library(dplyr)
library(stringr)
library(data.table)
library(patchwork)
library(ggsci)
load("./3-celltype/sce.all.pd1.celltype.Rdata")
load(file = "./1-qc/mycol.Rdata")

# 亚群再分析
table(sce.all.pd1.harmony@meta.data$celltype)
episub = subset(x = sce.all.pd1.harmony, celltype == "Epithelial")
table(episub@meta.data$orig.ident)

## 创建Seurat数据
counts = episub@assays$RNA@layers$counts
rownames(counts) = rownames(episub)
colnames(counts) = colnames(episub)

meta.data = episub@meta.data
episub = CreateSeuratObject(counts = counts, meta.data = meta.data)
rm(list = c("counts", "meta.data", "sce.all.pd1.harmony"))

## 添加分组信息
table(episub@meta.data$orig.ident)
table(episub@meta.data$group)

## 标准化，归一化，PCA
episub = NormalizeData(episub, normalization.method =  "LogNormalize", scale.factor = 1e4)
episub = FindVariableFeatures(episub, selection.method = "vst", nfeatures = 2000) 
episub = ScaleData(episub) 
episub = RunPCA(object = episub, features = VariableFeatures(episub))

## umap
episub = RunUMAP(episub, dims = 1:20, reduction = "pca")
DimPlot(episub, reduction = "umap", group.by = "orig.ident", label = F, cols = col_vec)
DimPlot(episub, reduction = "umap", group.by = "group", label = F, cols = col_vec)

save(episub, file = "./5-episub/episub.Rdata")

# 计算marker基因
Idents(episub) = episub@meta.data$group
sce.markers = FindMarkers(object = episub, ident.1 = "pCR", ident.2 = "non-pCR", only.pos = FALSE, min.pct = 0.01, thresh.use = 0.01)
write.csv(sce.markers, file = "./5-episub/episub.degs.csv")



