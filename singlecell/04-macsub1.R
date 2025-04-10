
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
macsub = subset(x = sce.all.pd1.harmony, celltype == "Macrophage")
table(macsub@meta.data$orig.ident)

## 创建Seurat数据
counts = macsub@assays$RNA@layers$counts
rownames(counts) = rownames(macsub)
colnames(counts) = colnames(macsub)

meta.data = macsub@meta.data
macsub = CreateSeuratObject(counts = counts, meta.data = meta.data)
rm(list = c("counts", "meta.data", "sce.all.pd1.harmony"))

## 添加分组信息
table(macsub@meta.data$orig.ident)
table(macsub@meta.data$group)

## 标准化，归一化，PCA
macsub = NormalizeData(macsub, normalization.method =  "LogNormalize", scale.factor = 1e4)
macsub = FindVariableFeatures(macsub, selection.method = "vst", nfeatures = 2000) 
macsub = ScaleData(macsub) 
macsub = RunPCA(object = macsub, features = VariableFeatures(macsub))

## umap
macsub = RunUMAP(macsub, dims = 1:20, reduction = "pca")
DimPlot(macsub, reduction = "umap", group.by = "orig.ident", label = F, cols = col_vec)
DimPlot(macsub, reduction = "umap", group.by = "group", label = F, cols = col_vec)

save(macsub, file = "./4-macsub/macsub.Rdata")

# 计算marker基因
Idents(macsub) = macsub@meta.data$group
sce.markers = FindMarkers(object = macsub, ident.1 = "pCR", ident.2 = "non-pCR", only.pos = FALSE, min.pct = 0, thresh.use = 0)
write.csv(sce.markers, file = "./4-macsub/macsub.degs.csv")



