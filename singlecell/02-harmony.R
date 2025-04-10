
# 准备环境
rm(list = ls()) 
options(stringsAsFactors = F) 

library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(stringr)
library(data.table)
library(harmony)

load(file = "./1-qc/sce.all.pd1.Rdata")
load(file = "./1-qc/mycol.Rdata")

# 标准化数据
sce.all.pd1 = NormalizeData(sce.all.pd1, normalization.method = "LogNormalize", scale.factor = 1e4) 
sce.all.pd1 = FindVariableFeatures(sce.all.pd1)
sce.all.pd1 = ScaleData(sce.all.pd1)
sce.all.pd1 = RunPCA(sce.all.pd1, features = VariableFeatures(object = sce.all.pd1))
sce.all.pd1 = RunHarmony(sce.all.pd1, group.by.vars = "orig.ident")
sce.all.pd1 = RunUMAP(sce.all.pd1, dims = 1:15, reduction = "harmony")

DimPlot(sce.all.pd1, reduction = "umap", group.by = "orig.ident", label = F, cols = col_vec)
DimPlot(sce.all.pd1, reduction = "umap", group.by = "group", label = F, cols = col_vec)

sce.all.pd1.harmony = sce.all.pd1

# 分群
sce.all.pd1.harmony = FindNeighbors(sce.all.pd1.harmony, reduction = "harmony", dims = 1:20) 
sce.all.pd1.harmony = FindClusters(sce.all.pd1.harmony, resolution = 1, algorithm = 1)

# 可视化
DimPlot(sce.all.pd1.harmony, reduction = "umap", group.by = "RNA_snn_res.1", cols = col_vec)

save(sce.all.pd1.harmony, file = "./2-harmony/sce.all.pd1.harmony.Rdata")
