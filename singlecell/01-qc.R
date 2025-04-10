
## 设置环境
rm(list = ls()) 
options(stringsAsFactors = F) 
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(stringr)
library(data.table)
library(ggsci)
load(file = "./1-qc/sce.all.Rdata")
load(file = "./1-qc/mycol.Rdata")

## 计算线粒体基因比例
mito_genes = rownames(sce.all)[grep("^MT-", rownames(sce.all))] # 人和鼠的基因名字稍微不一样 
mito_genes #13个线粒体基因
sce.all = PercentageFeatureSet(sce.all, "^MT-", col.name = "percent_mito")

## 计算核糖体基因比例
ribo_genes = rownames(sce.all)[grep("^RP[SL]", rownames(sce.all))]
ribo_genes
sce.all = PercentageFeatureSet(sce.all, "^RP[SL]", col.name = "percent_ribo")

## 计算红血细胞基因比例
rownames(sce.all)[grep("^HB[^(P)]", rownames(sce.all))]
sce.all = PercentageFeatureSet(sce.all, "^HB[^(P)]", col.name = "percent_hb")

## 可视化细胞的上述比例情况
feats = c("nFeature_RNA", "nCount_RNA")
VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2, cols = col_vec) + NoLegend()

feats = c("percent_mito", "percent_ribo", "percent_hb")
VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3, same.y.lims = T, cols = col_vec) + NoLegend()


# 根据上述指标，过滤低质量细胞/基因
## 过滤指标1:最少表达基因数的细胞&最少表达细胞数的基因
sce.all = JoinLayers(sce.all)
selected_c = WhichCells(sce.all, expression = (nFeature_RNA > 500) & (nFeature_RNA < 8000) & (nCount_RNA < 50000))
selected_f = rownames(sce.all)[Matrix::rowSums(sce.all@assays$RNA@layers$counts > 0) > 200]
sce.all.filt = subset(sce.all, features = selected_f, cells = selected_c)

dim(sce.all) 
dim(sce.all.filt) 

table(sce.all@meta.data$orig.ident) 
table(sce.all.filt@meta.data$orig.ident) 

## 过滤指标2:线粒体/核糖体基因比例(根据上面的violin图)
selected_mito = WhichCells(sce.all.filt, expression = percent_mito < 50)
selected_ribo = WhichCells(sce.all.filt, expression = percent_ribo > 3)
selected_hb = WhichCells(sce.all.filt, expression = percent_hb < 1)

length(selected_mito)
length(selected_ribo)
length(selected_hb)

sce.all.filt = subset(sce.all.filt, cells = selected_mito)
sce.all.filt = subset(sce.all.filt, cells = selected_ribo)
sce.all.filt = subset(sce.all.filt, cells = selected_hb)
dim(sce.all.filt)

table(sce.all.filt$orig.ident) 

## 过滤特定基因
sce.all.filt = sce.all.filt[!grepl("^MT-", rownames(sce.all.filt), ignore.case = T), ]
sce.all.filt = sce.all.filt[!grepl("^RP[SL]", rownames(sce.all.filt), ignore.case = T), ]
dim(sce.all.filt) 


## 可视化过滤后的情况
feats = c("nFeature_RNA", "nCount_RNA")
VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2, cols = col_vec) + NoLegend()

feats = c("percent_mito", "percent_ribo", "percent_hb")
VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3, cols = col_vec) + NoLegend()

## 添加meta
sce.all.filt@meta.data$group = ifelse(sce.all.filt@meta.data$orig.ident %in% c("P12T_PD1", "P18T_PD1", "P31T_PD1", "P31T_U", "P26T_PD1_C"),
                                      "non-pCR", "pCR")

sce.all.filt@meta.data$group = factor(sce.all.filt@meta.data$group, levels = c("non-pCR", "pCR"))

## 保存数据
save(sce.all.filt, file = "./1-qc/sce.all.filt.Rdata")

## subset
sce.all.pd1 = subset(sce.all.filt, treatment == "PD1")
save(sce.all.pd1, file = "./1-qc/sce.all.pd1.Rdata")








