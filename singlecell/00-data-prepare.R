
# 设置环境
rm(list = ls()) 
options(stringsAsFactors = F) 

library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(readr)
library(stringr)
library(ggsci)

# 批量读取数据
## 设置数据路径与样本名称
assays = dir("./0-RawData/")
dir = paste0("./0-RawData/", assays)

## 按文件顺序给样本命名，名称不要以数字开头，中间不能有空格 
samples_name = assays

## 批量读取
scRNAlist = list()
for(i in 1:length(dir)){
  counts = Read10X(data.dir = dir[i])
  scRNAlist[[i]] = CreateSeuratObject(counts, project = samples_name[i], min.cells = 10, min.features = 500)
  scRNAlist[[i]] = RenameCells(scRNAlist[[i]], add.cell.id = samples_name[i]) 
  if(T){scRNAlist[[i]][["percent.mt"]] = PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")}
  if(T){scRNAlist[[i]][["percent.rb"]] = PercentageFeatureSet(scRNAlist[[i]], pattern = "^RP[SL]")}
  if(T){HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  HB.genes = CaseMatch(HB.genes, rownames(scRNAlist[[i]]))
  scRNAlist[[i]][["percent.HB"]]<-PercentageFeatureSet(scRNAlist[[i]], features=HB.genes)}
}

names(scRNAlist) = samples_name

## 合并
sce.all = merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
sce.all

## 添加metadata
sce.all@meta.data$treatment = ifelse(str_detect(string = sce.all@meta.data$orig.ident, "PD1_C"), "PD1_C", 
                                     ifelse(str_detect(string = sce.all@meta.data$orig.ident, "PD1"), "PD1", "Untreated"))
sce.all@meta.data$treatment = factor(sce.all@meta.data$treatment, levels = c("Untreated", "PD1", "PD1_C"))

table(sce.all@meta.data$treatment)
save(sce.all, file = "./1-qc/sce.all.Rdata")

## 配色
col_vec = c(pal_npg()(8),  pal_aaas()(8), pal_nejm()(8), pal_lancet()(8), pal_jama()(7), pal_bmj()(8), pal_jco()(8))[sample(55,55)]
col_vec = unique(col_vec)
save(col_vec, file = "./1-qc/mycol.Rdata")

