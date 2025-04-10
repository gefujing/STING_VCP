
# 设置环境
rm(list = ls()) 
options(stringsAsFactors = F) 

library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(stringr)
library(data.table)
library(ggheatmap)
library(ggsci)

# 设置数据
load(file = "./2-harmony/sce.all.pd1.harmony.Rdata")
load(file = "./1-qc/mycol.Rdata")

# 检测分群相关性
DimPlot(object = sce.all.pd1.harmony, group.by = "seurat_clusters", reduction = "umap", label = T, cols = col_vec)
ggsave(filename = "./3-celltype/1-cluster.umap.pdf", width = 5.5, height = 4.5)

DimPlot(object = sce.all.pd1.harmony, group.by = "orig.ident", reduction = "umap", cols = col_vec) 
DimPlot(object = sce.all.pd1.harmony, group.by = "group", reduction = "umap", cols = col_vec)

# 检查总体基因情况
genes_to_check = c('PTPRC', # immune cell
                   'CD3D', 'CD3E', 'CD4','CD8A', # T Cells 
                   'CD19', 'CD79A', 'MS4A1', # B cells
                   'IGHG1', 'MZB1', 'SDC1', # Plasma cells
                   'CD68', 'CD163', 'CD14', 'C1QA', 'C1QB', 'ITGAM', 'AIF1',# macrophages
                   'TPSAB1', 'TPSB2', # mast cells,
                   'RGS5', 'CD73', 'CD105', 'CD44', # perivascular cellhttp://biotrainee.vip:12133/graphics/plot_zoom_png?width=561&height=787
                   'CD14', 'S100A9', 'S100A8', 'MMP19', # monocyte
                   'FCGR3A', 'FGFBP2', 'CX3CR1', 'KLRB1', 'NCR1', # NK cells
                   'LAMP3', 'IDO1','IDO2',## DC3 
                   'CD1E','CD1C', # DC2
                   'FGF7','MME', 'ACTA2', ## human Fibroblasts 
                   'DCN', 'LUM', 'GSN' , ## mouse PDAC Fibroblasts 
                   'MKI67' , 'TOP2A', 
                   'PECAM1', 'VWF',  'CDH5', ## Endothelial cells
                   'AMY1', 'AMY2A2', 'PRSS1',  ## Acinar cells
                   'EPCAM' , 'KRT19', 'KRT7', 'PROM1', 'ALDH1A1', 'CD24', # epithelial or tumor
                   'CHGB' ## Endocrine cells
)

genes_to_check
DotPlot(sce.all.pd1.harmony, features = unique(genes_to_check), group.by = "seurat_clusters", assay='RNA') + coord_flip()

# 检查总体基因情况
genes_to_check = c('PTPRC', # immune cell
                   'CD3D', 'CD3E', 'CD4','CD8A', # T Cells 
                   'CD79A', 'MS4A1', # B cells
                   'MZB1', 'SDC1', # Plasma cells
                   'CD68', 'CD163', 'CD14', 'C1QA', 'C1QB', 'ITGAM', 'AIF1',# macrophages
                   'TPSAB1', # mast cells,
                   'CD73', 'CD105', 'CD44', # perivascular cellhttp://biotrainee.vip:12133/graphics/plot_zoom_png?width=561&height=787
                   'CD14', 'S100A9', 'S100A8', 'MMP19', # monocyte
                   'FCGR3A', 'KLRB1', # NK cells
                   'LAMP3', 'IDO1','IDO2',## DC3 
                   'CD1E','CD1C', # DC2
                   'ACTA2', ## human Fibroblasts 
                   'DCN', 'LUM', 'GSN' , ## Fibroblasts 
                   'MKI67' , 'TOP2A', 
                   'PECAM1', 'VWF',  'CDH5', ## Endothelial cells
                   'EPCAM' , 'KRT19', 'CD24' # epithelial or tumor
)

genes_to_check
DotPlot(sce.all.pd1.harmony, features = unique(genes_to_check), group.by = "seurat_clusters", assay='RNA') + coord_flip()
ggsave(filename = "./3-celltype/2-cluster.markers2.pdf", width = 8, height = 7)

# 细胞注释
celltype = data.frame(ClusterID = 0:21,
                      celltype = 0:21) 
## 定义细胞亚群
celltype[celltype$ClusterID %in% c(0,6,9,10,13,14,16,20,21), 2] = "Epithelial"
celltype[celltype$ClusterID %in% c(1,2,4,11), 2] = "Tcell"
celltype[celltype$ClusterID %in% c(18), 2] = "ProlifeT"
celltype[celltype$ClusterID %in% c(3), 2] = "Bcell"
celltype[celltype$ClusterID %in% c(12), 2] = "Plasmacell" 
celltype[celltype$ClusterID %in% c(5,17), 2] = "Endothelial"
celltype[celltype$ClusterID %in% c(7,8), 2] = "Macrophage"
celltype[celltype$ClusterID %in% c(15,19), 2] = "Fibroblast"

## 写入细胞亚群
table(celltype$celltype)
sce.all.pd1.harmony@meta.data$celltype = "NA"

for(i in 1:nrow(celltype)){
  sce.all.pd1.harmony@meta.data[which(sce.all.pd1.harmony@meta.data$seurat_clusters == celltype$ClusterID[i]), 'celltype'] <- celltype$celltype[i]}
table(sce.all.pd1.harmony@meta.data$celltype)

## 去除双细胞
Idents(sce.all.pd1.harmony) = "seurat_clusters"
k = sce.all.pd1.harmony@meta.data$seurat_clusters %in% c(22)
k = rownames(sce.all.pd1.harmony@meta.data)[!k]

sce.all.pd1.harmony = subset(sce.all.pd1.harmony, cells = k)
save(sce.all.pd1.harmony, file = "./3-celltype/sce.all.pd1.celltype.Rdata")

# 查看细胞亚群
DimPlot(sce.all.pd1.harmony, reduction = "umap", group.by = "celltype", label = T, cols = col_vec) 
ggsave(filename = "./3-celltype/3-umap.celltype.pdf", width = 5.5, height = 4)

DimPlot(sce.all.pd1.harmony, reduction = "umap", group.by = "celltype", split.by = "group", label = T, cols = col_vec) 
ggsave(filename = "./3-celltype/4-umap.celltype.group.pdf", width = 8, height = 4)


