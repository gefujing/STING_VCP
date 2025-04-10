
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
library(pheatmap)
library(ggheatmap)
library(RColorBrewer)
library(ggthemes)
library(patchwork)
library(scRNAtoolVis)

# 设置数据
load(file = "./3-celltype/sce.all.pd1.harmony.Rdata")
load(file = "./1-qc/mycol.Rdata")
table(sce.all.pd1.harmony$celltype)

# 单细胞特征
## 细胞标志基因-文献
genes_to_check = c("EPCAM", "CD3D", "MZB1", "CD79A", "ACTA2", "CD68", "PECAM1", "MKI67")
DotPlot(sce.all.pd1.harmony, features = unique(genes_to_check), group.by = "celltype", assay='RNA')

VlnPlot(sce.all.pd1.harmony, features = unique(genes_to_check), group.by = "celltype", assay='RNA', pt.size = 0, ncol = 4, cols = col_vec)
ggsave(filename = "./3-celltype/5-celltype.markers.pdf", width = 10.5, height = 6)


## 组别比例
## 绘制堆叠条形图(KRAS)
cell.prop = as.data.frame(prop.table(table(sce.all.pd1.harmony$celltype, sce.all.pd1.harmony$group), 2))
colnames(cell.prop) = c("Celltype", "Group", "Proportion")
cell.prop$proportion2 = round(cell.prop$Proportion, 2)

th = theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) 
ggplot(cell.prop, aes(Group, Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  ggtitle("Celltypes of Group") +
  theme_bw() +
  guides(fill = guide_legend(title = NULL)) + 
  th +
  scale_fill_manual(values = col_vec)
  
ggsave(filename = "./3-celltype/6-cell.proportion.pdf", width = 3, height = 3)


