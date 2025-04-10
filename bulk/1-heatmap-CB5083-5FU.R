
# 设置环境
rm(list = ls())
library(stringr)
library(data.table)
library(FactoMineR)
library(factoextra) 
library(dplyr)
library(limma)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)


# 加载数据
exp = fread(input = "./0-RawData/1-CB5083-5FU/cb5083_5fu_fpkm_expression.csv", data.table = F)
exp = exp[!duplicated(exp$gene_name), ]
rownames(exp) = exp$gene_name
exp = exp[,-1]

# 分组
Group = c(rep("A", 3), rep("B", 3), rep("C", 3), rep("D", 3))
Group = factor(x = Group, levels = c("A", "B", "C", "D"))

# 差异基因
cg = read.csv(file = "./0-RawData/1-CB5083-5FU/deglist.csv")
cg = cg$gene_name

mat = exp[cg, ]
mat = as.matrix(mat)
mat = t(scale(t(mat)))

geneseq = as.data.frame(mat[,c(10,11,12)])
geneseq2 = rownames(geneseq)[order(rowSums(geneseq), decreasing = T)]
mat = mat[geneseq2, ]

# mat = as.data.frame(mat)
# mat$gene = rownames(mat)
# mat$seq = 1:nrow(mat)

write.csv(mat, file = "./1-CB5083-5FU/1-heatmap-matrix.csv")

ha = HeatmapAnnotation(Group = Group,
                       col = list(Group = c("A" = "#00AFBB", "B" = "#E7B800", "C" = "#A65E49", "D" = "#8E44AD")),
                       show_annotation_name = FALSE)

col_fun = colorRamp2(c(-1, 3), c("#411445", "yellow"))

p = Heatmap(mat,
            col = col_fun,
            show_row_names = F,
            cluster_columns = F,
            cluster_rows = F,
            top_annotation = ha)

p

# 基因标注
genes = read.csv("./0-RawData/1-CB5083-5FU/ABCDgenelist.csv", header = F)
colnames(genes) = "genes"

p + rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% genes$genes),
                                   labels = genes$genes, labels_gp = gpar(fontsize = 10)))






















