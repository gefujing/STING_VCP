
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
library(scRNAtoolVis)

load("./4-macsub/macsub.Rdata")
load(file = "./1-qc/mycol.Rdata")

# 加载数据
exp = AverageExpression(object = macsub, assays = "RNA", group.by = "orig.ident")
exp = as.data.frame(exp[["RNA"]])
exp = exp[apply(exp, 1, function(x) sum(x > 0) > 5), ]
exp = exp[,c(1,4,9,2,3,5,6,7,8)]

# 分组
Group = c("none-pCR", "none-pCR", "none-pCR", "pCR", "pCR", "pCR", "pCR", "pCR", "pCR")
Group = factor(x = Group, levels = c("none-pCR", "pCR"))

# 差异基因
cg = read.csv(file = "./4-macsub/macsub.degs.csv", header = T)
cg = cg$X
cg = intersect(cg, rownames(exp))

mat = exp[cg, ]
mat = as.matrix(mat)
mat = t(scale(t(mat)))

geneseq = as.data.frame(mat[,c(4:9)])
geneseq2 = rownames(geneseq)[order(rowSums(geneseq), decreasing = T)]
mat = mat[geneseq2, ]


# mat = as.data.frame(mat)
# mat$gene = rownames(mat)
# mat$seq = 1:nrow(mat)

write.csv(mat, file = "./4-macsub/macsub.heatmap.matrix.csv")

ha = HeatmapAnnotation(Group = Group,
                       col = list(Group = c("none-pCR" = "#00AFBB", "pCR" = "#E7B800")),
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
genes = read.csv("./4-macsub/heatmap_genelist.csv", header = F)
genes = genes$V1
genes = str_to_upper(genes)
genes = genes[c(2,5,6,9,11)]
genes = as.data.frame(genes)
colnames(genes) = "genes"

p + rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% genes$genes),
                                   labels = genes$genes, labels_gp = gpar(fontsize = 10)))


VlnPlot(macsub, features = genes$genes, group.by = "group", cols = col_vec, ncol = 6, pt.size = 0)
ggsave(filename = "./4-macsub/1-sce.markerss-vln.pdf", width = 9, height = 3) 


averageHeatmap(object = genes, markerGene = genes$genes, annoCol = T, myanCol = col_vector[1:10]) 





