
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
exp = fread(input = "./0-RawData/3-BMDM/genes_fpkm_expression.csv", data.table = F)
exp = exp[!duplicated(exp$gene_name), ]
rownames(exp) = exp$gene_name
exp = exp[,-1]
# exp = exp[apply(exp, 1, function(x) sum(x > 0) > 6), ]

# 分组
Group = c(rep("A", 4), rep("B", 4))
Group = factor(x = Group, levels = c("A", "B"))

# PCA
deg = read.csv(file = "./0-RawData/3-BMDM/BVSA_Gene_differential_expression.csv")
deg = deg[!duplicated(deg$gene_name),]

# 加change列, 标记上下调基因
k1 = (deg$pval < 0.05)&(deg$log2.fc. < -1)
k2 = (deg$pval < 0.05)&(deg$log2.fc. > 1)

deg$change = ifelse(test = k1,
                    yes = "down",
                    no = ifelse(test = k2,
                                yes = "up", 
                                no = "stable"))
table(deg$change)

write.csv(x = deg, file = "./3-BMDM/1-degs-bvsa.csv")

# D up 热图
cg = deg$gene_name[deg$significant != "no"]
mat = exp[cg, ]
mat = as.matrix(mat)
mat = t(scale(t(mat)))
mat = na.omit(mat)

geneseq = as.data.frame(mat[,c(5:8)])
geneseq2 = rownames(geneseq)[order(rowSums(geneseq), decreasing = T)]
mat = mat[geneseq2, ]


ha = HeatmapAnnotation(Group = Group,
                       col = list(Group = c("A" = "#00AFBB", "B" = "#E7B800")),
                       show_annotation_name = FALSE)

col_fun = colorRamp2(c(-2, -1, 2), c("#411445", "#411445", "yellow"))

p = Heatmap(mat,
            col = col_fun,
            show_row_names = F,
            cluster_columns = F,
            cluster_rows = F,
            top_annotation = ha)

p

# 基因标注
genes = c("Ccrl2", "Ccl3", "Il6", "Ifnb1", "Ccl4", "Ccl12", "Ifna2", "Ifna4", "Ifna1", "Ifna5")
genes = as.data.frame(genes)
p + rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% genes$genes),
                                   labels = genes$genes, labels_gp = gpar(fontsize = 10)))


