
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
Group = c(rep("A", 2), rep("B", 2), rep("C", 2), rep("D", 2))
Group = factor(x = Group, levels = c("A", "B", "C", "D"))

# PCA
dat = as.data.frame(t(exp))

dat.pca = PCA(dat, graph = FALSE)

fviz_pca_ind(dat.pca,
             col.ind = Group, # color by groups
             palette = c("#00AFBB", "#E7B800", "#A65E49", "#8E44AD"),
             addEllipses = T, # Concentration ellipses
             repel = T,
             legend.title = "Groups"
)

ggsave(filename = "./1-CB5083-5FU/1-pca.pdf", width = 4, height = 3)


# TOP100热图
cg = names(tail(sort(apply(exp, 1, sd)), 500))
mat = exp[cg, ]

annotation_col = data.frame(row.names = colnames(mat), Group = Group)
ann_colors = list(Group = c(A = "#00AFBB", B = "#E7B800", C = "#A65E49", D = "#8E44AD"))

pheatmap(mat,
         show_colnames = T,
         show_rownames = F,
         cluster_cols = F,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         scale = "row",
         color = colorRampPalette(colors = c("#411445", "#411445", "yellow"))(100),
         breaks = seq(-1, 1, length.out = 100))


# 取差异分析组
exp.ad = exp
Group.ad = c(rep("ABC", 9), rep("D", 3))
Group.ad = factor(x = Group.ad, levels = c("ABC", "D"))

# 差异分析
design = model.matrix(~Group.ad)
fit = lmFit(exp.ad, design)
fit = eBayes(fit)
deg = topTable(fit, coef = 2, number = Inf)

# 添加行名为一列
deg$symbol = rownames(deg)

# 加change列, 标记上下调基因
k1 = (deg$P.Value < 0.05)&(deg$logFC < -1)
k2 = (deg$P.Value < 0.05)&(deg$logFC > 1)

deg$change = ifelse(test = k1,
                    yes = "down",
                    no = ifelse(test = k2,
                                yes = "up", 
                                no = "stable"))
table(deg$change)

write.csv(x = deg, file = "./1-CB5083-5FU/1-degs-dvsall.csv")


# D up 热图
cg = deg$symbol[deg$change != "stable"]
mat = exp[cg, ]
mat = as.matrix(mat)
mat = t(scale(t(mat)))

ha = HeatmapAnnotation(Group = Group,
                       col = list(Group = c("A" = "#00AFBB", "B" = "#E7B800", "C" = "#A65E49", "D" = "#8E44AD")),
                       show_annotation_name = FALSE)

col_fun = colorRamp2(c(-2, -1, 2), c("#411445", "#411445", "yellow"))

p = Heatmap(mat,
            col = col_fun,
            show_row_names = F,
            cluster_columns = F,
            top_annotation = ha)

p

# 基因标注
genes = cg[1:20]
genes = genes[-6]
genes = as.data.frame(genes)


p + rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% genes$genes),
                                   labels = genes$genes, labels_gp = gpar(fontsize = 10)))
                  

# 差异分析2
exp.ad = exp[,c(1,2,7,8)]
Group.ad = Group[c(1,2,7,8)]
Group.ad = factor(x = Group.ad, levels = c("A", "D"))


design = model.matrix(~Group.ad)
fit = lmFit(exp.ad, design)
fit = eBayes(fit)
deg = topTable(fit, coef = 2, number = Inf)

# 添加行名为一列
deg$symbol = rownames(deg)

# 加change列, 标记上下调基因
k1 = (deg$P.Value < 0.05)&(deg$logFC < -1)
k2 = (deg$P.Value < 0.05)&(deg$logFC > 1)

deg$change = ifelse(test = k1,
                    yes = "down",
                    no = ifelse(test = k2,
                                yes = "up", 
                                no = "stable"))
table(deg$change)

write.csv(x = deg, file = "./1-CB5083-5FU/1-degs-dvsa.csv")


# D up 热图
cg = deg$symbol[deg$change != "stable"]
mat = exp[cg, ]
mat = as.matrix(mat)
mat = t(scale(t(mat)))

geneseq = as.data.frame(mat[,c(7,8)])
geneseq2 = rownames(geneseq)[order(rowSums(geneseq), decreasing = T)]
mat = mat[geneseq2, ]

mat = as.data.frame(mat)
mat$gene = rownames(mat)
mat$seq = 1:nrow(mat)

write.csv(mat, file = "./1-CB5083-5FU/1-heatmap-matrix.csv")

ha = HeatmapAnnotation(Group = Group,
                       col = list(Group = c("A" = "#00AFBB", "B" = "#E7B800", "C" = "#A65E49", "D" = "#8E44AD")),
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
genes = cg[1:20]
genes = genes[-6]
genes = as.data.frame(genes)


p + rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% genes$genes),
                                   labels = genes$genes, labels_gp = gpar(fontsize = 10)))
















# 差异分析3
exp.ac = exp[,c(1,2,3,7,8,9)]
Group.ac = Group[c(1,2,3,7,8,9)]
Group.ac = factor(x = Group.ad, levels = c("A", "C"))


design = model.matrix(~Group.ac)
fit = lmFit(exp.ac, design)
fit = eBayes(fit)
deg = topTable(fit, coef = 2, number = Inf)

# 添加行名为一列
deg$symbol = rownames(deg)

# 加change列, 标记上下调基因
k1 = (deg$P.Value < 0.05)&(deg$logFC < -0.5)
k2 = (deg$P.Value < 0.05)&(deg$logFC > 0.5)

deg$change = ifelse(test = k1,
                    yes = "down",
                    no = ifelse(test = k2,
                                yes = "up", 
                                no = "stable"))
table(deg$change)

write.csv(x = deg, file = "./1-CB5083-5FU/1-degs-cvsa.csv")

# C up 热图
cg = deg$symbol[deg$change != "stable"]
mat = exp[cg, ]
mat = as.matrix(mat)
mat = t(scale(t(mat)))

geneseq = as.data.frame(mat[,c(10,11,12)])
geneseq2 = rownames(geneseq)[order(rowSums(geneseq), decreasing = T)]
mat = mat[geneseq2, ]
write.csv(mat, file = "./1-CB5083-5FU/1-heatmap-matrix.csv")

ha = HeatmapAnnotation(Group = Group,
                       col = list(Group = c("A" = "#00AFBB", "B" = "#E7B800", "C" = "#A65E49", "D" = "#8E44AD")),
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
genes = cg[1:20]
genes = genes[-6]
genes = as.data.frame(genes)


p + rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% genes$genes),
                                   labels = genes$genes, labels_gp = gpar(fontsize = 10)))





