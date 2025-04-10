
# import packages
rm(list = ls())
library(edgeR)
library(GenomicAlignments)
library(Rsubread)
library(Biostrings)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(GO.db)
library(tidyverse)

# read counts files
con1 = read.csv(file = "./0-RawData/2-CRISPR/slc19a1.con1.csv")
con1$UID = paste0(con1$UID1, "_", con1$UID2)
rownames(con1) = con1$UID

con2 = read.csv(file = "./0-RawData/2-CRISPR/slc19a1.con2.csv")
con2$UID = paste0(con2$UID1, "_", con2$UID2)
rownames(con2) = con2$UID

cgamp1 = read.csv(file = "./0-RawData/2-CRISPR/slc19a1.cgamp1.csv")
cgamp1$UID = paste0(cgamp1$UID1, "_", cgamp1$UID2)
rownames(cgamp1) = cgamp1$UID

cgamp2 = read.csv(file = "./0-RawData/2-CRISPR/slc19a1.cgamp2.csv")
cgamp2$UID = paste0(cgamp2$UID1, "_", cgamp2$UID2)
rownames(cgamp2) = cgamp2$UID

# contrust deglist object
uidall = union(con1$UID, union(con2$UID, union(cgamp1$UID, cgamp2$UID)))

con1 = con1[uidall,]
con2 = con2[uidall,]
cgamp1 = cgamp1[uidall,]
cgamp2 = cgamp2[uidall,]

# counts
counts = data.frame(row.names = uidall,
                    con1 = con1$counts,
                    con2 = con2$counts,
                    cgamp1 = cgamp1$counts,
                    cgamp2 = cgamp2$counts)

counts[is.na(counts)] = 0

# samples
group = factor(c("con", "con", "cgamp", "cgamp"), levels = c("con", "cgamp"))

samples = data.frame(group = group,
                     sampleName = colnames(counts),
                     biorep = c(1,1,2,2))

genes = con1[,c(1,2,6)]

# object
d = DGEList(counts = counts, samples = samples, genes = genes)
d

# Remove control guides
d.raw = d
d = d[!d$genes$Symbol %in% grep("non", d$genes$Symbol, value = TRUE), ]
d = d[!d$genes$Symbol %in% grep("safe", d$genes$Symbol, value = TRUE), ]

# Permissive filtering
th1 = as.vector(cpm(1, median(d$samples$lib.size)))
th1 

th2 = 2 
keep.exprs = rowSums(cpm(d) > th1) >= th2
table(keep.exprs)

d.filtered.p = d[keep.exprs, , keep.lib.sizes = FALSE]
genes_sgrna = d.filtered.p$genes$Symbol
length(unique(genes_sgrna))

table(table(genes_sgrna))

# TMM normalised counts
yy = d.filtered.p[, order(d.filtered.p$samples$group, colSums(d.filtered.p$counts))]
yy$samples$norm.factors <- normLibSizes(yy$counts + 100, method = "TMM")
yy$samples[, 1:3]
yy.norm = yy  # we continue downstream analysis with this object

# Create the design matrix
design = model.matrix(~group, data = yy.norm$samples)
design

# Fit a model
yy.norm = estimateDisp(yy.norm, design)
bcv = sqrt(yy.norm$common.dispersion)

fit = glmFit(yy.norm, design)
lrt = glmLRT(fit, coef = "groupcgamp")
lrt2.toxa = topTags(lrt, n = Inf, sort.by = "PValue")$table

de = lrt2.toxa[complete.cases(lrt2.toxa), ]

colnames(de)[colnames(de) == "logFC"] = "log2FoldChange"
colnames(de)[colnames(de) == "PValue"] = "pvalue"
colnames(de)[colnames(de) == "Symbol"] = "gene_symbol"
colnames(de)[colnames(de) == "FDR"] = "adjusted_pvalue"
de$diffabundant = "Pass p-value cuttoff"
write.csv(de, file = "./2-CRISPR/1-seq-degs-slc19a1.csv")

# mean
de2 = aggregate(de$log2FoldChange, by = list(de$gene_symbol), FUN = mean)

# Differential abundance analysis at gene level
genesymbols = yy.norm$genes[, 2]
genesymbollist = list()

unq = unique(genesymbols)
unq = unq[!is.na(unq)]

for (i in unq) {
  sel = genesymbols == i & !is.na(genesymbols)
  if (sum(sel) > 3)
    genesymbollist[[i]] = which(sel)
}

fry.res = fry(yy.norm, index = genesymbollist, design, contrast = "groupcgamp")
degs = fry.res

degs$Symbol = rownames(degs)
degs$logfc = de2$x[match(rownames(degs), de2$Group.1)]
write.csv(degs, file = "./2-CRISPR/2-gene-degs-slc19a1.csv")


# vocanol
pdata = degs
for_label = pdata[c("TMEM173", "IRF3", "TBK1", "LRRC8A", "VCP"),]

ggplot(data = pdata) + 
  geom_point(aes(x = logfc, y = -log10(PValue.Mixed), color = logfc, size = -log10(PValue.Mixed))) + 
  geom_text_repel(data =  for_label, aes(x = logfc, y = -log10(PValue.Mixed), label = Symbol),
                  nudge_x = 0.5, nudge_y = 0.2, segment.curvature = -0.1,
                  segment.ncp = 3, direction = "y", hjust = "left", max.overlaps = 200 ) + 
  scale_color_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                        values = seq(0, 1, 0.2)) +
  scale_fill_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                       values = seq(0, 1, 0.2)) +
  geom_vline(xintercept = c(log2(1)), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 4) + 
  theme_bw() + 
  ggtitle(label = "Volcano Plot")


# rank plot
degs = read.csv(file = "./2-CRISPR/2-gene-degs-slc19a1.csv")
rownames(degs) = degs$X
degs = degs[,-1]


pdata = degs[,c(5,7,8)]
pdata = pdata[pdata$PValue.Mixed < 0.05, ]
pdata = pdata[order(pdata$logfc, decreasing = T),]
pdata$rank = 1:nrow(pdata)

kk = c("SLC19A1", "TBK1", "VCP", "VPS4B")

lab_gene = pdata[kk,]
lab_gene$Discription = c(rep("STING", 2), rep("MY", 2))

p = ggplot(data = pdata, aes(x = rank, y = logfc)) +
  geom_point(size = 1) +
  geom_hline(yintercept = 3, lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -3, lty = 4, col = "black", lwd = 0.8) +
  theme_bw()
p


p1 = p + 
  geom_point(data = lab_gene, size = 2, aes(color = Discription)) +
  geom_text_repel(aes(label = Symbol), data = lab_gene, 
                  max.overlaps = Inf, color = "black",
                  nudge_x = 2000, direction="y", hjust = 0,
                  fontface = "italic", size = 3)+
  scale_color_manual(values = c("#3b9a9c", "#dd0a35"))

p1
ggsave(p1, filename = "./2-CRISPR/1-slc19a1-rank.pdf", width = 3, height = 3.5)


