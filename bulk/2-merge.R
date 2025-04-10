

# 设置环境
rm(list = ls())
library(eulerr)

# 准备数据
crispr.scl19a1 = read.csv(file = "./2-CRISPR/2-gene-degs-slc19a1.csv")
crispr.scl19a1 = crispr.scl19a1[(crispr.scl19a1$PValue.Mixed < 0.05), ]
crispr.scl19a1 = crispr.scl19a1[order(crispr.scl19a1$logfc),]
crispr.scl19a1 = crispr.scl19a1[1:300,]
crispr.scl19a1 = crispr.scl19a1$Symbol

crispr.lrrc8a = read.csv(file = "./2-CRISPR/2-gene-degs-lrrc8a.csv")
crispr.lrrc8a = crispr.lrrc8a[(crispr.lrrc8a$PValue.Mixed < 0.05), ]
crispr.lrrc8a = crispr.lrrc8a[order(crispr.lrrc8a$logfc),]
crispr.lrrc8a = crispr.lrrc8a[1:300,]
crispr.lrrc8a = crispr.lrrc8a$Symbol

mergene = intersect(crispr.scl19a1, crispr.lrrc8a)

write.csv(mergene, file = "./2-CRISPR/3-merge-genes.csv")

# plot
dat = c("crispr.scl19a1" = 212, 
        "crispr.lrrc8a" = 212, 
        "crispr.scl19a1&crispr.lrrc8a" = 88)

plot(euler(dat),
     quantities = c(212,212,88))


