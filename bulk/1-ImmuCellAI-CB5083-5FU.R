
# 准备环境
rm(list = ls())
library(ggpubr)
library(ImmuCellAImouse)
library(RColorBrewer)
library(ggsignif)
source(file = "./0-RawData/1-CB5083-5FU/count_tpm_fpkm.R")


# 加载数据
exp = fread(input = "./0-RawData/1-CB5083-5FU/cb5083_5fu_fpkm_expression.csv", data.table = F)
exp = exp[!duplicated(exp$gene_name), ]
rownames(exp) = exp$gene_name
exp = exp[,-1]

# fpkm转tpm
tpms = apply(exp, 2, fpkmToTpm)
colSums(tpms)
exp[1:3, 1:4]
tpms[1:3, 1:4]
tpms = as.data.frame(tpms)

# 免疫浸润分析
tme = ImmuCellAI_mouse(sample = exp,
                       data_type = "rnaseq",  #"rnaseq"/"microarray"
                       group_tag = 0,  #是否有分组信息, 如果没有则填"0"
                       customer=FALSE)

tme = tme$abundance
write.csv(tme, file = "./1-CB5083-5FU/5-immucellai.csv")

# 可视化
pdata = as.data.frame(tme[,-37])
pdata$group = factor(x = c(rep("A", 3), rep("B", 3), rep("C", 3), rep("D", 3)), levels = c("A", "B", "C", "D"))
pdata = pivot_longer(pdata, cols = 1:36, names_to = "celltype", values_to = "proportion")


ggplot(pdata, aes(celltype, proportion, fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  facet_wrap(.~celltype, scales = "free", ncol = 9)+
  scale_fill_aaas()+
  geom_signif(comparisons = list(c("A", "D")))

ggsave(filename = "./1-CB5083-5FU/4-immucellai.pdf", width = 16, height = 10)














