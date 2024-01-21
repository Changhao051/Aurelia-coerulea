
###analysis of similarity across various cell types 
library(Seurat)
library(pheatmap)
setwd("~")
scrna_harmony <- readRDS("Aurelia.rds")
table(scrna_harmony$RNA_snn_res.1.4)
av <- AverageExpression(scrna_harmony,
                        group.by = "RNA_snn_res.1.4",
                        assays = "RNA")
av=av[[1]]
head(av)
#Select the 2000 genes with the highest standard deviation
cg = names(tail(sort(apply(av,1,sd)),2000))
#view the expression matrix of these 1000 genes in each cluster
view(av[cg,])
#View the correlation matrix of cell populations
view(cor(av[cg,],method = "spearman"))
#pheatmap
pheatmap::pheatmap(cor(av[cg,],method ="spearman"))
