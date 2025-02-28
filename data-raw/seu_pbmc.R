library(Seurat)
seu_pbmc <- SeuratData::LoadData("pbmc3k")
gene_set_1 <- Matrix::rowSums(seu_pbmc@assays$RNA$counts > 0) >= 10L
gene_set_2 <- Matrix::rowMeans(seu_pbmc@assays$RNA$counts) >= 0.1
genes_keep <- rownames(seu_pbmc)[gene_set_1 & gene_set_2]
seu_pbmc <- subset(seu_pbmc, features = genes_keep)
save(seu_pbmc, file = "../data/seu_pbmc.rda", compress = "xz")
