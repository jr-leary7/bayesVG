library(Seurat)
seu_brain <- SeuratData::LoadData("stxBrain", type = "anterior1")
gene_set_1 <- Matrix::rowSums(seu_brain@assays$Spatial$counts > 0) >= 10L
gene_set_2 <- Matrix::rowMeans(seu_brain@assays$Spatial$counts) >= 0.1
genes_keep <- rownames(seu_brain)[gene_set_1 & gene_set_2]
seu_brain <- subset(seu_brain, features = genes_keep)
spot_set_1 <- seu_brain$nCount_Spatial >= 500L
spot_set_2 <- seu_brain$nFeature_Spatial >= 1000L
spots_keep <- colnames(seu_brain)[spot_set_1 & spot_set_2]
seu_brain <- subset(seu_brain, cells = spots_keep)
seu_brain <- PercentageFeatureSet(seu_brain, 
                                  pattern = "^mt-", 
                                  col.name = "percent_mito")
seu_brain <- subset(seu_brain, subset = percent_mito < 20)
save(seu_brain, file = "../data/seu_brain.rda", compress = "xz")
