load(system.file("data/seu_pbmc.rda", package = "bayesVG"))
set.seed(312)
genes_test <- sample(rownames(seu_pbmc), size = 100L)
seu_pbmc <- subset(seu_pbmc, features = genes_test)

withr::with_output_sink(tempfile(), {
  seu_pbmc <- findVariableFeaturesBayes(seu_pbmc, 
                                        n.cells.subsample = 500L,
                                        algorithm = "meanfield", 
                                        save.model = TRUE)
  hvg_metadata <- seu_pbmc@assays$RNA@meta.data
  hvg_fit <- extractModel(seu_pbmc)
})

test_that("HVG model", {
  expect_s4_class(seu_pbmc, "Seurat")
  expect_s3_class(hvg_metadata, "data.frame")
  expect_equal(ncol(hvg_metadata), 17)
  expect_equal(nrow(hvg_metadata), 100)
  expect_s3_class(hvg_fit, "brmsfit")
})
