# load scRNA-seq data & subset 
load(system.file("data/seu_pbmc.rda", package = "bayesVG"))
set.seed(312)
genes_test_sc <- sample(rownames(seu_pbmc), size = 100L)
seu_pbmc <- subset(seu_pbmc, features = genes_test_sc)

# fit scRNA-seq model & extract output 
withr::with_output_sink(tempfile(), {
  seu_pbmc <- findVariableFeaturesBayes(seu_pbmc,
                                        n.cells.subsample = 500L,
                                        algorithm = "meanfield",
                                        n.cores.per.chain = 1L,
                                        save.model = TRUE) %>% 
              classifyHVGs(n.HVG = 50L)
})
hvg_metadata <- seu_pbmc@assays$RNA@meta.data
hvg_fit <- extractModel(seu_pbmc)

# load spatial data & preprocess 
load(system.file("data/seu_brain.rda", package = "bayesVG"))
seu_brain <- Seurat::SCTransform(seu_brain, 
                                 assay = "Spatial", 
                                 variable.features.n = 500L, 
                                 vst.flavor = "v2", 
                                 seed.use = 312, 
                                 verbose = FALSE)

# fit spatial model & extract output 
withr::with_output_sink(tempfile(), {
  seu_brain <- findSpatiallyVariableFeaturesBayes(sp.obj = seu_brain, 
                                                  naive.hvgs = Seurat::VariableFeatures(seu_brain),
                                                  kernel = "matern", 
                                                  kernel.smoothness = 1.5, 
                                                  n.cores = 1L, 
                                                  save.model = TRUE) %>% 
                classifySVGs(n.SVG = 100L)
})
svg_metadata <- seu_brain@assays$SCT@meta.features
svg_fit <- extractModel(seu_brain)

test_that("HVG model", {
  expect_s4_class(seu_pbmc, "Seurat")
  expect_s3_class(hvg_metadata, "data.frame")
  expect_equal(ncol(hvg_metadata), 20)
  expect_equal(nrow(hvg_metadata), 100)
  expect_s3_class(hvg_fit, "brmsfit")
})

test_that("SVG model", {
  expect_s4_class(seu_brain, "Seurat")
  expect_s3_class(svg_metadata, "data.frame")
  expect_equal(ncol(svg_metadata), 10)
  expect_equal(nrow(svg_metadata), 11464)
  expect_s3_class(svg_fit, "CmdStanVB")
})
