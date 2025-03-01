# load scRNA-seq data & subset 
load(system.file("data/seu_pbmc.rda", package = "bayesVG"))
set.seed(312)
genes_test_sc <- sample(rownames(seu_pbmc), size = 100L)
seu_pbmc <- subset(seu_pbmc, features = genes_test_sc)

# fit scRNA-seq model & extract output 
seu_pbmc <- findVariableFeaturesBayes(seu_pbmc,
                                      n.cells.subsample = 500L,
                                      algorithm = "meanfield",
                                      n.cores.per.chain = 1L,
                                      save.model = TRUE) %>% 
            classifyHVGs(n.HVG = 50L)
hvg_metadata <- seu_pbmc@assays$RNA@meta.data
hvg_fit <- extractModel(seu_pbmc)
hvg_plot <- plotHVGs(seu_pbmc)

# load spatial data & preprocess 
load(system.file("data/seu_brain.rda", package = "bayesVG"))
seu_brain <- suppressWarnings({ 
  Seurat::SCTransform(seu_brain, 
                      assay = "Spatial", 
                      variable.features.n = 500L, 
                      vst.flavor = "v2", 
                      seed.use = 312, 
                      verbose = FALSE)
})

# fit spatial model & extract output 
seu_brain <- findSpatiallyVariableFeaturesBayes(seu_brain, 
                                                naive.hvgs = Seurat::VariableFeatures(seu_brain),
                                                kernel = "matern", 
                                                kernel.smoothness = 1.5, 
                                                n.cores = 1L, 
                                                save.model = TRUE) %>% 
             classifySVGs(n.SVG = 100L)
svg_metadata <- seu_brain@assays$SCT@meta.features
svg_fit <- extractModel(seu_brain)
svg_plot <- plotSVGs(seu_brain)

# compute naive gene statistics
gene_stats_naive <- computeNaiveGeneStatistics(seu_pbmc, use.norm = TRUE)

# convert seu_brain to spatialexperiment 
spe_brain <- convertToSpatialExperiment(seu_brain, sample.id = "anterior1")

# run HVG model tests 
test_that("HVG model", {
  expect_s4_class(seu_pbmc, "Seurat")
  expect_s3_class(hvg_metadata, "data.frame")
  expect_equal(ncol(hvg_metadata), 20)
  expect_equal(nrow(hvg_metadata), 100)
  expect_s3_class(hvg_fit, "brmsfit")
  expect_s3_class(hvg_plot, "ggplot")
})

# run SVG model tests
test_that("SVG model", {
  expect_s4_class(seu_brain, "Seurat")
  expect_s3_class(svg_metadata, "data.frame")
  expect_equal(ncol(svg_metadata), 10)
  expect_equal(nrow(svg_metadata), 11464)
  expect_s3_class(svg_fit, "CmdStanVB")
  expect_s3_class(svg_plot, "ggplot")
})

# run naive gene statistics tests 
test_that("naive gene statistics", {
  expect_s3_class(gene_stats_naive, "data.frame")
  expect_equal(ncol(gene_stats_naive), 4)
  expect_equal(nrow(gene_stats_naive), 100)
})

# run spatialexperiment conversion tests 
test_that("spatialexperiment conversion", {
  expect_s4_class(spe_brain, "SpatialExperiment")
  expect_equal(ncol(spe_brain), 2444)
  expect_equal(nrow(spe_brain), 11464)
})
