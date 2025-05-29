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
hvg_metadata <- getBayesianGeneStats(seu_pbmc)
hvg_fit <- extractModel(seu_pbmc)
hvg_plot <- plotHVGs(seu_pbmc)

# load spatial data & preprocess 
load(system.file("data/seu_brain.rda", package = "bayesVG"))
seu_brain <- suppressWarnings({ 
  Seurat::NormalizeData(seu_brain, verbose = FALSE) %>% 
  Seurat::FindVariableFeatures(nfeatures = 1000L, verbose = FALSE)
})

# fit each kernel to spatial coordinates matrix 
spatial_mtx <- scale(as.matrix(dplyr::select(Seurat::GetTissueCoordinates(seu_brain), -cell)))
M <- nrow(spatial_mtx)
k <- 20
kmeans_centers <- stats::kmeans(spatial_mtx, centers = k, iter.max = 100L)$centers
dists_centers <- as.matrix(stats::dist(kmeans_centers))
lscale <- stats::median(dists_centers[upper.tri(dists_centers)])
phi_exp_quad <- phi_matern <- phi_periodic <- matrix(0, nrow = M, ncol = k)
for (i in seq(k)) {
  d2 <- rowSums((spatial_mtx - matrix(kmeans_centers[i, ], nrow = M, ncol = 2, byrow = TRUE))^2)
  phi_exp_quad[, i] <- expQuadKernel(d2, length.scale = lscale)
  phi_matern[, i] <- maternKernel(d2, length.scale = lscale, nu = 2.5)
  phi_periodic[, i] <- periodicKernel(d2, length.scale = lscale, period = 100L)
}

# fit spatial model, extract output, cluster SVGs, & run enrichment on SVG modules
seu_brain <- findSpatiallyVariableFeaturesBayes(seu_brain, 
                                                naive.hvgs = Seurat::VariableFeatures(seu_brain),
                                                lscale.estimator = "variogram", 
                                                kernel = "matern", 
                                                kernel.smoothness = 1.5, 
                                                n.cores = 1L, 
                                                save.model = TRUE) %>% 
             classifySVGs(n.SVG = 300L)
seu_brain <- suppressWarnings({
  Seurat::ScaleData(seu_brain, verbose = FALSE) %>% 
  Seurat::RunPCA(features = Seurat::VariableFeatures(seu_brain), 
                 npcs = 20L, 
                 verbose = FALSE, 
                 seed.use = 312) %>% 
  Seurat::FindNeighbors(reduction = "pca", 
                        dims = 1:20, 
                        k.param = 20L, 
                        nn.method = "annoy",
                        annoy.metric = "cosine", 
                        verbose = FALSE) %>% 
  Seurat::FindClusters(resolution = 0.5,
                       random.seed = 312, 
                       verbose = FALSE) %>% 
  Seurat::RunUMAP(reduction = "pca", 
                  dims = 1:20, 
                  n.neighbors = 20L, 
                  n.components = 2L, 
                  metric = "cosine", 
                  seed.use = 312, 
                  verbose = FALSE)
})
svg_metadata <- getBayesianGeneStats(seu_brain)
svg_fit <- extractModel(seu_brain)
svg_plot <- plotSVGs(seu_brain)
svg_clusters <- clusterSVGsBayes(seu_brain,
                                 svgs = Seurat::VariableFeatures(seu_brain), 
                                 n.clust = 3L, 
                                 n.cores = 1L)
seu_brain <- scoreSpatialModules(seu_brain, 
                                 svg.clusters = svg_clusters, 
                                 n.cores = 1L)
enrich_res <- enrichSpatialModules(svg_clusters, species = "mmusculus")


# compute naive gene statistics
gene_stats_naive <- computeNaiveGeneStatistics(seu_pbmc, use.norm = TRUE)

# convert seu_brain to SpatialExperiment from Seurat
spe_brain <- suppressWarnings(convertToSpatialExperiment(seu_brain, sample.id = "anterior1"))

# run downstream plotting utilities 
p1 <- plotSpatialExpression(seu_brain, gene.plot = "Nrgn")
p2 <- plotSpatialExpression(spe_brain, gene.plot = "Nrgn")
p3 <- plotTissueImage(seu_brain)
p4 <- plotSpatialAttributes(seu_brain, attribute.plot = "seurat_clusters")
p5 <- plotModuleScores(seu_brain, 
                       module.plot = "1", 
                       plot.type = "spatial")
p6 <- plotModuleScores(seu_brain, 
                       module.plot = "1", 
                       plot.type = "embedding", 
                       embedding.name = "umap")
p7 <- plotModuleScores(seu_brain, 
                       module.plot = "1", 
                       plot.type = "violin", 
                       violin.group = "seurat_clusters")

# run HVG tests 
test_that("HVG model", {
  expect_s4_class(seu_pbmc, "Seurat")
  expect_s3_class(hvg_metadata, "data.frame")
  expect_equal(ncol(hvg_metadata), 20)
  expect_equal(nrow(hvg_metadata), 100)
  expect_s3_class(hvg_fit, "brmsfit")
  expect_s3_class(hvg_plot, "ggplot")
})

# run kernel tests 
test_that("kernels", {
  expect_type(phi_exp_quad, "double")
  expect_equal(ncol(phi_exp_quad), 20)
  expect_equal(nrow(phi_exp_quad), ncol(seu_brain))
  expect_type(phi_matern, "double")
  expect_equal(ncol(phi_matern), 20)
  expect_equal(nrow(phi_matern), ncol(seu_brain))
  expect_type(phi_periodic, "double")
  expect_equal(ncol(phi_periodic), 20)
  expect_equal(nrow(phi_periodic), ncol(seu_brain))
})

# run SVG tests
test_that("SVG model", {
  expect_s4_class(seu_brain, "Seurat")
  expect_s3_class(svg_metadata, "data.frame")
  expect_equal(ncol(svg_metadata), 18)
  expect_equal(nrow(svg_metadata), 1000)
  expect_s3_class(svg_fit, "CmdStanVB")
  expect_s3_class(svg_plot, "ggplot")
  expect_type(svg_clusters, "list")
  expect_s3_class(svg_clusters$cluster_df, "data.frame")
  expect_type(svg_clusters$pca_embedding, "double")
  expect_s3_class(svg_clusters$model_fit, "CmdStanVB")
  expect_type(svg_clusters$log_likelihood, "double")
  expect_type(svg_clusters$BIC, "double")
  expect_type(seu_brain$svg_cluster_1_UCell, "double")
  expect_type(seu_brain$svg_cluster_2_UCell, "double")
  expect_type(seu_brain$svg_cluster_3_UCell, "double")
  expect_s3_class(enrich_res, "data.frame")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
  expect_s3_class(p5, "ggplot")
  expect_s3_class(p6, "ggplot")
  expect_s3_class(p7, "ggplot")
})

# run naive gene statistics tests 
test_that("naive gene statistics", {
  expect_s3_class(gene_stats_naive, "data.frame")
  expect_equal(ncol(gene_stats_naive), 4)
  expect_equal(nrow(gene_stats_naive), 100)
})

# run spatialexperiment conversion tests 
test_that("SpatialExperiment conversion", {
  expect_s4_class(spe_brain, "SpatialExperiment")
  expect_equal(ncol(spe_brain), 2444)
  expect_equal(nrow(spe_brain), 11464)
})
