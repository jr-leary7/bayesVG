# load scRNA-seq data & subset
load(system.file("data/seu_pbmc.rda", package = "bayesVG"))
set.seed(312)
genes_test_sc <- sample(rownames(seu_pbmc), size = 100L)
seu_pbmc <- subset(seu_pbmc, features = genes_test_sc)
sce_pbmc <- suppressWarnings(Seurat::as.SingleCellExperiment(seu_pbmc))

# fit scRNA-seq model & extract output
seu_pbmc <- findVariableFeaturesBayes(seu_pbmc,
                                      n.cells.subsample = 500L,
                                      algorithm = "meanfield",
                                      n.cores.per.chain = 1L,
                                      save.model = TRUE) %>%
            classifyHVGs(n.HVG = 50L)
hvg_metadata <- getBayesianGeneStats(seu_pbmc)
hvg_fit_seu <- extractModel(seu_pbmc)
hvg_plot <- plotHVGs(seu_pbmc)
sce_pbmc <- findVariableFeaturesBayes(sce_pbmc,
                                      n.cells.subsample = 500L,
                                      algorithm = "meanfield",
                                      n.cores.per.chain = 1L,
                                      save.model = TRUE) %>%
            classifyHVGs(n.HVG = 50L)
hvg_fit_sce <- extractModel(sce_pbmc)

# load spatial data & preprocess + convert to SpatialExperiment
load(system.file("data/seu_brain.rda", package = "bayesVG"))
seu_brain <- suppressWarnings({
  Seurat::NormalizeData(seu_brain, verbose = FALSE)
})
seu_brain_nb <- seu_brain

# convert seu_brain to SpatialExperiment from Seurat
spe_brain <- suppressWarnings(
  convertToSpatialExperiment(seu_brain, 
                             sample.id = "anterior1", 
                             scale.coords = TRUE)
)

# fit each kernel to spatial coordinates matrix
spatial_mtx <- coop::scaler(as.matrix(dplyr::select(Seurat::GetTissueCoordinates(seu_brain), -cell)))
M <- nrow(spatial_mtx)
k <- 20
kmeans_centers <- stats::kmeans(spatial_mtx, centers = k, iter.max = 100L, nstart = 10L)$centers
dists_centers <- fields::rdist(kmeans_centers)
lscale <- stats::median(dists_centers[upper.tri(dists_centers)])
phi_exp_quad <- phi_matern <- phi_periodic <- matrix(0, nrow = M, ncol = k)
for (i in seq(k)) {
  d2 <- rowSums((spatial_mtx - matrix(kmeans_centers[i, ], nrow = M, ncol = 2, byrow = TRUE))^2)
  phi_exp_quad[, i] <- expQuadKernel(d2, length.scale = lscale)
  phi_matern[, i] <- maternKernel(d2, length.scale = lscale, nu = 2.5)
  phi_periodic[, i] <- periodicKernel(d2, length.scale = lscale, period = 100L)
}
phi_exp_quad <- qr.Q(qr(phi_exp_quad, LAPACK = TRUE))
phi_matern <- qr.Q(qr(phi_matern, LAPACK = TRUE))
phi_periodic <- qr.Q(qr(phi_periodic, LAPACK = TRUE))

# fit spatial model, extract output, cluster SVGs, & run enrichment on SVG modules
naive_hvgs_seu <- getNaiveHVGs(seu_brain, n.hvg = 750L)
naive_hvgs_spe <- getNaiveHVGs(spe_brain, n.hvg = 750L)
seu_brain <- findSpatiallyVariableFeaturesBayes(seu_brain,
                                                naive.hvgs = naive_hvgs_seu,
                                                lscale.estimator = "kmeans",
                                                kernel = "matern",
                                                kernel.smoothness = 1.5,
                                                n.cores = 1L,
                                                save.model = TRUE) %>%
             classifySVGs(n.SVG = 300L)
seu_brain_nb <- findSpatiallyVariableFeaturesBayes(seu_brain_nb,
                                                   naive.hvgs = naive_hvgs_seu,
                                                   likelihood = "nb",
                                                   lscale.estimator = "kmeans",
                                                   kernel = "matern",
                                                   kernel.smoothness = 1.5,
                                                   n.cores = 1L,
                                                   save.model = TRUE) %>%
                classifySVGs(n.SVG = 300L)
spe_brain <- findSpatiallyVariableFeaturesBayes(spe_brain,
                                                naive.hvgs = naive_hvgs_spe,
                                                lscale.estimator = "kmeans",
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
svg_metadata_seu <- getBayesianGeneStats(seu_brain)
svg_metadata_seu_nb <- getBayesianGeneStats(seu_brain_nb)
svg_metadata_spe <- getBayesianGeneStats(spe_brain)
svg_fit_seu <- extractModel(seu_brain)
svg_fit_seu_nb <- extractModel(seu_brain_nb)
svg_fit_spe <- extractModel(spe_brain)
svg_plot <- plotSVGs(seu_brain)
svg_clusters <- clusterSVGsBayes(seu_brain,
                                 svgs = Seurat::VariableFeatures(seu_brain),
                                 n.clust = 2L,
                                 n.cores = 1L)
seu_brain <- scoreSpatialModules(seu_brain,
                                 svg.clusters = svg_clusters,
                                 n.cores = 1L)
enrich_res <- enrichSpatialModules(svg_clusters, species = "mmusculus")

# compute naive gene statistics
gene_stats_naive <- computeNaiveGeneStatistics(seu_pbmc, use.norm = TRUE)

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
  expect_s4_class(sce_pbmc, "SingleCellExperiment")
  expect_equal(ncol(hvg_metadata), 20)
  expect_equal(nrow(hvg_metadata), 100)
  expect_s3_class(hvg_fit_seu, "brmsfit")
  expect_s3_class(hvg_fit_sce, "brmsfit")
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
  expect_s4_class(seu_brain_nb, "Seurat")
  expect_s4_class(spe_brain, "SpatialExperiment")
  expect_s3_class(svg_metadata_seu, "data.frame")
  expect_s3_class(svg_metadata_seu_nb, "data.frame")
  expect_s3_class(svg_metadata_spe, "data.frame")
  expect_equal(ncol(svg_metadata_seu), 12)
  expect_equal(nrow(svg_metadata_seu), 750)
  expect_equal(ncol(svg_metadata_seu_nb), 12)
  expect_equal(nrow(svg_metadata_seu_nb), 750)
  expect_equal(ncol(svg_metadata_spe), 10)
  expect_equal(nrow(svg_metadata_spe), 750)
  expect_s3_class(svg_fit_seu, "CmdStanVB")
  expect_s3_class(svg_fit_seu_nb, "CmdStanVB")
  expect_s3_class(svg_fit_spe, "CmdStanVB")
  expect_s3_class(svg_plot, "ggplot")
  expect_type(svg_clusters, "list")
  expect_s3_class(svg_clusters$cluster_df, "data.frame")
  expect_type(svg_clusters$pca_embedding, "double")
  expect_s3_class(svg_clusters$model_fit, "CmdStanVB")
  expect_type(svg_clusters$log_likelihood, "double")
  expect_type(svg_clusters$BIC, "double")
  expect_type(seu_brain$svg_cluster_1_UCell, "double")
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
  expect_type(naive_hvgs_seu, "character")
  expect_length(naive_hvgs_seu, 750)
  expect_type(naive_hvgs_spe, "character")
  expect_length(naive_hvgs_spe, 750)
})

# run spatialexperiment conversion tests
test_that("SpatialExperiment conversion", {
  expect_s4_class(spe_brain, "SpatialExperiment")
  expect_equal(ncol(spe_brain), 2444)
  expect_equal(nrow(spe_brain), 9369)
})
