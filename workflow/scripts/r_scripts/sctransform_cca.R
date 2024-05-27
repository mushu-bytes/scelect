library(SeuratObject)
library(Seurat)
library(remotes)
library(sctransform)
library(glmGamPoi)
library(ggplot2)

path <- "/data/immune.rds"
data <- readRDS(path)
data <- SplitObject(data, split.by = "development_stage")

data <- list(data$`31-year-old human stage`,
             data$`30-year-old human stage`)
normalized_data <- lapply(X = data, FUN = SCTransform, method = "glmGamPoi")

features <- SelectIntegrationFeatures(object.list = normalized_data,
                                      nfeatures = 3000)
pre_integration <- PrepSCTIntegration(object.list = normalized_data,
                                      anchor.features = features)

# npcs must be less than the number of cells in the assay
pre_integration <- lapply(X = pre_integration,
                          FUN = RunPCA,
                          features = features)

# Very difficult to integrate across datasets with varying cells
anchors <- FindIntegrationAnchors(object.list = pre_integration,
                                  normalization.method = "SCT",
                                  anchor.features = features,
                                  dims = 1:30,
                                  reduction = "cca",
                                  k.anchor = 20)

ccat_sct <- IntegrateData(anchorset = anchors,
                          normalization.method = "SCT",
                          dims = 1:30)
ccat_sct <- RunPCA(ccat_sct)
resolution_range <- seq(from = 0, to = 1, by = 0.1)
ccat_sct <- FindNeighbors(ccat_sct, dims = 1:30)
ccat_sct  <- RunUMAP(ccat_sct, reduction = "pca", dims = 1:30)
ccat_sct <- FindClusters(ccat_sct, resolution = resolution_range)

# Adjust the contrast in the plot
dim_plot <- DimPlot(object = ccat_sct, reduction = "umap")

ggsave(filename = "/plots/cca_umap_plot.png", plot = dim_plot)