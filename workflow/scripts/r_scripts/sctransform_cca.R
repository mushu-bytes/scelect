library(SeuratObject)
library(Seurat)
library(remotes)
library(sctransform)
library(glmGamPoi)
library(ggplot2)
library(cluster)

args <- commandArgs(trailingOnly = TRUE)

input_path <- args[[1]]
output_path <- args[[2]]
dataset_key <- args[[3]]

data <- readRDS(input_path)
data <- SplitObject(data, split.by = dataset_key)

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

cca_sct <- IntegrateData(anchorset = anchors,
                         normalization.method = "SCT",
                         dims = 1:30)
cca_sct <- RunPCA(cca_sct)
resolution_range <- seq(from = 0, to = 1, by = 0.1)
cca_sct <- FindNeighbors(cca_sct, dims = 1:30)
cca_sct  <- RunUMAP(cca_sct, reduction = "pca", dims = 1:30)
cca_sct <- FindClusters(cca_sct, resolution = resolution_range)

# Adjust the contrast in the plot
dim_plot <- DimPlot(object = cca_sct, reduction = "umap")

ggsave(filename = glue("{output_path}/cca_umap_plot.png"), plot = dim_plot)