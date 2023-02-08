library(Seurat)
library(dplyr)
library(ggplot2)

source("process_seu.R")
source("explore_seu.R")
source("enrichment.R")

datadir <- "data/"
outdir <- "outs/"

seu <- readRDS(paste0(datadir, "seu_filtered.RDS")) %>%
  process_seu(
    feature_type = "var",
    nfeatures = 3000,
    ndims = 30,
    n_neighbours = 30,
    resolution = 1
  )
head(seu)
seu <- FindClusters(seu, resolution = 0.6)
colnames(seu@meta.data)[match("seurat_clusters", colnames(seu@meta.data))] <- "clusters_lowres"

seu <- FindClusters(seu, resolution = 1.4)
colnames(seu@meta.data)[match("seurat_clusters", colnames(seu@meta.data))] <- "clusters_highres"

seu@meta.data <- seu@meta.data %>%
  rename(clusters_medres = clusters_umap) %>%
  select(1:7, clusters_lowres, clusters_medres, clusters_highres)

saveRDS(seu, paste0(datadir, "seu_processed.RDS"))

explore_seu(seu = seu, shiny.dir = paste0(outdir, "Shota_Jul2022_resolutions"))

markers_enrichment <- function(seu, idents, id, dir_markers, dir_enrichment, organism){
  Idents(seu) <- idents
  message("Finding markers for ", idents)
  markers <- FindAllMarkers(seu)
  write.csv(markers, paste0(dir_markers, "/markers_", id, ".csv"))
  for(cl in levels(Idents(seu))){
    message("Running enrichment for cluster ", cl)
    res <- enrichment(
      genes = markers$gene[markers$cluster == cl],
      background = rownames(seu)[rowSums(seu@assays$RNA@counts[, Idents(seu) == cl]) > 0],
      file = paste0(dir_enrichment, "/enrichment_", id, "_cl_", cl, ".csv"),
      organism = organism
    )
  }
  invisible(seu)
}

for(id in c("lowres", "medres", "highres")){
  markers_enrichment(
    seu = seu,
    idents = paste0("clusters_", id),
    id = id,
    dir_markers = outdir,
    dir_enrichment = paste0(outdir, "/enrichment/", id),
    organism = "hsapiens"
  )
}
