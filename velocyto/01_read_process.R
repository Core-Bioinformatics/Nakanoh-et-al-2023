library(Seurat)
library(dplyr)
library(ggplot2)

source("QCplots_seu.R")
source("filter_seu.R")
source("process_seu.R")
source("run_ClustAssess.R")
source("explore_seu.R")

cellranger_dir <- "3.cellranger/aggregated/"
matrix_path <- "/outs/count/filtered_feature_bc_matrix/"
contents_path <- "../1.RawData/SLX-21814.H7T7CDRX2.s_1.contents.csv"
datadir <- "data/"
outdir <- "outs/"
umap.plot.groups <- c("sample", "clusters_umap")

seu <- Read10X(paste0(cellranger_dir, matrix_path)) %>%
  CreateSeuratObject(project = "SC")
contents <- read.csv(contents_path)

seu@meta.data <- seu@meta.data %>%
  mutate(
    id = substr(colnames(seu), 18, 18),
    sample = factor(contents$Sample.name[as.integer(id)], levels = contents$Sample.name),
    .after = orig.ident
  ) %>%
  mutate(
    MT_RNA = PercentageFeatureSet(seu, pattern = "^MT-")[, 1],
    RB_RNA = PercentageFeatureSet(seu, pattern = "^RP[SL]")[, 1],
  )

saveRDS(seu, file = paste0(datadir, "seu_initial.RDS"))


seu <- seu %>%
  QCplots_seu(plots.pdf = paste0(outdir, "01_QC_initial.pdf"), group.by = "sample", pdf_width = 14) %>%
  filter_seu( 
    organism = "hsapiens",
    nCount_min = 5000,
    nCount_max = 50000,
    nFeature_min = 2500,
    nFeature_max = 7500,
    MT_max = 15,
    RB_max = 35,
    save.object = paste0(datadir, "seu_filtered.RDS")
  ) %>%
  QCplots_seu(plots.pdf = paste0(outdir, "02_QC_filtered.pdf"), group.by = "sample", pdf_width = 14) %>%
  process_seu(
    feature_type = "var",
    processing.plots = paste0(outdir, "03_processing.pdf"),
    umap.plots = paste0(outdir, "04_UMAP.pdf"),
    umap.plot.groups = umap.plot.groups,
    save.object = paste0(datadir, "seu_processed_firstpass.RDS"),
    nfeatures = 2000,
    ndims = 30,
    resolution = 0.8
  ) %>%
  explore_seu(
    shiny.dir = paste0(outdir, "ShinyCell"),
    cloupe.dir = outdir,
    cellphonedb.dir = NULL,
    marker.csv = paste0(outdir, "markers.csv"),
    presence.plots.pdf = NULL,
    verbose = TRUE
  ) %>%
  run_ClustAssess(
    plots = paste0(outdir, "05_ClustAssess"),
    save_objects = paste0(datadir, "ClustAssess_objects.RData"),
    n_features_sequence = seq(from = 500, to = 5000, by = 500)
  )


