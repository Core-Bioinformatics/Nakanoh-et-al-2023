suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

root_dir='a.Analysis'
processed_seurat = readRDS(paste0(root_dir,'/data/seu_processed.RDS'))

sub = processed_seurat@meta.data$clusters_medres %in% c(0,2,9,12,13)
subset_seu = processed_seurat[, sub]


hvg = VariableFeatures(subset_seu)


subset_seu = RunPCA(subset_seu, features = hvg)

subset_seu <- RunUMAP(subset_seu, dims = seq_len(30))


saveRDS(subset_seu, 'seu_subset_5clusters.rds')

