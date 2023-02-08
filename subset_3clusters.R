suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))

root_dir='a.Analysis'

processed_seurat = readRDS(paste0(root_dir,'/data/seu_processed.RDS'))

meta_sub = processed_seurat@meta.data[(processed_seurat@meta.data$clusters_medres %in% c(4,1,10)) & (processed_seurat@meta.data$sample %in% c('SN_ABCP2d','SN_ABCP3d')),]

subset_seu = processed_seurat[, rownames(meta_sub)]


hvg = VariableFeatures(subset_seu)

subset_seu = RunPCA(subset_seu, features = hvg)

subset_seu <- RunUMAP(subset_seu, dims = seq_len(30))


saveRDS(subset_seu, 'seu_subset_3clusters.rds')

