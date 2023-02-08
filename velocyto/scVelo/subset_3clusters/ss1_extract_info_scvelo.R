suppressMessages(library(ggplot2))
suppressMessages(library(monocle3))
suppressMessages(library(Seurat))


seu = readRDS('seu_subset_3clusters.rds')

hvg_list = VariableFeatures(seu)
write.csv(hvg_list,'hvg3000_seurat_subset_4_1_10.csv')


pca_df = as.data.frame(Embeddings(seu, reduction = "pca"))
umap_df = as.data.frame(Embeddings(seu,reduction='umap'))

write.csv(pca_df,'pca_df__subset_4_1_10.csv')
write.csv(umap_df,'umap_df__subset_4_1_10.csv')


meta = seu@meta.data


write.csv(meta,'meta__subset_4_1_10.csv')






















