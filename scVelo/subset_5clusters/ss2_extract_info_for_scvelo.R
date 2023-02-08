suppressMessages(library(Seurat))


seu = readRDS('seu_subset_clusters.rds')

hvg_list = VariableFeatures(seu)

write.csv(hvg_list,'hvg3000_seurat.csv')


pca_df = as.data.frame(Embeddings(seu, reduction = "pca"))
umap_df = as.data.frame(Embeddings(seu,reduction='umap'))

write.csv(pca_df,'pca_df__clusters_0_2_9_12_13.csv')
write.csv(umap_df,'umap_df__clusters_0_2_9_12_13.csv')
