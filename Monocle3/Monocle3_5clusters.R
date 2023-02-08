suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(SeuratWrappers))
suppressMessages(library(monocle3))

seu = readRDS('seu_subset_5clusters.rds')

seu@reductions$umap = seu@reductions$pca

seu@reductions$umap@cell.embeddings = seu@reductions$umap@cell.embeddings[,1:3]

seu@meta.data$monocle3_clusters = seu@meta.data$clusters_medres
seu@meta.data$monocle3_partitions = as.factor(rep(1,dim(seu@meta.data)[1]))

cds <- as.cell_data_set(seu, default.reduction='umap')

data_cds <- learn_graph(cds, use_partition = TRUE, 
                        close_loop =FALSE, verbose=TRUE)



threshold90 = read.csv('meta_gene_expression_threshold90.csv',row.names=1)

roots_df = threshold90[which(threshold90$threshold90 %in% c(4)),]

order_data_cds <- order_cells(data_cds, 
                        reduction_method = "UMAP",
                       root_cells=rownames(roots_df))


meta = seu@meta.data

meta$monocle_pseudotime = order_data_cds@principal_graph_aux@listData[['UMAP']][['pseudotime']]

write.csv(meta,'meta_5clusters_monocle_pseudotime.csv')

