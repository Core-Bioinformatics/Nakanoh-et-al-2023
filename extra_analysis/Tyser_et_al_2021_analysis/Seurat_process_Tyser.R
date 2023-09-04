suppressMessages(library(ggplot2))
suppressMessages(library(Seurat))

dir_data='human-gastrula-shiny'
raw_counts = readRDS(paste0(dir_data,'/raw_reads.rds'))
meta = readRDS(paste0(dir_data,'/umap.rds'))

rownames(raw_counts) = meta$cell_name
t_raw_counts = t(raw_counts)
head(t_raw_counts)

rownames(meta) = meta$cell_name

#############################

raw_seurat <- CreateSeuratObject(counts = t_raw_counts, project = "tyser",meta.data =meta,min.cells=3)

# remove MT and RP genes
all.index = 1:nrow(raw_seurat)
MT.index <- grep(pattern = "^MT", x = rownames(raw_seurat), value = FALSE)
RP.index = grep(pattern = "^RP[SL][[:digit:]]", x = rownames(raw_seurat), value = FALSE)


raw_seurat = raw_seurat[!((all.index %in% MT.index) | (all.index %in% RP.index)  ), ]
raw_seurat

#############################


variable.features.max = 4000  ### from the paper they used 4000

filtered_seu <- SCTransform(raw_seurat, verbose = FALSE, return.only.var.genes=FALSE,
                   variable.features.n = variable.features.max)

filtered_seu <- RunPCA(filtered_seu, 
                       verbose = FALSE, approx = FALSE)
filtered_seu <- RunUMAP(filtered_seu, dims = 1:30, verbose = FALSE)

#############################

umap_df = meta[,c('X0','X1')]

colnames(umap_df) = c('og_UMAP1','og_UMAP2')
head(umap_df,3)

umap_df$UMAP_1 = filtered_seu[["umap"]]@cell.embeddings[,1]

umap_df$UMAP_2 = filtered_seu[["umap"]]@cell.embeddings[,2]

filtered_seu[["umap"]]@cell.embeddings = as.matrix(umap_df)

new_columns=c('orig.ident','nCount_RNA','nFeature_RNA',
              'cluster_id',
              'sub_cluster','cell_name','nCount_SCT','nFeature_SCT')

filtered_seu@meta.data = filtered_seu@meta.data[,new_columns]


saveRDS(filtered_seu,'tyer_seu_filtered_1195cells_SCT_ogUMAP.rds')



















