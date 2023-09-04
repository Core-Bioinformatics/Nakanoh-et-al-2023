suppressMessages(library(Seurat))


raw_seu =readRDS('saved_object/filtered_feature_bc_matrix_seu.rds')

meta = read.csv('MFE56636-meta_fixed_index.csv',row.names=1)

raw_seu = subset(raw_seu, cells=rownames(meta))


meta$orig.ident = raw_seu@meta.data$orig.ident

meta$nCount_RNA = raw_seu@meta.data$nCount_RNA

meta$nFeature_RNA = raw_seu@meta.data$nFeature_RNA

raw_seu@meta.data = meta


############################


# remove MT and RP genes
all.index = 1:nrow(raw_seu)
MT.index <- grep(pattern = "^MT", x = rownames(raw_seu), value = FALSE)
RP.index = grep(pattern = "^RP[SL][[:digit:]]", x = rownames(raw_seu), value = FALSE)


raw_seu= raw_seu[!((all.index %in% MT.index) | (all.index %in% RP.index)  ), ]


###########################

variable.features.max = 3000  

filtered_seu <- SCTransform(raw_seu, verbose = FALSE, return.only.var.genes=FALSE,
                   variable.features.n = variable.features.max)

filtered_seu <- RunPCA(filtered_seu, 
                       verbose = FALSE, approx = FALSE)
filtered_seu <- RunUMAP(filtered_seu, dims = 1:30, verbose = FALSE)

##########################

umap_df = meta[,c('UMAP_1','UMAP_2')]


colnames(umap_df) = c('og_UMAP1','og_UMAP2')

umap_df$UMAP_1 = filtered_seu[["umap"]]@cell.embeddings[,1]

umap_df$UMAP_2 = filtered_seu[["umap"]]@cell.embeddings[,2]


filtered_seu[["umap"]]@cell.embeddings = as.matrix(umap_df)

##########################

saveRDS(filtered_seu,'saved_object/Zhai2022_seu_filtered_56636cells_SCT_ogUMAP.rds')














