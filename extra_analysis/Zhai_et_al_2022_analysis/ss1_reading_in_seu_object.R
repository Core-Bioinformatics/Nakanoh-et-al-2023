suppressMessages(library(Seurat))

bc_expression <- Read10X(data.dir = 'filtered_feature_bc_matrix')


bc_seu = CreateSeuratObject(counts =bc_expression,min.cells=3)


saveRDS(bc_seu,'saved_object/filtered_feature_bc_matrix_seu.rds')




























