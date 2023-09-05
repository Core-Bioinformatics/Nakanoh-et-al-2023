suppressMessages(library(Seurat))
suppressMessages(library(scmap))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(data.table))

dir_ref='path/to/reference/data/ref.rds'

dir_query='path/to/query/data/query.rds'

num_features=250


ref_sce = readRDS(dir_ref)

query_sce = readRDS(dir_query)

############################################

rowData(ref_sce)$feature_symbol = rownames(ref_sce)
rowData(query_sce)$feature_symbol =rownames(query_sce)

common_gene <- intersect(rownames(ref_sce),rownames(query_sce))
length(common_gene)


ref_common <- ref_sce[common_gene,]
query_common <- query_sce[common_gene,]
identical(rownames(ref_common),rownames(query_common))

###########################################

ref_common <- selectFeatures(ref_common,suppress_plot=FALSE,n_features =num_features)

query_common <- selectFeatures(query_common,suppress_plot=FALSE,n_features =num_features)


##########################################

set.seed(1)

ref_common <- indexCell(ref_common)

scmapCell_results <- scmapCell(query_common,list(reference=metadata(ref_common)$scmap_cell_index))


#########################################

cell_projections <- apply(scmapCell_results$reference$cells, 2, function(x) return(colnames(ref_common)[x]))
cell_projections <- as.data.frame(cell_projections)
                       
                          
nearest_neigbors <- scmapCell_results$reference$cells
nearest_neigbors <- as.data.frame(nearest_neigbors)

                          
cosine_similarities <- scmapCell_results$reference$similarities
cosine_similarities <- as.data.frame(cosine_similarities)                         
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
