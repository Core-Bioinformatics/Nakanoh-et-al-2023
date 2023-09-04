suppressMessages(library(Seurat))

dir_data='/path/to/data'

seurat = readRDS(paste0(dir_data,'/saved_object/Zhai2022_seu_filtered_56636cells_SCT.rds'))

write.csv(seurat[["pca"]]@cell.embeddings,'saved_object/pca_56636cells_50dim_zhai2022.csv')



