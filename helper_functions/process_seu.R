process_seu <- function(
  seu,
  feature_type = c("var", "abn"),
  regress_cell_cycle = TRUE,
  processing.plots = NULL,
  umap.plots = NULL,
  umap.plot.groups = "clusters_umap",
  use.harmony = FALSE,
  harmony.plots = NULL,
  harmony.group.by.vars = "orig.ident",
  harmony.plot.groups = "clusters_harmony",
  save.object = NULL,
  nfeatures = 3000, # change the number of features
  ndims = 25, # PCA components to use
  n_neighbours = 30,
  resolution = 0.8, # adjust clustering resolution to change # of clusters
  raster = FALSE
){
  library(dplyr)
  library(Seurat)
  
  seu <- seu %>%
    NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = nfeatures)
  
  feature_type <- feature_type[1]
  if(feature_type == "var"){
    features <- VariableFeatures(seu)
  }else{
    features <- names(sort(Matrix::rowSums(seu@assays$RNA@counts), decreasing=TRUE)[seq_len(nfeatures)])
  }
  
  if(regress_cell_cycle){
    seu <- CellCycleScoring(seu, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
    vars.to.regress <- c("S.Score", "G2M.Score")
  }else{
    vars.to.regress <- NULL
  }
  
  seu <- seu %>%
    ScaleData(features = features, vars.to.regress = vars.to.regress) %>%
    RunPCA(features = features) %>%
    FindNeighbors(dims = seq_len(ndims),
                  k.param = n_neighbours) %>%
    FindClusters(resolution = resolution) %>%
    RunUMAP(dims = seq_len(ndims))
  seu[["clusters_umap"]] <- seu[["seurat_clusters"]]
  
  if(!is.null(processing.plots)){
    pdf(processing.plots)
    print(VariableFeaturePlot(seu))
    print(DimPlot(seu, reduction = "pca"))
    print(ElbowPlot(seu, ndims = min(ndims * 2, 50)))
    dev.off()
  }
  
  if(!is.null(umap.plots)){
    pdf(umap.plots)
    for(grp in umap.plot.groups){
      print(DimPlot(seu, reduction = "umap", group.by = grp, raster = raster) +
              guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5))))
    }
    dev.off()
  }
  
  if(use.harmony){
    seu <- seu %>%
      harmony::RunHarmony(group.by.vars = harmony.group.by.vars) %>%
      FindNeighbors(reduction = "harmony", dims = seq_len(ndims)) %>%
      FindClusters(reduction = "harmony", resolution = resolution) %>%
      RunUMAP(reduction = "harmony", dims = seq_len(ndims))
    seu[["clusters_harmony"]] <- seu[["seurat_clusters"]]
    
    if(!is.null(harmony.plots)){
      pdf(harmony.plots)
      for(grp in harmony.plot.groups){
        print(DimPlot(seu, reduction = "harmony", group.by = grp, raster = raster) +
                guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5))))
      }
      dev.off()
    }
  }
  
  if(!is.null(save.object)){
    saveRDS(seu, save.object)
  }
  
  invisible(seu)
}
