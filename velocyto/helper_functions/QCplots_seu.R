QCplots_seu <- function(seu, plots.pdf = NULL, group.by = "orig.ident", pdf_width = 7){
  library(Seurat)
  
  categories <- grep("RNA$", colnames(seu@meta.data), value = TRUE)
  
  if(!is.null(plots.pdf)){
    pdf(plots.pdf, width = pdf_width)
  }
  
  print(VlnPlot(seu, categories, ncol = length(categories), pt.size = 0, group.by = group.by))
  print(VlnPlot(seu, categories, ncol = length(categories), pt.size = 0, log = TRUE, group.by = group.by))
  print(FeatureScatter(seu, "nCount_RNA", "nFeature_RNA", group.by = group.by))
  print(FeatureScatter(seu, "nCount_RNA", "MT_RNA", group.by = group.by))
  print(FeatureScatter(seu, "nCount_RNA", "RB_RNA", group.by = group.by))
  
  if(!is.null(plots.pdf)){
    dev.off()
  }
  
  invisible(seu)
}
