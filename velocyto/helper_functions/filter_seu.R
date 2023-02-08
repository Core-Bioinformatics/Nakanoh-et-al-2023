filter_seu <- function(
  seu, 
  organism = "mmusculus",
  nCount_min = 0,
  nCount_max = Inf,
  nFeature_min = 0,
  nFeature_max = Inf,
  MT_max = 100,
  RB_max = 100,
  save.object = NULL,
  rb_genes_mm10_path = "data/list_ribo_su_gene_mouse.txt"
){
  library(Seurat)
  
  ncells.before <- ncol(seu)
  
  seu <- subset(seu, subset = 
                  nCount_RNA > nCount_min & 
                  nCount_RNA < nCount_max &
                  nFeature_RNA > nFeature_min & 
                  nFeature_RNA < nFeature_max &
                  MT_RNA < MT_max &
                  RB_RNA < RB_max)
  
  if(organism == "hsapiens"){
    genes.to.exclude <- c(grep("^MT-", rownames(seu), value = TRUE),
                          grep("^RP[SL]", rownames(seu), value = TRUE))
  }else if(organism == "mmusculus"){
    rb_genes_mm10 <- read.delim(file = rb_genes_mm10_path, header = FALSE)[, 1]
    genes.to.exclude <- c(grep("^mt-", rownames(seu), value = TRUE), rb_genes_mm10)
  }else{
    genes.to.exclude <- NULL
  }
  seu <- seu[setdiff(rownames(seu), genes.to.exclude), ]
  
  ncells.after <- ncol(seu)
  message("Filtered seurat object; ",
          "before: ", ncells.before, " cells; ", 
          "after: ", ncells.after, " cells ", 
          "(", round(ncells.after / ncells.before, 3) * 100, "%)")
  
  if(!is.null(save.object)){
    saveRDS(seu, save.object)
  }
  
  invisible(seu)
}
