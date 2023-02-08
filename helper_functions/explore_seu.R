explore_seu <- function(
    seu, 
    shiny.dir = NULL,
    cloupe.dir = NULL,
    cloupe.prefix = "",
    cloupe.proj = "umap",
    cloupe.categories = "seurat_clusters",
    cellphonedb.path = NULL,
    cellphonedb.cell_type = "seurat_clusters",
    marker.csv = NULL,
    overwrite.marker.csv = FALSE,
    marker.categories = "seurat_clusters",
    presence.plots.pdf = NULL, 
    markers.per.cluster.to.plot = Inf,
    verbose = TRUE
){
  invisible(seu)
  
  if(!is.null(shiny.dir)){
    if(verbose) message("Creating shiny app in ", shiny.dir)
    scConf <- ShinyCell::createConfig(seu)
    ShinyCell::makeShinyApp(seu, scConf, shiny.dir = shiny.dir)
  }
  
  if(!is.null(cloupe.dir)){
    if(verbose) message("Creating cloupe csvs in ", cloupe.dir)
    source(paste0(substr(getwd(), 1, gregexpr("userName", getwd())[[1]][1] + 5), "cloupe_custom_csvs.R"))
    if(!is.null(cloupe.proj)){
      export.cloupe.projection(
        seu, 
        paste0(cloupe.dir, "/", cloupe.prefix, "cloupe_", cloupe.proj, "_projection.csv"), 
        cloupe.proj
      )
    }
    if(!is.null(cloupe.categories)){
      export.cloupe.category(
        seu, 
        paste0(cloupe.dir, "/", cloupe.prefix, "cloupe_", paste(cloupe.categories, collapse="_"), "_category.csv"), 
        cloupe.categories
      )
    }
  }
  
  if(!is.null(cellphonedb.path)){
    if(verbose) message("Creating cellphonedb metadata in ", cellphonedb.path)
    cellphonedb_meta <- tibble::tibble(
      Cell = colnames(seu),
      cell_type = seu@meta.data[, cellphonedb.cell_type]
    )
    write.table(cellphonedb_meta, file = cellphonedb.path, 
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  
  if(!is.null(marker.csv)){
    if(overwrite.marker.csv | !file.exists(marker.csv)){
      if(verbose) message("Generating markers and saving in ", marker.csv)
      Idents(seu) <- marker.categories
      markers <- FindAllMarkers(seu)
      write.csv(x = markers, file = marker.csv, row.names = FALSE)
    }else{
      if(verbose) message("Found existing marker table ", marker.csv, ", reading it")
      markers <- read.csv(marker.csv)
    }
    
    if(!is.null(cloupe.dir)){
      export.cloupe.markers(markers, paste0(cloupe.dir, "/cloupe_features.csv"))
    }
    
    if(!is.null(presence.plots.pdf)){
      if(verbose) message("Creating presence plots in ", presence.plots.pdf)
      markers <- markers %>% 
        dplyr::group_by(cluster) %>% 
        dplyr::slice_head(n = markers.per.cluster.to.plot)
      
      df.all <- cbind(
        seu@reductions$umap@cell.embeddings,
        seu@meta.data, 
        t(as.matrix(seu@assays$RNA@data[match(markers$gene, rownames(seu)),]))
      ) %>%
        tidyr::pivot_longer(cols = markers$gene, names_to = "gene", values_to = "exp") %>%
        dplyr::left_join(markers, by = "gene")
      
      pdf(presence.plots.pdf)
      for(i in seq_len(nrow(markers))){
        df <- dplyr::filter(df.all, gene == markers$gene[i]) %>% dplyr::arrange(exp)
        p <- ggplot(df) +
          theme_minimal() +
          geom_point(mapping = aes(x = UMAP_1, y = UMAP_2, colour = exp), size = 1) +
          scale_colour_gradientn(colours = c('#E6E6E6', RColorBrewer::brewer.pal(9, "YlOrRd"))) +
          labs(title=paste0(markers$gene[i], " (cluster ", markers$cluster[i], ")"))
        print(p)
      }
      dev.off()
    }
  }
  
  if(verbose) message("Done!")
  invisible(seu)
}
