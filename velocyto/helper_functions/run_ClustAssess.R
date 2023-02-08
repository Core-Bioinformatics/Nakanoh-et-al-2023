run_ClustAssess <- function(
    seu,
    plots = NULL,
    save_objects = NULL,
    assay = "RNA",
    n_features_sequence = seq(from = 500, to = 2000, by = 500),
    n_neigh_sequence_conn_comps = c(1:4, seq(from = 5, to = 50, by = 5)),
    n_neigh_sequence_importance = seq(from = 5, to = 50, by = 5),
    resolution_range = seq(from = 0.1, to = 2, by = 0.1),
    additional_gridsearch_ecs_thresholds = c(0.99, 0.95),
    n_repetitions = 100,
    n_cores = 16
){
  library(Seurat)
  library(ClustAssess)
  
  n_features <- max(n_features_sequence)
  
  most_abundant_genes = rownames(seu@assays[[assay]])[
    order(Matrix::rowSums(seu@assays$RNA), decreasing = TRUE)][seq_len(n_features)]
  
  seu <- FindVariableFeatures(seu, nfeatures = n_features)
  var_features = VariableFeatures(seu)
  
  seu <- ScaleData(seu, features = union(var_features, most_abundant_genes))
  
  ma_hv_steps = sapply(n_features_sequence, function(x) { 
    length(intersect(most_abundant_genes[seq_len(x)], var_features))
  })
  ma_hv_genes_intersection <- intersect(most_abundant_genes, var_features)
  
  pca_feature_stability_object = c(
    get_feature_stability(
      data_matrix = seu@assays[[assay]]@scale.data,
      feature_set = most_abundant_genes,
      steps = n_features_sequence,
      n_repetitions = n_repetitions,
      feature_type = "MA",
      graph_reduction_type = "PCA",
      npcs = 30,
      min_dist = 0.3,
      n_neighbors = 30,
      metric = "cosine",
      ncores = n_cores,
      ecs_thresh = 1,
      algorithm = 1
    ),
    get_feature_stability(
      data_matrix = seu@assays[[assay]]@scale.data,
      feature_set = var_features,
      steps = n_features_sequence,
      n_repetitions = n_repetitions,
      feature_type = "HV",
      graph_reduction_type = "PCA",
      npcs = 30,
      min_dist = 0.3,
      n_neighbors = 30,
      metric = "cosine",
      ncores = n_cores,
      ecs_thresh = 1,
      algorithm = 1
    ),
    get_feature_stability(
      data_matrix = seu@assays[[assay]]@scale.data,
      feature_set = ma_hv_genes_intersection,
      steps = ma_hv_steps,
      n_repetitions = n_repetitions,
      feature_type = "MA_HV",
      graph_reduction_type = "PCA",
      npcs = 30,
      min_dist = 0.3,
      n_neighbors = 30,
      metric = "cosine",
      ncores = n_cores,
      ecs_thresh = 1,
      algorithm = 1
    )
  )
  if(!is.null(save_objects)){
    save(list = ls(all.names = TRUE), file = save_objects, envir = environment())
  }
  if(!is.null(plots)){
    pdf(paste0(plots, "_1_feature_stability.pdf"), height = 12, width = 12)
    print(plot_feature_stability_boxplot(pca_feature_stability_object))
    print(plot_feature_stability_ecs_incremental(pca_feature_stability_object))
    print(plot_feature_stability_mb_facet(pca_feature_stability_object))
    print(plot_feature_stability_ecs_facet(pca_feature_stability_object))
    dev.off()
  }
  
  nn_conn_comps_object = c(
    get_nn_conn_comps(
      object = seu@reductions$pca@cell.embeddings,
      n_neigh_sequence = n_neigh_sequence_conn_comps,
      n_repetitions = n_repetitions,
      graph_reduction_type = "UMAP",
      ncores = n_cores,
      min_dist = 0.3,
      n_neighbors = 30,
      metric = "cosine"
    ),
    get_nn_conn_comps(
      object = seu@assays[[assay]]@scale.data,
      n_neigh_sequence = n_neigh_sequence_conn_comps,
      n_repetitions = n_repetitions,
      graph_reduction_type = "PCA",
      ncores = n_cores,
      nv = 30
    )
  )
  if(!is.null(save_objects)){
    save(list = ls(all.names = TRUE), file = save_objects, envir = environment())
  }
  
  nn_importance_object = mapply(
    c,
    get_nn_importance(
      object = seu@assays[[assay]]@scale.data,
      n_neigh_sequence = n_neigh_sequence_importance,
      n_repetitions = n_repetitions,
      graph_reduction_type = "PCA",
      ecs_thresh = 1,
      ncores = n_cores,
      algorithm = 1,
      nv = 30
    ),
    get_nn_importance(
      object = seu@reductions$pca@cell.embeddings,
      n_neigh_sequence = n_neigh_sequence_importance,
      n_repetitions = n_repetitions,
      graph_reduction_type = "UMAP",
      ecs_thresh = 1,
      ncores = n_cores,
      algorithm = 1,
      min_dist = 0.3,
      n_neighbors = 30,
      metric = "cosine"
    ),
    SIMPLIFY = FALSE
  )
  if(!is.null(save_objects)){
    save(list = ls(all.names = TRUE), file = save_objects, envir = environment())
  }
  if(!is.null(plots)){
    pdf(paste0(plots, "_2_graph_construction.pdf"), height = 12, width = 12)
    print(plot_connected_comps_evolution(nn_conn_comps_object))
    print(plot_n_neigh_k_correspondence(nn_importance_object))
    print(plot_n_neigh_ecs(nn_importance_object))
    dev.off()
  }
  
  adj_matrix <- FindNeighbors(seu@reductions$pca@cell.embeddings, k.param = 35)[["snn"]]
  clustering_diff_obj <- get_clustering_difference(
    graph_adjacency_matrix = adj_matrix,
    resolution = resolution_range,
    n_repetitions = n_repetitions,
    ecs_thresh = 1,
    ncores = n_cores,
    algorithm = 1:4
  )
  if(!is.null(save_objects)){
    save(list = ls(all.names = TRUE), file = save_objects, envir = environment())
  }
  if(!is.null(plots)){
    pdf(paste0(plots, "_3_clustering_difference.pdf"), height = 20, width = 14)
    print(plot_clustering_difference_boxplot(clustering_diff_obj))
    print(plot_clustering_difference_facet(
      clustering_diff_obj, 
      seu@reductions$umap@cell.embeddings
    ))
    dev.off()
  }
  
  resolution_gridsearch = get_resolution_importance(
    embedding = seu@reductions$umap@cell.embeddings,
    resolution = resolution_range,
    n_neigh = 30,
    n_repetitions = n_repetitions,
    clustering_method = c(1, 2),
    graph_type = 2,
    ecs_thresh = 1,
    ncores = n_cores
  )
  
  for(ecs_thresh in additional_gridsearch_ecs_thresholds){
    assign(
      paste("resolution_gridsearch_threshold_", ecs_thresh),
      merge_partitions(
        resolution_gridsearch, 
        ecs_thresh = ecs_thresh, 
        ncores = n_cores
      )
    )
  }
  
  if(!is.null(save_objects)){
    save(list = ls(all.names = TRUE), file = save_objects, envir = environment())
  }
  
  if(!is.null(plots)){
    pdf(paste0(plots, "_4_resolution_gridsearch.pdf"), height = 12, width = 12)
    
    print(
      plot_k_resolution_corresp(resolution_gridsearch) + 
        ggtitle("resolution - k correspondence with ecs threshold = 1")
    )
    print(
      plot_k_n_partitions(resolution_gridsearch) + 
        ggtitle("k - # partitions correspondence with ecs threshold = 1")
    )
    
    for(ecs_thresh in additional_gridsearch_ecs_thresholds){
      print(
        plot_k_resolution_corresp(get(paste("resolution_gridsearch_threshold_", ecs_thresh))) + 
          ggtitle(paste("resolution - k correspondence with ecs threshold =", ecs_thresh))
      )
      print(
        plot_k_n_partitions(get(paste("resolution_gridsearch_threshold_", ecs_thresh))) + 
          ggtitle(paste("k - # partitions correspondence with ecs threshold =", ecs_thresh))
      )
    }
    
    dev.off()
  }
  
  invisible(seu)
}