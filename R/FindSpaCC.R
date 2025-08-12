#' @include build_expression_similarity and cosine_similarity functions
#'
NULL

##' Run FindSpaCC clustering algorithm
##'
##' @description Implements the SpaCC algorithm in R using reticulate to run the Python version. Requires the python "leidenalg", "igraph" and "leiden" modules to be installed. Returns a vector of partition indices.
##' Windows users can still this with devtools::install_github("rstudio/reticulate", ref = "86ebb56"); reticulate::use_condaenv("r-reticulate"); reticulate::conda_install("r-reticulate", "leidenalg", channel = "vtraag")
##' Windows users needs to install the FNN package
##' @param optional weights are used to apply "influence" or "confidence" to different nodes on a spatial graph. NULL indicates that all nodes have equal weights.
##' @param k The number of neighbors of each node when constructing a KNN graph. If set to NULL, the function will automatically try multiple k values internally and select the best one. The larger k is, the denser the graph becomes. The smaller the k, the sparser the graph.
##'
##' 1. Internal preprocessing steps
##' The NormalizeData, FindVariableFeatures, and ScaleData functions are used for standardization, selecting highly variable genes, and Z-score scaling to ensure the comparability of expression matrices and reduce dimensions.
##' The highly variable gene matrix was extracted from the normalized expression matrix.
##' build_expression_similarity() is used to calculate the expression similarity between cells and supports cosine/pearson/spearman; Internally, rows of cells are standardized.
##' cosine_similarity() is used to customize the efficient cosine similarity to avoid memory explosion of proxy::simil.
##'
##' 2. Construction of spatial graph
##' k represents the number of KNN neighbors in space. If the user does not specify it, it is fixed at 6.
##' nn_idx represents the nearest k neighbor indexes of each cell, returned by FNN::get.knn(spatial_coords, k).
##' spatial_adj is a sparse binary adjacency matrix
##' spatial_graph, dense 0/1 adjacency matrix, is only used for subsequent alignment with the dimension of the expression matrix.
##'
##' 3. Similarity fusion and threshold
##' combined_sim is used to blend similarity matrices, taking into account both expressive similarity and spatial adjacency, and has been symmetric
##' threshold is a binarized threshold. Take the 70th quantile of the non-zero value of the fusion matrix → retain the first 30% of the strongest edges.
##' adj_bin, the final 0/1 adjacency matrix (combined_sim > threshold)·1, is used for mapping.
##'
##' 4. Community detection
##' g is the igraph object, an undirected weighted graph, where nodes = cells and edge weights = fusion similarity.
##' resolution_parameter is the resolution parameter of the Leiden algorithm. It traverses {0.5, 0.8, 1.0, 1.2, 1.5} and selects the maximum modularity value.
##' best_partition is the optimal community label vector, with length equal to the number of cells and value equal to the community ID.
##' The highest modularization value of best_mod measures the quality of community division. The larger the Q, the tighter the internal structure.
##'
##' 5. Function return value
##' Integrate the SpaCC function
##'
FindSpaCC <- function(seo,
                      spatial_coords,
                      spatial_weight = NULL,
                      k = NULL) {
  seo <- NormalizeData(seo, normalization.method = "LogNormalize")
  seo <- FindVariableFeatures(seo, selection.method = "vst", nfeatures = 3000)
  seo <- ScaleData(seo)
  
  # Extract the expression matrix of highly variable genes
  expr = GetAssayData(object = seo, slot = "scale.data")
  dim(expr)
  expr <- expr[rowSums(expr) > 0, ]
  
  n_cells <- nrow(t(expr))
  
  # Cosine similarity of gene expression
  build_expression_similarity <- function(spatial_expr,
                                          spatial_coords = NULL,
                                          method = "cosine") {
    expr_filtered <- spatial_expr
    
    # Cell standardization
    expr_norm <- t(apply(expr_filtered, 1, function(x) {
      (x - mean(x)) / sd(x)
    }))
    
    # Expression similarity calculation
    if (method == "cosine") {
      # Cosine similarity (In-memory Efficient implementation)
      sim_matrix <- cosine_similarity(expr_norm)
    } else if (method == "pearson") {
      # Pearson related
      sim_matrix <- cor(expr_norm, method = "pearson")
    } else if (method == "spearman") {
      # Spearman related
      sim_matrix <- cor(expr_norm, method = "spearman")
    } else {
      stop("Unsupported similarity method. Choose from: cosine, pearson, spearman")
    }
    
    return(sim_matrix)
  }
  
  # Auxiliary function: In-memory efficient cosine similarity calculation
  cosine_similarity <- function(mat) {
    # Normalized matrix
    norm_mat <- mat / sqrt(rowSums(mat^2))
    sim <- crossprod(norm_mat)
    diag(sim) <- 1
    return(sim)
  }
  
  expr_sim <- build_expression_similarity(expr)
  expr_sim[1:10, 1:10]
  
  # Spatial k-neighborhood graph
  k <- if (is.null(k))
    6
  else
    k
  
  # It does not contain its own cells
  nn_idx <- FNN::get.knn(spatial_coords, k = k)$nn.index
  
  spatial_adj <- Matrix(0, n_cells, n_cells, sparse = TRUE)
  
  for (i in seq_len(n_cells)) {
    spatial_adj[i, nn_idx[i, ]] <- 1
  }
  
  diag(spatial_adj) <- 0
  spatial_graph <- as.matrix(spatial_adj)
  
  # Fusion similarity
  spatial_weight <- if (is.null(spatial_weight))
    0.3
  else
    spatial_weight
  
  combined_sim <- spatial_weight * spatial_graph +
    (1 - spatial_weight) * expr_sim
  
  # Symmetry & binarization
  combined_sim <- pmax(combined_sim, t(combined_sim))
  threshold <- quantile(combined_sim[combined_sim > 0], 0.7)
  adj_bin <- (combined_sim > threshold) * 1
  
  # Build igraph
  g <- graph_from_adjacency_matrix(adj_bin, mode = "undirected", weighted = TRUE)
  
  # Multi-resolution Leiden clustering
  best_partition <- NULL
  best_mod <- -Inf
  
  suppressWarnings({
    for (res in c(0.5, 0.8, 1.0, 1.2, 1.5)) {
      partition <- leiden::leiden(
        g,
        partition_type = "RBConfiguration",
        weights = E(g)$weight,
        resolution_parameter = res,
        seed = 42
      )
      
      # mod is used to measure "how much higher the sum of the weights of the edges within a graph after it is divided into several communities (clusters) is than in a random situation."
      # The larger the Q, the closer the community is within it, the sparser it is between communities, and the higher the quality of division.
      # Q ≈ 0: Within the community, there is no difference from randomness → Division is meaningless.
      # Q > 0.3 is usually regarded as the threshold for a "significant community structure".
      # Q ≈ 1: An extremely ideal division.
      mod <- modularity(g, partition, weights = E(g)$weight)
      
      if (mod > best_mod) {
        best_mod <- mod
        best_partition <- partition
      }
      
    }
  })
  
  list(clusters = best_partition, graph = g)
}