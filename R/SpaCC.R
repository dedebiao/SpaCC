#' @include FindSpaCC.R
#'
NULL

##' Run SpaCC clustering algorithm
##'
##' @description Implements the SpaCC algorithm in R using reticulate to run the Python version. Requires the python "leidenalg", "igraph" and "leiden" modules to be installed. Returns a vector of partition indices.
##' Windows users can still this with devtools::install_github("rstudio/reticulate", ref = "86ebb56"); reticulate::use_condaenv("r-reticulate"); reticulate::conda_install("r-reticulate", "leidenalg", channel = "vtraag")
##' @param optional weights are used to apply "influence" or "confidence" to different nodes on a spatial graph. NULL indicates that all nodes have equal weights.
##' @param k The number of neighbors of each node when constructing a KNN graph. If set to NULL, the function will automatically try multiple k values internally and select the best one. The larger k is, the denser the graph becomes. The smaller the k, the sparser the graph.
##'
##' Output the cluster label vector (integer or factor) returned by FindSpaCC().
##' Output an igraph object that includes: node = cell /spot, edge = spatial adjacency relationship, edge weight = spatial distance or expression similarity.
##' The output file function writes the edge list file (from, to, weight) in the running directory, which is convenient for visualization in Cytoscape/Gephi later.
SpaCC <- function(seo,
                  spatial_coords,
                  spatial_weight = NULL,
                  k = NULL) {
  result <- FindSpaCC(
    seo = seo,
    spatial_coords = spatial_coords,
    spatial_weight = spatial_weight,
    k = k
  )
  
SpaCC_graph <- function(result) {
    if (!"graph" %in% names(result))
      stop("FindSpaCC returns the missing 'graph' field in the object.")
    
    g <- result$graph
    
    # Force vertex names to be characters to avoid igraph warnings
    if (!is.character(V(g)$name)) {
      V(g)$name <- as.character(V(g)$name)
    }
    
    # Optional: Remove self-loop (Many analyses do not require self-loop)
    # g <- simplify(g, remove.multiple = FALSE, remove.loops = TRUE)
    
    g
  }
  g <- SpaCC_graph(result)
  edges <- as.data.frame(get.edgelist(g))
  colnames(edges) <- c("from", "to")
  edges$weight <- E(g)$weight
  write.csv(edges, 'Graph_connection_matrix.csv')
  seo[["SpaCC"]] <- result$clusters
  
  message("Number of clusters: ", length(unique(result$clusters)))
  print(table(result$clusters))
  
 seo
  suppressWarnings({
    return(seo)
  })
}