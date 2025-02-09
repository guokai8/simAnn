#' Threshold-based Clustering using Graph Community Detection
#'
#' @description
#' This function performs clustering on a similarity matrix based on a specified threshold.
#' The similarity matrix is first binarized (values \eqn{\geq}{>=} the threshold are set to 1),
#' an undirected graph is constructed from the binary matrix, and then a community detection algorithm is applied.
#'
#' @param sim_mat A numeric matrix representing pairwise similarities. Rows and columns should correspond to the same set of objects.
#' @param threshold A numeric value between 0 and 1. Pairs with similarity greater than or equal to this threshold are connected.
#' @param method A character string specifying the community detection method to use. Options include:
#'   \itemize{
#'     \item \code{"components"}: Uses connected components (the simplest approach).
#'     \item \code{"louvain"}: Uses the Louvain algorithm for community detection.
#'     \item \code{"walktrap"}: Uses the Walktrap algorithm based on random walks.
#'   }
#'   Default is \code{"components"}.
#' @param remove_self Logical. If \code{TRUE} (the default), self-loops are removed by setting the diagonal of \code{sim_mat} to zero.
#'
#' @return A vector of cluster memberships for each object.
#'
#' @examples
#' \dontrun{
#'   sim_mat <- matrix(runif(100, 0, 1), nrow = 10)
#'   diag(sim_mat) <- 1
#'   clusters <- threshold_cluster(sim_mat, threshold = 0.6, method = "components")
#' }
#'
#' @importFrom igraph graph_from_adjacency_matrix components cluster_louvain cluster_walktrap membership
#' @export
threshold_cluster <- function(sim_mat, threshold, 
                              method = c("components", "louvain", "walktrap"),
                              remove_self = TRUE) {
  method <- match.arg(method)
  
  if (!is.matrix(sim_mat)) {
    stop("sim_mat must be a matrix.")
  }
  if (!is.numeric(threshold) || threshold < 0 || threshold > 1) {
    stop("threshold must be a numeric value between 0 and 1.")
  }
  
  # Binarize the similarity matrix: values >= threshold are set to 1, otherwise 0.
  binary_mat <- ifelse(sim_mat >= threshold, 1, 0)
  if (remove_self) diag(binary_mat) <- 0
  
  # Construct an undirected graph from the binary matrix.
  g <- igraph::graph_from_adjacency_matrix(binary_mat, mode = "undirected", diag = FALSE)
  
  # Apply the chosen community detection method.
  if (method == "components") {
    comps <- igraph::components(g)
    clustering <- comps$membership
  } else if (method == "louvain") {
    cl <- igraph::cluster_louvain(g)
    clustering <- igraph::membership(cl)
  } else if (method == "walktrap") {
    cl <- igraph::cluster_walktrap(g)
    clustering <- igraph::membership(cl)
  } else {
    stop("Unsupported method provided.")
  }
  
  return(clustering)
}

#############################################
# quality_indices.R
#############################################

#' Calinski–Harabasz Index from MDS Coordinates
#'
#' @description
#' Computes the Calinski–Harabasz (CH) index for a clustering based on coordinates obtained
#' from classical multidimensional scaling (cmdscale) applied to the distance matrix defined as 1 - sim_mat.
#' A higher CH index indicates better cluster separation.
#'
#' @param sim_mat A numeric similarity matrix.
#' @param cl A vector of cluster labels.
#'
#' @return A numeric scalar representing the CH index.
#'
#' @details
#' The function first converts the similarity matrix to a distance matrix (using 1 - sim_mat),
#' performs cmdscale to obtain coordinates in (by default) 2 dimensions, and then computes the between-cluster
#' and within-cluster sum of squares to form the CH index.
#'
#' @examples
#' \dontrun{
#'   sim_mat <- matrix(runif(100, 0, 1), nrow = 10)
#'   diag(sim_mat) <- 1
#'   cl <- threshold_cluster(sim_mat, threshold = 0.7)
#'   ch <- ch_index(sim_mat, cl)
#' }
#'
#' @export
ch_index <- function(sim_mat, cl) {
  if (!is.matrix(sim_mat) || !is.numeric(sim_mat)) {
    stop("sim_mat must be a numeric matrix.")
  }
  if (nrow(sim_mat) != ncol(sim_mat)) {
    stop("sim_mat must be square.")
  }
  if (length(cl) != nrow(sim_mat)) {
    stop("Length of cl must equal the number of rows in sim_mat.")
  }
  
  # Convert similarity to distance (assumes sim_mat is normalized between 0 and 1)
  dist_mat <- as.dist(1 - sim_mat)
  
  # Use classical MDS to get coordinates (we choose 2 dimensions by default)
  coords <- cmdscale(dist_mat, k = min(nrow(sim_mat) - 1, 2))
  n <- nrow(coords)
  overall_mean <- colMeans(coords)
  clusters <- unique(cl)
  k <- length(clusters)
  
  if (k < 2) {
    return(NA_real_)
  }
  
  wss <- 0  # within-cluster sum of squares
  bss <- 0  # between-cluster sum of squares
  
  for (cluster in clusters) {
    cluster_coords <- coords[cl == cluster, , drop = FALSE]
    cluster_mean <- colMeans(cluster_coords)
    wss <- wss + sum(rowSums((cluster_coords - matrix(cluster_mean, nrow(cluster_coords),
                                                      ncol(cluster_coords), byrow = TRUE))^2))
    bss <- bss + nrow(cluster_coords) * sum((cluster_mean - overall_mean)^2)
  }
  
  ch <- (bss / (k - 1)) / (wss / (n - k))
  return(ch)
}

#' Davies–Bouldin Quality Index (Reciprocal)
#'
#' @description
#' Computes a quality measure based on the Davies–Bouldin (DB) index.
#' The DB index is defined as the average ratio of the sum of intra-cluster dispersions to the inter-cluster centroid distance.
#' Because lower DB values indicate better clustering, this function returns the reciprocal (1 / DB) so that higher values indicate better quality.
#'
#' @param sim_mat A numeric similarity matrix.
#' @param cl A vector of cluster labels.
#'
#' @return A numeric scalar representing the reciprocal of the DB index.
#'
#' @details
#' The function first converts the similarity matrix to a distance matrix (using 1 - sim_mat),
#' obtains coordinates via cmdscale (defaulting to 2 dimensions), computes the average dispersion for each cluster,
#' and then computes the DB index. The reciprocal is returned as the quality measure.
#'
#' @examples
#' \dontrun{
#'   sim_mat <- matrix(runif(100, 0, 1), nrow = 10)
#'   diag(sim_mat) <- 1
#'   cl <- threshold_cluster(sim_mat, threshold = 0.7)
#'   quality_db <- db_index(sim_mat, cl)
#' }
#'
#' @export
db_index <- function(sim_mat, cl) {
  if (!is.matrix(sim_mat) || !is.numeric(sim_mat)) {
    stop("sim_mat must be a numeric matrix.")
  }
  if (nrow(sim_mat) != ncol(sim_mat)) {
    stop("sim_mat must be square.")
  }
  if (length(cl) != nrow(sim_mat)) {
    stop("Length of cl must equal the number of rows in sim_mat.")
  }
  
  dist_mat <- as.dist(1 - sim_mat)
  coords <- cmdscale(dist_mat, k = min(nrow(sim_mat) - 1, 2))
  clusters <- unique(cl)
  k <- length(clusters)
  if (k < 2) return(NA_real_)
  
  # Compute centroids for each cluster
  centroids <- t(sapply(clusters, function(cluster) {
    colMeans(coords[cl == cluster, , drop = FALSE])
  }))
  
  # Compute average dispersion (S_i) for each cluster
  S <- sapply(clusters, function(cluster) {
    cluster_coords <- coords[cl == cluster, , drop = FALSE]
    centroid <- colMeans(cluster_coords)
    mean(sqrt(rowSums((cluster_coords - matrix(centroid, nrow(cluster_coords),
                                               ncol(cluster_coords), byrow = TRUE))^2)))
  })
  
  # Compute distance between centroids
  M <- as.matrix(dist(centroids))
  
  R <- numeric(k)
  for (i in 1:k) {
    R[i] <- max(sapply(1:k, function(j) {
      if (i == j) return(-Inf)
      (S[i] + S[j]) / M[i, j]
    }))
  }
  db <- mean(R)
  quality <- 1 / db  # Reciprocal: higher values indicate better clustering.
  return(quality)
}

#############################################
# select_best_threshold.R
#############################################

#' Select the Best Threshold for Clustering
#'
#' @description
#' Evaluates a range of threshold values to determine the optimal threshold for clustering a similarity matrix.
#' For each threshold, the function applies the clustering algorithm (via \code{threshold_cluster}) and computes
#' one of several quality metrics. The available quality metrics include:
#' \itemize{
#'   \item \code{"modularity"}: Uses the modularity value computed by \code{igraph::modularity}.
#'   \item \code{"silhouette"}: Uses the average silhouette width computed on a distance matrix defined as \code{1 - sim_mat}.
#'   \item \code{"ch_index"}: Uses the Calinski–Harabasz index computed from coordinates obtained by cmdscale.
#'   \item \code{"db_index"}: Uses the reciprocal of the Davies–Bouldin index (so that higher is better).
#' }
#'
#' @param sim_mat A numeric similarity matrix.
#' @param thresholds A numeric vector of threshold values to evaluate. Default is \code{seq(0.5, 0.9, by = 0.01)}.
#' @param cluster_method A character string specifying the clustering method to use with \code{threshold_cluster}. Options:
#'   \code{"components"}, \code{"louvain"}, or \code{"walktrap"}. Default is \code{"components"}.
#' @param quality_metric A character string specifying the quality metric to use for evaluating clustering results.
#'   Options are: \code{"modularity"}, \code{"silhouette"}, \code{"ch_index"}, or \code{"db_index"}. Default is \code{"modularity"}.
#' @param remove_self Logical. If \code{TRUE} (the default), self-loops are removed by setting the diagonal of \code{sim_mat} to 0.
#' @param verbose Logical. If \code{TRUE}, prints progress messages. Default is \code{TRUE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{best_threshold}{The threshold value that maximizes the selected quality metric.}
#'   \item{thresholds}{The evaluated threshold values.}
#'   \item{quality_values}{The computed quality metric values for each threshold.}
#' }
#'
#' @examples
#' \dontrun{
#'   sim_mat <- matrix(runif(100, 0, 1), nrow = 10)
#'   diag(sim_mat) <- 1
#'   result <- pickThreshold(sim_mat, thresholds = seq(0.5, 0.9, by = 0.01),
#'                                   cluster_method = "components", quality_metric = "modularity")
#'   best_th <- result$best_threshold
#' }
#'
#' @importFrom igraph graph_from_adjacency_matrix modularity
#' @importFrom cluster silhouette
#' @import ggplot2
#' @export
pickThreshold <- function(sim_mat,
                                  thresholds = seq(0.5, 0.9, by = 0.01),
                                  cluster_method = c("components", "louvain", "walktrap"),
                                  quality_metric = c("modularity", "silhouette", "ch_index", "db_index"),
                                  remove_self = TRUE,
                                  verbose = TRUE) {
  cluster_method <- match.arg(cluster_method)
  quality_metric <- match.arg(quality_metric)
  
  quality_values <- numeric(length(thresholds))
  
  for (i in seq_along(thresholds)) {
    th <- thresholds[i]
    
    # Perform clustering at the current threshold.
    clustering <- threshold_cluster(sim_mat, threshold = th,
                                    method = cluster_method,
                                    remove_self = remove_self)
    
    if (quality_metric == "modularity") {
      # Reconstruct the binary graph for modularity computation.
      binary_mat <- ifelse(sim_mat >= th, 1, 0)
      if (remove_self) diag(binary_mat) <- 0
      g <- igraph::graph_from_adjacency_matrix(binary_mat, mode = "undirected", diag = FALSE)
      quality_values[i] <- igraph::modularity(g, clustering)
      
    } else if (quality_metric == "silhouette") {
      # Compute the silhouette width using the distance matrix (1 - sim_mat).
      dist_mat <- as.dist(1 - sim_mat)
      if (length(unique(clustering)) < 2) {
        quality_values[i] <- NA_real_
      } else {
        sil <- cluster::silhouette(clustering, dist_mat)
        quality_values[i] <- mean(sil[, 3])
      }
      
    } else if (quality_metric == "ch_index") {
      quality_values[i] <- ch_index(sim_mat, clustering)
      
    } else if (quality_metric == "db_index") {
      quality_values[i] <- db_index(sim_mat, clustering)
    }
  }
  
  best_index <- which.max(quality_values)
  best_threshold <- thresholds[best_index]
  
  # Directly use ggplot2 (imported via DESCRIPTION) to plot the quality metric vs. threshold.
  df <- data.frame(threshold = thresholds, quality = quality_values)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = threshold, y = quality)) +
    ggplot2::geom_line(color = "blue") +
    ggplot2::geom_point(color = "blue") +
    ggplot2::labs(title = sprintf("Quality Metric (%s) vs. Threshold", quality_metric),
                  x = "Threshold", y = quality_metric) +
    ggplot2::theme_minimal()
  print(p)
  
  message(sprintf("Optimal threshold: %.3f with %s value: %.4f",
                  best_threshold, quality_metric, quality_values[best_index]))
  
  return(list(best_threshold = best_threshold,
              thresholds = thresholds,
              quality_values = quality_values))
}
