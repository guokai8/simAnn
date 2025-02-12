#' Cluster Semantically Similar Terms
#'
#' @description
#' Groups semantically similar terms based on similarity scores and ontological
#' relationships. The function supports multiple similarity calculation methods
#' and returns clustering results in multiple formats.
#'
#' @details
#' The function performs hierarchical clustering of terms based on their semantic
#' similarity. It considers both the similarity scores and the ontological structure
#' when forming clusters. Available similarity methods include term-based, information
#' content-based, and graph-based approaches.
#'
#' @param tree An OntologyTree object containing the ontology structure and annotations.
#' @param terms A character vector of term IDs/names or numeric indices to be clustered.
#' @param method Character string specifying the similarity calculation method:
#' \describe{
#'   \item{term}{Term overlap based methods}
#'   \item{ic}{Universal information content (\doi{10.1186/1471-2105-9-S5-S4})}
#'   \item{resnik}{Resnik's information content (\doi{10.1613/jair.514})}
#'   \item{lin}{Lin similarity (\doi{10.5555/645527.657297})}
#'   \item{faith}{FaITH similarity (\doi{10.1007/978-3-642-17746-0_39})}
#'   \item{rel}{Relevance similarity (\doi{10.1186/1471-2105-7-302})}
#'   \item{simic}{SimIC similarity (\doi{10.48550/arXiv.1001.0958})}
#'   \item{gogo}{GOGO similarity (\doi{10.1038/s41598-018-33219-y})}
#'   \item{wang}{Wang's method (\doi{10.1093/bioinformatics/btm087})}
#'   \item{anc}{Ancestor based similarity}
#' }
#'
#' @param ic_method For "faith", "rel", and "simic" methods:
#' \describe{
#'   \item{universal}{Calculate universal-based Information Content}
#'   \item{annotation}{Calculate annotation-based Information Content}
#'   \item{wang}{Wang's Information Content}
#'   \item{offspring}{Calculate offspring-based Information Content}
#' }
#'
#' @param submethod For "term" method only:
#' \describe{
#'   \item{kappa}{Cohen's kappa coefficient between annotation sets}
#'   \item{overlap}{Overlap coefficient: |A∩B| / min(|A|,|B|)}
#'   \item{jaccard}{Jaccard similarity: |A∩B| / |A∪B|}
#'   \item{dice}{Dice coefficient: 2|A∩B| / (|A|+|B|)}
#' }
#'
#' @param normType For Resnik similarity only:
#' \describe{
#'   \item{unif}{Normalize by log of max annotation count}
#'   \item{max}{Normalize by maximum IC value}
#'   \item{sum}{Normalize by maximum pairwise IC}
#'   \item{none}{No normalization}
#' }
#'
#' @param weights Named numeric vector of weights for GOGO/Wang methods.
#'        Names must match relation types (e.g., "is_a", "part_of"), values < 1.
#' @param threshold Numeric value for similarity threshold in clustering.
#'        Terms with similarity above this threshold will be grouped together.
#' @param useIgraph Logical indicating whether to use igraph package for calculations.
#'        Recommended for large graphs. Default is FALSE.
#' @param useCache Logical indicating whether to use cached IC values when applicable.
#'        Default is TRUE.
#' @param verbose Logical indicating whether to print progress messages.
#'        Default is FALSE.
#'
#' @return A list containing multiple representations of clustering results:
#' \describe{
#'   \item{clusters}{List of term clusters, where each element contains terms in a cluster}
#'   \item{membership}{Named vector mapping terms to their cluster assignments}
#'   \item{similarity}{Matrix of pairwise similarity scores between terms}
#'   \item{graph}{igraph object representing the clustering structure}
#'   \item{dendrogram}{Dendrogram object from hierarchical clustering}
#'   \item{stats}{Summary statistics for the clustering results}
#' }
#'
#' @examples
#' \dontrun{
#' # Load required packages
#' library(org.Hs.eg.db)
#'
#' # Build GO tree for human biological process terms
#' tree <- buildGOTree(species = "human", namespace = "BP")
#'
#' # Define edge weights for Wang's method
#' weights <- c("is_a" = 0.8, "part_of" = 0.6)
#'
#' # Generate clusters using Wang's method
#' clusters <- clusterST(tree,
#'                      terms = 1:30,
#'                      method = "wang",
#'                      weights = weights,
#'                      threshold = 0.2)
#'
#' # Example with term overlap method
#' clusters2 <- clusterST(tree,
#'                       terms = c("GO:0006915", "GO:0012501", "GO:0043065"),
#'                       method = "term",
#'                       submethod = "jaccard",
#'                       threshold = 0.5)
#'
#' }
#'
#' @references
#' \itemize{
#'   \item Wang JZ, et al (2007). A new method to measure the semantic similarity of GO terms.
#'   \item Yang H, et al (2012). GOGO: An improved algorithm to measure the semantic similarity.
#'   \item Resnik P (1995). Using information content to evaluate semantic similarity in a taxonomy.
#' }
#'
#' @seealso
#' \code{\link{simterm}} for calculating semantic similarity
#' \code{\link{plotCluster}} for visualizing clustering results
#' 
#' @return a list include the TermClusters 
#' @author Kai Guo
#' @export
clusterST <- function(
    tree,
    terms,
    method = c("term", "ic", "resnik", "lin","faith", "rel", "simic", "gogo", "wang", "anc"),
    threshold = 0.8,
    ic_method = NULL,
    submethod = c("kappa", "overlap", "jaccard", "dice"),
    normType = "max",
    weights = NULL,
    annotUniverse= NULL, 
    useCache = TRUE,
    useIgraph = TRUE,
    verbose = TRUE
) {
  # Input validation
  method <- match.arg(method)
  if(is.null(ic_method)){
    ic_method = defaultIC(tree)
  }
  if(!inherits(tree, "OntologyTree")) {
    stop("'tree' must be an OntologyTree object")
  }
  
  # Convert terms to indices if necessary
  termIds <- if(is.character(terms)) {
    term2id(tree, terms, strict = FALSE)
  } else {
    terms
  }
  icScores=calculateIC(tree, method = ic_method, verbose = verbose)
  if(ic_method == "annotation"){
    validTerms <- validateAnnotatedTerms(tree, termIds)
    termIdx <- termIds[validTerms]
  }
  if(length(termIds) == 0) {
    stop("No valid terms provided")
  }
  
  # Get term names from metadata
  if(is.null(tree@metadata)){
    term_names <- tree@termNames
    names(term_names) <- tree@termNames
  }else{
    term_names <- tree@metadata$TermAnnot
    names(term_names) <- tree@metadata$TermID
  }
  # Calculate similarity matrix
  sim_matrix <- simterm(
    tree = tree,
    terms = termIds,
    method = method,
    ic_method= ic_method,
    submethod = if(method == "term") match.arg(submethod),
    normType = normType,
    annotUniverse = annotUniverse,
    weights = weights,
    useCache = useCache,
    useIgraph = useIgraph,
    verbose = verbose
  )
  if(length(termIdx)!=length(termIds)){
    cat("Some Terms don't have ", ic_method, " IC values|\n")
    print("######################")
    print(tree@termNames[setdiff(termIds,termIdx)])
    print("######################")
  }
  result <- clusterSimilarTermsCPP(sim_matrix, termIds, 
                                    parentMap = tree@parentMap, 
                                    termNames = term_names, 
                                    icScores = icScores, 
                                    annotationsList = tree@annotations$list,
                                    annotationsNames =tree@annotations$names,
                                    threshold = threshold)
  
  attr(result, "method") <- method
  attr(result, "threshold") <- threshold
  attr(result, "timestamp") <- Sys.time()
  
  class(result) <- c("TermClusters", "list")
  return(result)
}
#' Summary method for TermClusters objects
#' @export
summary.TermClusters <- function(object, ...) {
  n_clusters <- length(object)
  all_children <- do.call(rbind, lapply(object, function(x) x$children))
  
  structure(list(
    n_clusters = n_clusters,
    n_total_terms = nrow(all_children) + n_clusters,
    avg_cluster_size = nrow(all_children) / n_clusters,
    similarity_stats = list(
      min = min(all_children$similarity),
      max = max(all_children$similarity),
      mean = mean(all_children$similarity),
      median = median(all_children$similarity)
    ),
    ic_stats = list(
      min = min(all_children$ic),
      max = max(all_children$ic),
      mean = mean(all_children$ic),
      median = median(all_children$ic)
    ),
    method = attr(object, "method"),
    threshold = attr(object, "threshold")
  ), class = "summary.TermClusters")
}

#' Print method for TermClusters summary
#' @author Kai Guo
#' @export
print.summary.TermClusters <- function(x, ...) {
  cat("Term Clustering Summary\n\n")
  cat("Method:", x$method, "\n")
  cat("Similarity threshold:", x$threshold, "\n")
  cat("Number of clusters:", x$n_clusters, "\n")
  cat("Total terms:", x$n_total_terms, "\n")
  cat("Average cluster size:", round(x$avg_cluster_size, 2), "\n\n")
  
  cat("Similarity statistics:\n")
  cat("  Min:", round(x$similarity_stats$min, 3), "\n")
  cat("  Max:", round(x$similarity_stats$max, 3), "\n")
  cat("  Mean:", round(x$similarity_stats$mean, 3), "\n")
  cat("  Median:", round(x$similarity_stats$median, 3), "\n\n")
  
  cat("Information Content statistics:\n")
  cat("  Min:", round(x$ic_stats$min, 3), "\n")
  cat("  Max:", round(x$ic_stats$max, 3), "\n")
  cat("  Mean:", round(x$ic_stats$mean, 3), "\n")
  cat("  Median:", round(x$ic_stats$median, 3), "\n")
}

#' Plot Similarity Heatmap Based on Clustering
#'
#' @description
#' Creates a heatmap visualization of term similarities, organized by clustering order.
#' The heatmap can display term names or IDs, includes cluster annotations, and
#' supports multiple color schemes. When available, uses the pheatmap package for
#' enhanced visualization features.
#'
#' @details
#' The heatmap visualization includes:
#' \itemize{
#'   \item Similarity scores represented by color intensity
#'   \item Cluster annotations on rows and columns
#'   \item Optional term names instead of IDs
#'   \item Customizable color schemes for both similarity values and cluster annotations
#' }
#'
#' @param cluster_result The result object from clusterST containing:
#' \itemize{
#'   \item similarity_matrix: Matrix of similarity scores
#'   \item clusters: List of cluster information
#'   \item cluster_matrix: Data frame with term mappings
#' }
#' @param show_names Logical indicating whether to show term names instead of IDs.
#'        Default: FALSE.
#' @param annotation_style Character string specifying how to display cluster annotations:
#' \describe{
#'   \item{combined}{Shows "Cluster N - Term Name" format}
#'   \item{separate}{Shows cluster number and term name in separate annotation columns}
#' }
#' @param color_scheme Character string specifying the color scheme for cluster annotation:
#' \describe{
#'   \item{default}{35 colors optimized for distinction}
#'   \item{brewer}{RColorBrewer palettes (Set1, Set2, Set3, Dark2, Accent)}
#'   \item{custom1}{20 colors focused on blues and earth tones}
#'   \item{custom2}{20 colors with muted tones}
#' }
#' @param custom_colors Character vector of custom colors for cluster annotation.
#'        If provided, overrides color_scheme.
#' @param cex Numeric value for scaling text size.
#'        Default: 0.8.
#' @param margins Numeric vector of length 2 specifying plot margins.
#'        Default: c(8, 8).
#' @param colors Color function or vector for similarity values.
#'        Default: blue-white-red gradient with 100 colors.
#' @param ... Additional parameters passed to pheatmap::pheatmap().
#'
#' @return Invisibly returns the reordered similarity matrix used for plotting.
#'
#' @examples
#' \dontrun{
#' # Basic usage with default settings
#' plotClusterHeatmap(cluster_result)
#'
#' # Show term names with combined annotation
#' plotClusterHeatmap(cluster_result,
#'                    show_names = TRUE,
#'                    annotation_style = "combined")
#'
#' # Use RColorBrewer color scheme
#' plotClusterHeatmap(cluster_result,
#'                    color_scheme = "brewer",
#'                    show_names = TRUE)
#'
#' # Custom colors for cluster annotation
#' custom_cols <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")
#' plotClusterHeatmap(cluster_result,
#'                    custom_colors = custom_cols,
#'                    annotation_style = "separate")
#'
#' # Modify heatmap appearance
#' plotClusterHeatmap(cluster_result,
#'                    cex = 1.0,
#'                    margins = c(10, 10),
#'                    colors = colorRampPalette(c("navy", "white", "darkred"))(50))
#' }
#'
#' @seealso
#' \code{\link[pheatmap]{pheatmap}} for additional plotting parameters
#' \code{\link{clusterST}} for generating cluster results
#'
#' @references
#' Raivo Kolde (2019). pheatmap: Pretty Heatmaps.
#' R package version 1.0.12. https://CRAN.R-project.org/package=pheatmap
#'
#' @importFrom stats hclust
#' @importFrom pheatmap pheatmap
#' @author Kai Guo
#' @export
plotClusterHeatmap <- function(cluster_result, 
                               show_names = FALSE,
                               annotation_style = c("combined", "separate"),
                               color_scheme = c("default", "brewer", "custom1", "custom2"),
                               custom_colors = NULL,
                               cex = 0.8,
                               margins = c(8, 8),
                               colors = colorRampPalette(c("#313695", "#FFFFFF", "#A50026"))(100),...) {
  
  annotation_style <- match.arg(annotation_style)
  color_scheme <- match.arg(color_scheme)
  
  # 
  color_schemes <- list(
    default = c("#A6761D","#D95F02","deepskyblue","#1B9E77","#E7298A","#7570B3","#E6AB02",
                "#FC9D9A","#83AF9B","#FE4365","#F9CDAD","#C8C8A9","#B6C29A","#8A977B",
                "#F4D000","#E58308","#DC5712","#A6CEE3","#1F78B4","#B2DF8A","#33A02C",
                "slateblue1","darkgreen","darkred","plum1","#FB9A99","#E31A1C",
                "darkmagenta","hotpink2","slategray4","magenta2","yellow4",
                "#66A61E","#E7298A","#7570B3","#E6AB02",'#e6194b','#3cb44b','#ffe119',
                '#4363d8','#f58231','#911eb4','#46f0f0','#f032e6','#bcf60c',
                '#fabebe','#008080','#e6beff','#9a6324','#fffac8','#800000',
                '#aaffc3','#808000','#ffd8b1','#000075','#808080'),
    brewer = c(RColorBrewer::brewer.pal(8, "Set2"),
               RColorBrewer::brewer.pal(8, "Set3"),
               RColorBrewer::brewer.pal(8, "Set1"),
               RColorBrewer::brewer.pal(8, "Dark2"),
               RColorBrewer::brewer.pal(8, "Accent")),
    custom1 = c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
                "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
                "#4DC7C7", "#9D7660", "#D7B5A6", "#748CB6", "#D4A6C8",
                "#D5A548", "#8C4646", "#B3B3B3", "#579575", "#839557"),
    custom2 = c("#1A5878", "#C44237", "#AD8941", "#E99093", "#50594B",
                "#2F4858", "#7BA696", "#A54657", "#946B54", "#E4B680",
                "#8E402A", "#749C75", "#6B6B4A", "#9B6B6B", "#E79C46",
                "#D36E3B", "#8F8FA3", "#D1A33D", "#436978", "#7E6B44")
  )
  
  # 
  sim_matrix <- cluster_result$similarity_matrix
  clusters <- cluster_result$clusters
  
  # 
  terms_and_groups <- list()
  group_names <- list()  # 
  group_counter <- 1
  
  # 
  for(cluster in clusters) {
    if(cluster$relation_type != "singleton") {
      curr_terms <- c(cluster$cluster_term, cluster$children)
      group_label <- if(annotation_style == "combined") {
        paste("Cluster", group_counter, "-", cluster$term_name)
      } else {
        paste("Cluster", group_counter)
      }
      
      # 
      group_names[[group_label]] <- cluster$term_name
      
      for(term in curr_terms) {
        if(!term %in% names(terms_and_groups)) {
          terms_and_groups[[term]] <- group_label
        }
      }
      group_counter <- group_counter + 1
    }
  }
  
  # 
  for(cluster in clusters) {
    if(cluster$relation_type == "singleton") {
      term <- cluster$cluster_term
      if(!term %in% names(terms_and_groups)) {
        terms_and_groups[[term]] <- "Singleton"
        group_names[["Singleton"]] <- "Individual Terms"
      }
    }
  }
  
  # 
  term_order <- names(terms_and_groups)
  group_info <- unname(unlist(terms_and_groups))
  
  # 
  ordered_matrix <- sim_matrix[term_order, term_order]
  
  # 
  if(show_names) {
    id_to_name <- setNames(
      c(cluster_result$cluster_matrix$parent_name, 
        cluster_result$cluster_matrix$child_name),
      c(cluster_result$cluster_matrix$parent, 
        cluster_result$cluster_matrix$child)
    )
    id_to_name <- id_to_name[!duplicated(names(id_to_name))]
    rownames(ordered_matrix) <- colnames(ordered_matrix) <- id_to_name[term_order]
  }
  
  # 
  unique_groups <- unique(group_info[group_info != "Singleton"])
  n_groups <- length(unique_groups) 
  
  if (!is.null(custom_colors)) {
    selected_colors <- custom_colors
  } else {
    selected_colors <- color_schemes[[color_scheme]]
  }
  if (n_groups > length(selected_colors)) {
    warning("Not enough colors provided. Recycling colors.")
    selected_colors <- rep(selected_colors, length.out = n_groups)
  }
  
  if(n_groups > 0) {
    cluster_colors <- selected_colors[1:n_groups]
    names(cluster_colors) <- unique_groups
    cluster_colors["Singleton"] <- "gray80"
  } else {
    cluster_colors <- c(Singleton = "gray80")
  }
  
  # 
  if(requireNamespace("pheatmap", quietly = TRUE)) {
    if(annotation_style == "combined") {
      annotation_row <- data.frame(
        Cluster = factor(group_info),
        row.names = rownames(ordered_matrix)
      )
      
      annotation_colors <- list(
        Cluster = cluster_colors
      )
      
    } else {  # separate style
      # 
      annotation_row <- data.frame(
        Group = factor(group_info),
        Name = unname(sapply(group_info, function(x) group_names[[x]])),
        row.names = rownames(ordered_matrix)
      )
      
      annotation_colors <- list(
        Group = cluster_colors
      )
    }
    
    pheatmap::pheatmap(ordered_matrix,
                       color = colors,
                       annotation_row = annotation_row,
                       annotation_col = annotation_row,
                       annotation_colors = annotation_colors,
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       show_rownames = TRUE,
                       show_colnames = TRUE,
                       fontsize = cex * 10,
                       fontsize_row = cex * 8,
                       fontsize_col = cex * 8,...)
  } else {
    # 
    image(1:nrow(ordered_matrix), 
          1:ncol(ordered_matrix), 
          ordered_matrix,
          col = colors,
          axes = FALSE,
          xlab = "",
          ylab = "")
    
    axis(1, 1:nrow(ordered_matrix), 
         labels = rownames(ordered_matrix), 
         las = 2, 
         cex.axis = cex)
    axis(2, 1:ncol(ordered_matrix), 
         labels = colnames(ordered_matrix), 
         las = 2, 
         cex.axis = cex)
    
    box()
  }
  
  # 
  invisible(ordered_matrix)
}

#' Print method for TermClusters objects
#' @author Kai Guo
#' @export
print.TermClusters <- function(x, ...) {
  cat("Term Clustering Results\n")
  cat("Method:", attr(x, "method"), "\n")
  cat("Threshold:", attr(x, "threshold"), "\n")
  cat("\nClustering Statistics:\n")
  cat("  Total clusters:", x$stats$total_clusters, "\n")
  cat("  Parent-child clusters:", x$stats$parent_child_clusters, "\n")
  cat("  Similarity-based clusters:", x$stats$similarity_clusters, "\n")
  cat("  Singleton terms:", x$stats$singleton_terms, "\n")
  cat("  Mean similarity:", round(x$stats$mean_similarity, 3), "\n")
  
  cat("\nAvailable components:\n")
  cat("  $clusters - Detailed clustering results\n")
  cat("  $summary - Summary data frame\n")
  cat("  $cluster_matrix - Detailed parent-child relationships\n")
  cat("  $similarity_matrix - Original similarity matrix\n")
  cat("  $stats - Statistical summary\n")
}
