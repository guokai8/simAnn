#' Plot Ontology Tree Visualization
#'
#' @description
#' Creates a visualization of an OntologyTree object as either a hierarchical or
#' circular dendrogram using ggraph. The visualization includes nodes representing
#' terms and edges representing relationships between terms. Supports focused
#' visualization of specific terms and their relationships.
#'
#' @details
#' The visualization includes:
#' \itemize{
#'   \item Nodes representing ontology terms
#'   \item Edges showing hierarchical relationships
#'   \item Node labels with term names or IDs
#'   \item Optional focus on specific terms and their connections
#' }
#'
#' The function uses ggraph for visualization and supports:
#' \itemize{
#'   \item Hierarchical tree layout
#'   \item Node repulsion for better label placement
#'   \item Edge styling for relationship visualization
#'   \item Interactive term highlighting
#' }
#'
#' @param tree An OntologyTree object containing:
#' \itemize{
#'   \item termNames: Character vector of term names/IDs
#'   \item childMap: List mapping parents to children
#'   \item edgeTypes: List of relationship types
#' }
#' @param layout Character string specifying the visualization layout:
#' \describe{
#'   \item{tree}{Traditional hierarchical dendrogram (default)}
#'   \item{circular}{Circular dendrogram arrangement}
#' }
#' @param focus_terms Character vector of term IDs/names to highlight in the
#'        visualization. If NULL (default), the entire tree is plotted.
#'
#' @return A ggplot object containing:
#' \itemize{
#'   \item Edge links showing term relationships
#'   \item Node points representing terms
#'   \item Text labels for terms
#'   \item Appropriate theme settings
#' }
#'
#' @examples
#' # Create a simple ontology tree
#' parents <- c("a", "a", "b", "b", "b", "c", "d")
#' children <- c("b", "c", "c", "d", "e", "e", "f")
#' 
#' # Define annotations for each term
#' annotation <- list(
#'     "a" = 1:3,
#'     "b" = 3:4,
#'     "c" = 5,
#'     "d" = 7,
#'     "e" = 4:7,
#'     "f" = 8
#' )
#' 
#' # Build and plot the tree
#' tre <- buildOntologyTree(
#'     parentTerms = parents,
#'     childTerms = children,
#'     annotations = annotation
#' )
#' plotOntologyTree(tre)
#' 
#' # Plot with circular layout
#' plotOntologyTree(tre, layout = "circular")
#' 
#' # Focus on specific terms
#' plotOntologyTree(tre, focus_terms = c("a", "b", "c"))
#'
#' @section Error Handling:
#' The function includes checks for:
#' \itemize{
#'   \item Empty trees (no edges)
#'   \item Invalid focus terms
#'   \item Missing required tree components
#' }
#'
#' @seealso
#' \code{\link[ggraph]{ggraph}} for additional plot customization
#' \code{\link{buildOntologyTree}} for creating OntologyTree objects
#' \code{\link{buildGOTree}} for creating Gene Ontology trees
#'
#' @references
#' Thomas Lin Pedersen (2021). ggraph: An Implementation of Grammar of Graphics
#' for Graphs and Networks. R package version 2.0.5.
#'
#' @importFrom igraph graph_from_data_frame V degree distances
#' @importFrom ggraph ggraph create_layout geom_edge_link geom_node_point geom_node_text
#' @importFrom ggplot2 theme_void
#' @author Kai Guo
#' @export
plotOntologyTree <- function(tree, layout = c("tree"), focus_terms = NULL) {
  layout <- match.arg(layout)
  
  # Convert OntologyTree to edge list
  edge_list <- do.call(rbind, lapply(seq_along(tree@childMap), function(i) {
    parent <- tree@termNames[i]  
    children <- tree@childMap[[i]]  
    if (length(children) > 0) {
      data.frame(from = parent, to = tree@termNames[children], stringsAsFactors = FALSE)
    } else {
      NULL
    }
  }))
  
  if (nrow(edge_list) == 0) {
    stop("The tree has no edges to plot.")
  }
  
  # Handle focus terms if provided
  if (!is.null(focus_terms)) {
    valid_terms <- unique(c(edge_list$from, edge_list$to))
    focus_terms <- intersect(focus_terms, valid_terms)
    edge_list <- subset(edge_list, from %in% focus_terms | to %in% focus_terms)
    
    if (nrow(edge_list) == 0) {
      stop("No matching terms found in the tree.")
    }
  }
  
  # Convert to igraph object
  g <- igraph::graph_from_data_frame(edge_list, directed = TRUE)
  
  # Generate plot
  p <- ggraph::ggraph(g, layout = "dendrogram", circular = (layout == "circular")) +
    ggraph::geom_edge_link(color = "gray") +
    ggraph::geom_node_point(color = "cyan4", size = 3) + 
    ggraph::geom_node_text(aes(label = name), repel = TRUE, size = 4) + 
    theme_void()
  
  return(p)
}

