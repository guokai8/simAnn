#' @title Plot Ontology Tree
#' @description Visualizes an OntologyTree as either a hierarchical tree or a circular dendrogram.
#'
#' @param tree An OntologyTree object representing the hierarchical ontology structure.
#' @param layout Character, either "tree" (default) for a hierarchical dendrogram or "circular" for a circular dendrogram.
#' @param focus_terms Optional character vector specifying terms to highlight in the tree. If NULL, the entire tree is plotted.
#'
#' @return A ggplot object displaying the ontology tree structure.
#'
#' @importFrom igraph graph_from_data_frame V
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text
#' @importFrom ggplot2 scale_color_identity theme_void
#' @export
#'
#' @examples
#' \dontrun{
#' # Plot the entire ontology tree
#' plotOntologyTree(tree)
#'
#' # Plot a circular version of the ontology tree
#' plotOntologyTree(tree, layout = "circular")
#'
#' # Highlight specific terms in the ontology tree
#' plotOntologyTree(tree, focus_terms = c("GO:0008150", "GO:0003674"))
#' }
plotOntologyTree <- function(tree, layout = c("tree", "circular"), focus_terms = NULL) {
  layout <- match.arg(layout)
  
  # Convert OntologyTree to edge list
  edge_list <- do.call(rbind, lapply(seq_along(tree@childMap), function(i) {
    parent <- tree@termNames[i]  # Get parent name
    children <- tree@childMap[[i]]  # Get child indices
    if (length(children) > 0) {
      data.frame(from = parent, to = tree@termNames[children], stringsAsFactors = FALSE)
    } else {
      NULL
    }
  }))
  
  if (nrow(edge_list) == 0) {
    stop("The tree has no edges to plot.")
  }
  
  # Filter for specific focus_terms
  if (!is.null(focus_terms)) {
    valid_terms <- unique(c(edge_list$from, edge_list$to))
    focus_terms <- intersect(focus_terms, valid_terms)  # Keep only terms in the tree
    edge_list <- subset(edge_list, from %in% focus_terms | to %in% focus_terms)
    
    if (nrow(edge_list) == 0) {
      stop("No matching terms found in the tree.")
    }
  }
  
  # Convert to igraph object
  g <- igraph::graph_from_data_frame(edge_list, directed = TRUE)
  
  # Choose layout
  p <- ggraph::ggraph(g, layout = "dendrogram", circular = (layout == "circular"))
  
  # Highlight focus_terms in red, others in blue
  node_colors <- ifelse(V(g)$name %in% focus_terms, "red", "blue")
  
  # Generate plot
  p <- p +
    geom_edge_link(color = "gray") +
    geom_node_point(aes(color = node_colors), size = 3) +
    geom_node_text(aes(label = name), repel = TRUE, size = 4) +
    scale_color_identity() +  # Directly use color mapping
    theme_void()
  
  return(p)
}
