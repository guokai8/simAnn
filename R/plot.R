#' @title Plot Ontology Tree
#' @description Visualizes an OntologyTree as a hierarchical tree or circular dendrogram.
#'
#' @param tree An OntologyTree object representing the hierarchical ontology structure.
#' @param layout Character, either "tree" (default) for a hierarchical dendrogram.
#' @param focus_terms Optional character vector specifying terms to highlight in the tree. If NULL, the entire tree is plotted.
#'
#' @return A ggplot object displaying the ontology tree structure.
#'
#' @importFrom igraph graph_from_data_frame V degree distances
#' @importFrom ggraph ggraph create_layout geom_edge_link geom_node_point geom_node_text
#' @importFrom ggplot2 theme_void
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

