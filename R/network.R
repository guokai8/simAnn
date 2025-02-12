#' Plot Network Visualization of Term Clusters
#'
#' @description
#' Creates an interactive network visualization of term clusters showing hierarchical
#' relationships. The visualization displays parent-child relationships between terms
#' using different line styles for relation types and colors for cluster membership.
#' The function provides extensive customization options for node appearance, layout,
#' and labeling.
#'
#' @details
#' The network visualization represents:
#' \itemize{
#'   \item Nodes: GO terms (parents and children)
#'   \item Edges: Relationships between terms
#'   \item Colors: Different clusters
#'   \item Shapes: Distinguish between parent and child terms
#'   \item Line styles: Different types of relationships
#' }
#'
#' @param cluster_result A cluster result object from clusterST containing:
#' \itemize{
#'   \item clusters: List of term clusters
#'   \item graph: igraph object with cluster structure
#'   \item membership: Cluster assignments
#' }
#' @param vertex_size Numeric value for base node size in the plot.
#'        Larger values create bigger nodes. Default: 10.
#' @param show_names Logical indicating whether to display term names instead of
#'        GO IDs in node labels. Default: FALSE.
#' @param layout Character string specifying the layout algorithm:
#' \describe{
#'   \item{fr}{Fruchterman-Reingold force-directed layout}
#'   \item{circle}{Circular layout with nodes arranged in a circle}
#' }
#' @param custom_colors Character vector of colors for clusters. If provided,
#'        these colors will be used instead of the default palette.
#' @param parent_shape Character string specifying the shape for parent nodes
#'        (e.g., "square", "circle", "rectangle", "csquare", "crectangle").
#' @param child_shape Character string specifying the shape for child nodes.
#' @param relation_types Character vector specifying which types of relations to display.
#'        NULL shows all relation types. Options include:
#' \itemize{
#'   \item "is_a"
#'   \item "part_of"
#'   \item "regulates"
#'   \item "parent-child"
#' }
#' @param label_size Numeric multiplier for node label size. Default: 0.7.
#' @param label_dist Numeric value for distance of labels from nodes. Default: 1.5.
#' @param repel_force Numeric value for force of node repulsion in layout. Default: 1.
#' @param node_alpha Numeric value between 0 and 1 for node color transparency.
#'        Default: 0.8.
#' @param show_legend Logical indicating whether to display the plot legend.
#'        Default: TRUE.
#' @param ... Additional parameters passed to plot.igraph for fine-tuning the
#'        visualization.
#'
#' @return Invisibly returns the igraph object used for plotting, which can be
#'         used for further customization or analysis.
#'
#' @examples
#' \dontrun{
#' # Load required packages and build GO tree
#' library(org.Hs.eg.db)
#' tree <- buildGOTree(species = "human", namespace = "BP")
#'
#' # Generate clusters using Wang's method
#' weights <- c("is_a" = 0.8, "part_of" = 0.6)
#' clusters <- clusterST(tree,
#'                      terms = 1:30,
#'                      method = "wang",
#'                      weights = weights,
#'                      threshold = 0.2)
#'
#' # Basic plot with term names
#' plotClusterNetwork(clusters,
#'                   show_names = TRUE)
#'
#' # Customized visualization
#' plotClusterNetwork(clusters,
#'    show_names = TRUE,
#'    layout = "fr",
#'    vertex_size = 12,
#'    label_size = 0.6,
#'    label_dist = 2,
#'    repel_force = 1.2,
#'    node_alpha = 0.85,
#'    parent_shape = "square",
#'    child_shape = "circle")
#'
#' # Use custom color palette
#' nature_colors <- c("#2166AC", "#4393C3", "#92C5DE", "#D6604D")
#' plotClusterNetwork(clusters,
#'    relation_types = c("is_a", "part_of"),
#'    custom_colors = nature_colors,
#'    show_names = TRUE,
#'    show_legend = TRUE)
#'
#' # Focus on specific relation types
#' plotClusterNetwork(clusters,
#'    relation_types = "is_a",
#'    vertex_size = 15,
#'    show_names = TRUE,
#'    layout = "circle")
#' }
#'
#' @seealso
#' \code{\link{clusterST}} for generating term clusters
#' \code{\link[igraph]{plot.igraph}} for additional plotting parameters
#'
#' @references
#' Csardi G, Nepusz T (2006). The igraph software package for complex network
#' research. InterJournal, Complex Systems, 1695.
#'
#' @importFrom igraph plot.igraph V E
#' @return a igraph figure
#' @author Kai Guo
#' @export
plotClusterNetwork <- function(cluster_result,
                               vertex_size = 10,
                               show_names = FALSE,
                               layout = c("fr", "circle"), 
                               custom_colors = NULL,
                               parent_shape = "square",
                               child_shape = "circle",
                               relation_types = NULL,
                               label_size = 0.7,
                               label_dist = 1.5,
                               repel_force = 1,
                               node_alpha = 0.8,
                               show_legend = FALSE,
                               ...,
                               plot_genes = FALSE,    # 
                               top_n_genes = 5,       # 
                               selected_genes = NULL  # 
) {
  require("igraph")  
  # 
  nature_colors <- c("#2166AC", "#4393C3", "#92C5DE", "#D6604D", 
                     "#B2182B", "#053061", "#67001F", "#F4A582",
                     "#67001F", "#F4A582", "#762A83","#A6761D","#D95F02","deepskyblue","#1B9E77","#E7298A","#7570B3","#E6AB02",
                     "#FC9D9A","#83AF9B","#FE4365","#F9CDAD","#C8C8A9","#B6C29A","#8A977B",
                     "#F4D000","#E58308","#DC5712","#A6CEE3","#1F78B4","#B2DF8A","#33A02C",
                     "slateblue1","darkgreen","darkred","plum1","#FB9A99","#E31A1C",
                     "darkmagenta","hotpink2","slategray4","magenta2","yellow4",
                     "#66A61E","#E7298A","#7570B3","#E6AB02",'#e6194b','#3cb44b','#ffe119',
                     '#4363d8','#f58231','#911eb4','#46f0f0','#f032e6','#bcf60c',
                     '#fabebe','#008080','#e6beff','#9a6324','#fffac8','#800000',
                     '#aaffc3','#808000','#ffd8b1','#000075','#808080')
  
  layout <- match.arg(layout)
  # 
  filtered_clusters <- if(!is.null(relation_types)) {
    subset <- cluster_result$cluster_matrix$parent %in% 
      cluster_result$clusters[sapply(cluster_result$clusters, 
                                     function(x) x$relation_type %in% relation_types)]$cluster_term
    cluster_result$cluster_matrix[subset, ]
  } else {
    cluster_result$cluster_matrix
  }
  
  # 
  if(!plot_genes) {
    
    all_terms <- unique(c(filtered_clusters$parent, filtered_clusters$child))
    nodes <- data.frame(
      name = all_terms,
      # 
      type = ifelse(all_terms %in% filtered_clusters$parent,
                    "parent","child"),
      stringsAsFactors = FALSE
    )
    
    # 
    edges <- data.frame(
      from = filtered_clusters$parent,
      to = filtered_clusters$child,
      weight = abs(filtered_clusters$similarity) + 0.1,  # 
      relation_type = filtered_clusters$relation_type,
      stringsAsFactors = FALSE
    )
    
  } else {

    cm <- filtered_clusters
    
    # 
    cm <- cm[!is.na(cm$parent) & cm$parent!="", ]
    
    # 
    parent2genes <- tapply(cm$GeneID, cm$parent, function(x){

      unq <- unique(unlist(lapply(x, function(s){
        gs <- strsplit(s, "[,;]")[[1]]
        trimws(gs)
      })))
      unq[unq!=""]
    })
    

    filter_genes <- function(gvec) {
      if(!is.null(selected_genes)) {
        # 
        gvec <- intersect(gvec, selected_genes)
      }
      # 
      if(length(gvec) > top_n_genes) {
        gvec <- gvec[1:top_n_genes]
      }
      gvec
    }
    
    parent2genes <- lapply(parent2genes, filter_genes)
    
    # 
    keep_parents <- names(parent2genes)[sapply(parent2genes, length) > 0]
    parent2genes <- parent2genes[keep_parents]
    
    # 
    parent_nodes <- keep_parents
    gene_nodes <- unique(unlist(parent2genes))
    
    all_terms <- c(parent_nodes, gene_nodes)
    nodes <- data.frame(
      name = all_terms,
      type = ifelse(all_terms %in% parent_nodes, "parent", "child"),
      stringsAsFactors = FALSE
    )
    
    # 
    edges_list <- list()
    idx <- 1
    for(p in parent_nodes) {
      gset <- parent2genes[[p]]
      if(length(gset) > 0) {
        edges_list[[idx]] <- data.frame(
          from = rep(p, length(gset)),
          to   = gset,
          weight = 1,  # 
          relation_type = "parent-gene",
          stringsAsFactors = FALSE
        )
        idx <- idx + 1
      }
    }
    edges <- do.call(rbind, edges_list)
    
    if(nrow(edges)==0) {
      warning("No edges (no genes) to display under the given top_n_genes / selected_genes filter.")
      return(invisible(NULL))
    }
  }
  
  # 
  if(nrow(edges) == 0) {
    warning("No edges to display for the specified relation types or no valid data.")
    return(invisible(NULL))
  }
  
  # 
  g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
  
  # 
  V(g)$size <- ifelse(V(g)$type == "parent", vertex_size * 1.5, vertex_size)
  
  # 
  if(show_names) {
    id_to_name <- c(
      setNames(cluster_result$cluster_matrix$parent_name,
               cluster_result$cluster_matrix$parent),
      setNames(cluster_result$cluster_matrix$child_name,
               cluster_result$cluster_matrix$child)
    )
    id_to_name <- id_to_name[!duplicated(names(id_to_name))]
    
    V(g)$label <- sapply(V(g)$name, function(n) {
      if(n %in% names(id_to_name)) {
        if(!is.na(id_to_name[n])) {
          id_to_name[n]
        } else {
          n
        }
      } else {
        n
      }
    })
  } else {
    V(g)$label <- V(g)$name
  }
  
  # 
  if(is.null(custom_colors)) {
    # 
    uniq_parents <- unique(edges$from)
    cluster_colors <- nature_colors[seq_along(uniq_parents)]
    names(cluster_colors) <- uniq_parents
  } else {
    #
    uniq_parents <- unique(edges$from)
    cluster_colors <- rep(custom_colors, length.out = length(uniq_parents))
    names(cluster_colors) <- uniq_parents
  }
  
  # 
  V(g)$color <- sapply(V(g)$name, function(n) {
    if(n %in% names(cluster_colors)) {
      # parent
      adjustcolor(cluster_colors[n], alpha.f = node_alpha)
    } else {
      # child
      parent <- edges$from[edges$to == n]
      if(length(parent) > 0) parent <- parent[1]
      if(!is.na(parent) && parent %in% names(cluster_colors)) {
        adjustcolor(cluster_colors[parent], alpha.f = node_alpha * 0.7)
      } else {
        adjustcolor("gray80", alpha.f = node_alpha)
      }
    }
  })
  
  # 
  E(g)$color <- sapply(seq_len(nrow(edges)), function(i) {
    p <- edges$from[i]
    if(p %in% names(cluster_colors)) {
      adjustcolor(cluster_colors[p], alpha.f = 0.3)
    } else {
      "gray80"
    }
  })
  
  # 
  E(g)$lty <- sapply(edges$relation_type, function(rt) {
    switch(rt,
           "parent-child" = 1,
           "similarity-based" = 2,
           "singleton" = 3,
           "parent-gene" = 1,
           1)
  })
  
  # 
  E(g)$width <- ifelse("weight" %in% colnames(edges), edges$weight*2, 1)
  
  # 
  if(layout == "fr") {
    l <- layout_with_fr(g,
                        weights = E(g)$weight,
                        niter = 500)
  } else {
    l <- layout_in_circle(g)
  }
  
  # 
  V(g)$shape <- ifelse(V(g)$type == "parent", parent_shape, child_shape)
  g <- igraph::simplify(g)
  plot(g,
       layout = l * 1.5,
       vertex.shape = V(g)$shape,
       vertex.frame.color = adjustcolor("gray40", alpha.f = 0.3),
       vertex.label.dist = label_dist,
       vertex.label.color = "black", 
       vertex.label.family = "sans",
       vertex.label.cex = label_size,
       edge.curved = FALSE,
       margin = c(0, 0, 0, 0),
       ...)
  
  if(show_legend && length(unique(edges$relation_type)) > 1) {
    legend("bottomright",
           legend = unique(edges$relation_type),
           lty = sapply(unique(edges$relation_type), function(rt) {
             switch(rt,
                    "parent-child" = 1,
                    "similarity-based" = 2,
                    "singleton" = 3,
                    "parent-gene" = 1,
                    1)
           }),
           title = "Relation Types",
           cex = 0.8,
           bg = adjustcolor("white", alpha.f = 0.8))
  }
  
  invisible(g)
}
