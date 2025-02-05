#' Plot network visualization of term clusters
#' 
#' @description
#' Creates a network visualization where each cluster's parent term is connected to its children.
#' Different relation types are shown with different line styles and colors are used to distinguish clusters.
#' Nodes can be sized and shaped differently for parents vs children.
#' 
#' @param cluster_result The result object from clusterSimilarTerms
#' @param vertex_size Base size for nodes (default: 10)
#' @param show_names Whether to show term names instead of IDs (default: FALSE) 
#' @param layout Layout algorithm: "fr" (Fruchterman-Reingold) or "circle"
#' @param custom_colors Optional vector of colors for clusters
#' @param parent_shape Shape of parent nodes ("square", "circle", etc.)
#' @param child_shape Shape of child nodes
#' @param relation_types Types of relations to show (NULL shows all)
#' @param label_size Size multiplier for node labels (default: 0.7)
#' @param label_dist Distance of labels from nodes (default: 1.5)
#' @param repel_force Force of node repulsion (default: 1)
#' @param node_alpha Node color transparency (default: 0.8)
#' @param show_legend Display legend or not
#' @param ... Additional parameters passed to plot.igraph
#' @return Invisibly returns the igraph object
#' @import igraph
#' @examples
#' \dontrun{
#' # Generate clusters
#' clu <- clusterSimilarTerms(dag5, 1:30, 
#'                           method = "wang",
#'                           weights = weights,
#'                           threshold = 0.2)
#'                           
#' # Basic plot with term names
#' plotClusterNetwork(clu, show_names = TRUE)
#' 
#' # Customized plot
#' plotClusterNetwork(clu,
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
#' # Show only parent-child relations with custom colors
#' nature_palette <- c("#2166AC", "#4393C3", "#92C5DE", "#D6604D")
#' plotClusterNetwork(clu,
#'    relation_types = "parent-child",
#'    custom_colors = nature_palette,
#'    show_names = TRUE)
#' }
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
  library(igraph)
  
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
                        niter = 500,
                        area = vcount(g)^2.3,
                        repulserad = vcount(g)^2.8 * repel_force)
  } else {
    l <- layout_in_circle(g)
  }
  
  # 
  V(g)$shape <- ifelse(V(g)$type == "parent", parent_shape, child_shape)
  
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
