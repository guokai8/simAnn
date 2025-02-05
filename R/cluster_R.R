#' Process similar terms clustering
#' 
#' @description
#' Groups semantically similar terms based on similarity scores and parent-child relationships.
#' Returns results in multiple formats as a list.
#' 
#' @param tree An OntologyTree object
#' @param terms A character vector of term IDs/names, or numeric indices
#' @param method Similarity method to use:
#'   * "term": Term overlap based methods
#'   * "ic": Universal information content
#'   * "resnik": Resnik's information content
#'   * "faith": FaITH similarity
#'   * "rel": Relevance similarity
#'   * "simic": SimIC similarity
#'   * "gogo": GOGO similarity
#'   * "wang": Wang's method
#'   * "anc": Ancestor based similarity
#' @param threshold Similarity threshold for clustering (default: 0.8)
#' @param ic_method For "faith","rel",and "simic" methods:
#'   * "universal": Calculate universal-based Information Content
#'   * "annotation": Calculate annotation-based Information Content
#'   * "wang": Wang's Information Content
#'   * "offspring": Calculate offspring-based Information Content
#' @param submethod For "term" method only:
#'   * "kappa": Cohen's kappa
#'   * "overlap": Overlap coefficient
#'   * "jaccard": Jaccard similarity
#'   * "dice": Dice coefficient
#' @param normType For Resnik similarity only (unif/max/sum/none)
#' @param weights Named vector of weights for GOGO/Wang methods
#' @param normType For Resnik similarity only:
#'   * "unif": By log of max annotation count
#'   * "max": By maximum IC value
#'   * "sum": By maximum pairwise IC
#'   * "none": No normalization
#' @param useCache Use cached IC values when applicable
#' @param useIgraph Use igraph-based calculation
#' @param verbose Print progress messages
#' @return A list containing multiple representations of clustering results
#' @export
clusterSimilarTerms <- function(
    tree,
    terms,
    method = c("term", "ic", "resnik", "faith", "rel", "simic", "gogo", "wang", "anc"),
    threshold = 0.8,
    ic_method = c("annotation","universal","wang","offspring"),
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
  ic_method <- match.arg(ic_method)
  if(!inherits(tree, "OntologyTree")) {
    stop("'tree' must be an OntologyTree object")
  }
  
  # Convert terms to indices if necessary
  terms_vec <- if(is.character(terms)) {
    term2id(tree, terms, strict = FALSE)
  } else {
    terms
  }
  validTerms <- validateAnnotatedTerms(tree, terms_vec)
  terms_vec <- terms_vec[validTerms]
  if(length(terms_vec) == 0) {
    stop("No valid terms provided")
  }
  
  # Get term names from metadata
  term_names <- tree@metadata$name
  if(is.null(term_names)) {
    warning("Term names not found in tree metadata")
    term_names <- tree@termNames
  }
  names(term_names) <- tree@metadata$id
  
  # Calculate similarity matrix
  sim_matrix <- simterm(
    tree = tree,
    terms = terms_vec,
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
  
  # Initialize results
  parent_clusters <- list()
  singleton_clusters <- list()
  processed_terms <- integer(0)
  
  # Process terms and build clusters
  for(i in seq_along(terms_vec)) {
    current_term_id <- terms_vec[i]
    if(current_term_id %in% processed_terms) next
    
    current_term <- tree@termNames[current_term_id]
    if(current_term %in% names(parent_clusters)) next
    
    similarities <- sim_matrix[i, ]
    similar_indices <- which(similarities >= threshold & seq_along(similarities) != i)
    
    if(length(similar_indices) > 0) {
      similar_term_ids <- terms_vec[similar_indices]
      
      all_term_ids <- c(current_term_id, similar_term_ids)
      parent_found <- FALSE
      parent_term_id <- NULL
      
      for(term1 in all_term_ids) {
        for(term2 in all_term_ids) {
          if(term1 != term2 && term2 %in% tree@parentMap[[term1]]) {
            parent_term_id <- term2
            parent_found <- TRUE
            break
          }
        }
        if(parent_found) break
      }
      
      if(!parent_found) {
        ic_scores <- sapply(all_term_ids, function(id) 
          calculateIC(tree, method = ic_method, verbose = verbose)[id])
        parent_term_id <- all_term_ids[which.min(ic_scores)]
      }
      
      child_term_ids <- setdiff(all_term_ids, parent_term_id)
      child_similarities <- sim_matrix[match(parent_term_id, terms_vec), 
                                       match(child_term_ids, terms_vec)]
      
      parent_term <- tree@termNames[parent_term_id]
      child_terms <- tree@termNames[child_term_ids]
      
      # Extract GeneIDs for the cluster
      cluster_gene_ids <- extract_genes_from_terms(tree, all_term_ids)
      
      if(parent_term %in% names(parent_clusters)) {
        existing_cluster <- parent_clusters[[parent_term]]
        new_children <- setdiff(child_terms, existing_cluster$children)
        
        if(length(new_children) > 0) {
          new_similarities <- sim_matrix[match(parent_term_id, terms_vec),
                                         match(new_children, colnames(sim_matrix))]
          
          all_children <- c(existing_cluster$children, new_children)
          all_similarities <- c(existing_cluster$child_similarities, new_similarities)
          
          parent_clusters[[parent_term]] <- list(
            cluster_term = parent_term,
            mean_similarity = mean(all_similarities),
            max_similarity = max(all_similarities),
            relation_type = if(parent_found) "parent-child" else existing_cluster$relation_type,
            children = all_children,
            term_name = term_names[parent_term],
            child_similarities = all_similarities,
            child_names = term_names[all_children],
            GeneID = unique(c(existing_cluster$GeneID, cluster_gene_ids))
          )
        }
      } else {
        parent_clusters[[parent_term]] <- list(
          cluster_term = parent_term,
          mean_similarity = mean(child_similarities),
          max_similarity = max(child_similarities),
          relation_type = if(parent_found) "parent-child" else "similarity-based",
          children = child_terms,
          term_name = term_names[parent_term],
          child_similarities = child_similarities,
          child_names = term_names[child_terms],
          GeneID = cluster_gene_ids
        )
      }
      
      processed_terms <- c(processed_terms, all_term_ids)
      
    } else {
      current_term <- tree@termNames[current_term_id]
      singleton_gene_ids <- extract_genes_from_terms(tree, current_term_id)
      
      singleton_clusters[[current_term]] <- list(
        cluster_term = current_term,
        mean_similarity = 0,
        max_similarity = 0,
        relation_type = "singleton",
        children = current_term,
        term_name = term_names[current_term],
        child_similarities = 0,
        child_names = term_names[current_term],
        GeneID = singleton_gene_ids
      )
      
      processed_terms <- c(processed_terms, current_term_id)
    }
  }
  
  clusters <- list()
  
  for (name in names(parent_clusters)) {
    clusters[[length(clusters) + 1]] <- parent_clusters[[name]]
  }
  
  for (name in names(singleton_clusters)) {
    clusters[[length(clusters) + 1]] <- singleton_clusters[[name]]
  }
  
  is_singleton <- sapply(clusters, function(x) x$relation_type == "singleton")
  mean_sims <- sapply(clusters, function(x) if(is.na(x$mean_similarity)) -1 else x$mean_similarity)
  parent_child <- sapply(clusters, function(x) x$relation_type == "parent-child")
  
  order_idx <- order(is_singleton, -parent_child, -mean_sims)
  clusters <- clusters[order_idx]
  
  summary_df <- data.frame(
    cluster_term = vapply(clusters, function(x) x$cluster_term, character(1)),
    mean_similarity = vapply(clusters, function(x) x$mean_similarity, numeric(1)),
    max_similarity = vapply(clusters, function(x) x$max_similarity, numeric(1)),
    relation_type = vapply(clusters, function(x) x$relation_type, character(1)),
    children = I(lapply(clusters, function(x) x$children)),
    GeneID = vapply(clusters, function(x) paste(x$GeneID, collapse = ","), character(1)),
    term_name = vapply(clusters, function(x) x$term_name, character(1)),
    stringsAsFactors = FALSE
  )
  
  cluster_rows <- vector("list", length = 0)
  for(x in clusters) {
    if(x$relation_type == "singleton") {
      cluster_rows[[length(cluster_rows) + 1]] <- c(
        parent = x$cluster_term,
        parent_name = x$term_name,
        child = x$children[1],
        child_name = x$child_names[1],
        similarity = 0,
        relation_type = x$relation_type,
        GeneID = paste(x$GeneID, collapse = ",")
      )
    } else if(length(x$children) > 0) {
      for(i in seq_along(x$children)) {
        cluster_rows[[length(cluster_rows) + 1]] <- c(
          parent = x$cluster_term,
          parent_name = x$term_name,
          child = x$children[i],
          child_name = x$child_names[i],
          similarity = x$child_similarities[i],
          relation_type = x$relation_type,
          GeneID = paste(x$GeneID, collapse = ",")
        )
      }
    }
  }
  
  cluster_matrix <- if(length(cluster_rows) > 0) {
    result <- as.data.frame(do.call(rbind, cluster_rows), stringsAsFactors = FALSE)
    result$similarity <- as.numeric(result$similarity)
    rownames(result) <- NULL
    result
  } else {
    result <- data.frame(
      parent = character(0),
      parent_name = character(0),
      child = character(0),
      child_name = character(0),
      similarity = numeric(0),
      relation_type = character(0),
      GeneID = character(0),
      stringsAsFactors = FALSE
    )
    rownames(result) <- NULL
    result
  }
  
  stats <- list(
    total_clusters = length(clusters),
    parent_child_clusters = sum(sapply(clusters, function(x) x$relation_type == "parent-child")),
    similarity_clusters = sum(sapply(clusters, function(x) x$relation_type == "similarity-based")),
    singleton_terms = sum(is_singleton),
    mean_similarity = mean(sapply(clusters, function(x) x$mean_similarity), na.rm = TRUE)
  )
  
  result <- list(
    clusters = clusters,
    summary = summary_df,
    cluster_matrix = cluster_matrix,
    similarity_matrix = sim_matrix,
    stats = stats
  )
  
  attr(result, "method") <- method
  attr(result, "threshold") <- threshold
  attr(result, "timestamp") <- Sys.time()
  
  class(result) <- c("term_clusters", "list")
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

#' Plot similarity heatmap based on clustering
#' 
#' @description
#' Creates a heatmap visualization of term similarities, organized by clustering order
#' 
#' @param cluster_result The result object from clusterSimilarTerms
#' @param show_names Logical, whether to show term names instead of IDs (default: FALSE)
#' @param annotation_style One of "combined" or "separate" for cluster name display
#' @param cex Font size for labels (default: 0.8)
#' @param margins Numeric vector of length 2 for plot margins (default: c(8, 8))
#' @param colors Color palette for the heatmap
#' @return Invisibly returns the reordered similarity matrix
#' @importFrom stats hclust
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

#' Print method for term_clusters objects
#' @export
print.term_clusters <- function(x, ...) {
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
