#' Update annotations for OntologyTree object
#' 
#' @param tree OntologyTree object
#' @param new_annot Named list of genes per term or data frame
#' @param use_background Whether to use original annotations
#' @return Updated OntologyTree object
#' @export

updateAnnotationsR <- function(tree, new_annot, use_background = FALSE) {
  if (!inherits(tree, "OntologyTree"))
    stop("'tree' must be an OntologyTree object")
  
  # Convert data frame to list if needed
  annot_list <- if (is.data.frame(new_annot)) {
    split(strsplit(new_annot$GeneID, ","), new_annot$Annot)
  } else {
    new_annot
  }
  
  # Keep only valid terms (those that exist in tree@termNames)
  valid_terms <- intersect(names(annot_list), tree@termNames)
  if (length(valid_terms) == 0)
    stop("No valid terms found")
  annot_list <- annot_list[valid_terms]
  
  # Convert each element to a character vector (if not already)
  annot_list <- lapply(annot_list, as.character)
  
  # Get term indices for valid terms using the R function term2id (assumed to return 1-based indices)
  term_indices <- term2id(tree, valid_terms)
  annot_map<-lapply(tree@annotations$list, function(x)tree@annotations$names[x])
  if (use_background) {
    # For each valid term, update the existing annotation by taking the union with the new gene IDs.
    for (i in seq_along(term_indices)) {
      idx <- term_indices[i]
      # Get the current gene IDs for the term from the tree (convert indices to gene IDs)
      old_ids <- tree@annotations$names[ annot_map[[idx]] ]
      # Compute the union of the old gene IDs and the new gene IDs
      new_ids <- union(old_ids, annot_list[[ valid_terms[i] ]])
      # Temporarily store the new gene IDs (as characters) in the mapping
      annot_map[[idx]] <- new_ids
    }
  } else {
    # For each valid term, replace the existing annotation with the new gene IDs.
    for (i in seq_along(term_indices)) {
      idx <- term_indices[i]
      annot_map[[idx]] <- annot_list[[ valid_terms[i] ]]
    }
  }
  
  # Now, re-map all gene IDs (which are currently stored as character vectors in annot_map)
  # to indices in a unified annotation set.
  all_genes <- unique(unlist(lapply(annot_map, as.character)))
  
  # Create an index mapping (a named vector): names are gene IDs, values are indices (1-based)
  gene_index_map <- setNames(seq_along(all_genes), all_genes)
  
  # Replace each element in annot_map with the corresponding indices according to gene_index_map
  annot_map <- lapply(annot_map, function(ids) unname(gene_index_map[ as.character(ids) ]))
  
  # Update the tree's annotations slots:
  tree@annotations$list <- annot_map
  tree@annotations$names <- all_genes
  
  return(tree)
}


#' Convert data frame to annotation list
#' 
#' @param df Data frame with Term and GeneID columns
#' @return Named list of gene vectors
#' @author Kai Guo
#' @export 
df2annolist <- function(df) {
  if(!all(c("Term", "GeneID") %in% colnames(df))) {
    stop("Data frame must contain 'Term' and 'GeneID' columns")
  }
  annot_list <- split(strsplit(df$GeneID, "/"), df$Term)
  lapply(annot_list, unlist)
}
#' Cluster similar terms with optional term restriction and background
#' 
#' @param tree OntologyTree object
#' @param terms Terms to analyze
#' @param new_annot Named list of genes per term or data frame with Term/GeneID columns
#' @param use_background Whether to use original annotations
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
#' @param ... Additional parameters passed to clusterST
#' @return Term clusters object
#' @author Kai Guo
#' @export
clusterSTW <- function(tree, terms, new_annot, use_background = FALSE, ...) {
  
  # Process input annotations
  annot_list <- if(is.data.frame(new_annot)) {
    split(strsplit(new_annot$GeneID, ","), new_annot$Annot)
  } else {
    new_annot
  }
  
  # Update tree annotations
  tree <- updateAnnotations(tree, annot_list, use_background)
  
  # Restrict terms if not using background
  if(!use_background) {
    available_terms <- names(annot_list)
    if(is.character(terms)) {
      terms <- intersect(terms, available_terms)
    } else {
      term_names <- tree@termNames[terms]
      valid_idx <- term_names %in% available_terms
      terms <- terms[valid_idx]
    }
    if(length(terms) == 0) stop("No valid terms found in annotations")
  }
  
  # Run clustering 
  clusterST(tree = tree, terms = terms, method = c("term", "ic", "resnik", "faith", "rel", "simic", "gogo", "wang", "anc"),
            threshold = 0.8,
            ic_method = c("annotation","universal","wang","offspring"),
            submethod = c("kappa", "overlap", "jaccard", "dice"),
            normType = "max",
            weights = NULL,
            annotUniverse= NULL, 
            useCache = TRUE,
            useIgraph = TRUE,
            verbose = TRUE,
            ...)
}
