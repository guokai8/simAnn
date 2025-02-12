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
    lapply(split(strsplit(new_annot$GeneID, ","), new_annot$Annot),unlist)
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

#' Cluster Similar Terms with Wang's Method and Weighted Edges
#'
#' This function clusters similar terms using Wang's semantic similarity method,
#' supporting weighted edge types in the ontology graph. It can analyze either
#' direct gene sets or enrichment analysis results.
#'
#' @param tree An OntologyTree object containing the Gene Ontology hierarchy.
#' @param terms Character vector of GO terms to analyze.
#' @param new_annot Either a named list of gene sets per term, or a data frame 
#'        containing enrichment results with Term and GeneID columns.
#' @param use_background Logical. Whether to use original annotations from the ontology.
#' @param method Character string specifying the similarity calculation method:
#' \describe{
#'   \item{term}{Term overlap based methods (kappa, overlap, jaccard, dice)}
#'   \item{ic}{Universal information content}
#'   \item{resnik}{Resnik's semantic similarity}
#'   \item{lin}{Lin's semantic similarity}
#'   \item{faith}{FaITH semantic similarity}
#'   \item{rel}{Relevance similarity}
#'   \item{simic}{SimIC similarity}
#'   \item{wang}{Wang's semantic similarity method}
#'   \item{anc}{Ancestor-based similarity}
#' }
#' @param ic_method Character string specifying IC calculation method for faith/rel/simic:
#' \describe{
#'   \item{universal}{Universal-based Information Content}
#'   \item{annotation}{Annotation-based Information Content}
#'   \item{wang}{Wang's Information Content}
#'   \item{offspring}{Offspring-based Information Content}
#' }
#' @param submethod For "term" method, specifies similarity coefficient:
#' \describe{
#'   \item{kappa}{Cohen's kappa coefficient}
#'   \item{overlap}{Overlap coefficient}
#'   \item{jaccard}{Jaccard similarity}
#'   \item{dice}{Dice coefficient}
#' }
#' @param normType For Resnik similarity, specifies normalization:
#' \describe{
#'   \item{unif}{By log of max annotation count}
#'   \item{max}{By maximum IC value}
#'   \item{sum}{By maximum pairwise IC}
#'   \item{none}{No normalization}
#' }
#' @param weights Named numeric vector of weights for edge types (e.g., is_a, part_of).
#'        Values should be < 1.
#' @param useIgraph Logical. Whether to use igraph-based calculation.
#' @param useCache Logical. Whether to use cached IC values.
#' @param verbose Logical. Whether to print progress messages.
#' @param ... Additional parameters passed to underlying clustering function.
#'
#' @return A term clusters object containing the clustering results.
#'
#' @examples
#' # Example 1: Using gene sets
#' # Define edge weights
#' weights <- c("is_a" = 0.8, "part_of" = 0.6)
#' 
#' # Define gene sets
#' geneset <- list(
#'   "GO:0000086" = c("TAF2", "FOXM1", "RRM1", "ATR"),
#'   "GO:0000422" = c("LRBA", "EIF2AK1", "SPATA18"),
#'   "GO:0001678" = c("SYBU", "KCNB1", "PLA2G6", "LIN28A"),
#'   "GO:0001704" = c("SNAI1", "GPI", "AHDC1"),
#'   "GO:0006282" = c("TAF2", "FOXM1", "BRD8", "ATR", "PARPBP")
#' )
#' 
#' # Cluster terms using Wang's method
#' clu <- clusterSTW(tree, 
#'                   terms = names(geneset),
#'                   new_annot = geneset,
#'                   method = "wang",
#'                   weights = weights,
#'                   threshold = 0.2)
#'
#' # Example 2: Using enrichment results
#' # Assuming 'res' is enrichment results from richR
#' \dontrun{
#' clu <- clusterSTW(tree,
#'                   terms = res$Annot,
#'                   new_annot = as.data.frame(res),
#'                   method = "wang",
#'                   weights = weights,
#'                   threshold = 0.2)
#' }
#'
#' @references
#' Wang et al. (2007) A new method to measure the semantic similarity of GO terms.
#' Bioinformatics, 23(10):1274-1281. \doi{10.1093/bioinformatics/btm087}
#'
#' @seealso
#' \code{\link{clusterST}} for basic term clustering
#'
#' @importFrom stats p.adjust
#' @author Kai Guo
#' @export
#'
clusterSTW <- function(tree, terms, new_annot, method = c("term", "ic", "resnik", "lin","faith", "rel", "simic", "gogo", "wang", "anc"), 
                       threshold = 0.5,
                       ic_method = c("annotation","universal","wang","offspring"),
                       submethod = c("kappa", "overlap", "jaccard", "dice"),
                       normType = "max",
                       weights = NULL,
                       annotUniverse= NULL, 
                       useCache = TRUE,
                       useIgraph = TRUE,
                       verbose = TRUE,
                       use_background = FALSE, ...) {
  
  method <- match.arg(method)
  ic_method <- match.arg(ic_method)
  submethod <- match.arg(submethod)
  # Process input annotations
  annot_list <- if(is.data.frame(new_annot)) {
    lapply(split(strsplit(new_annot$GeneID, ","), new_annot$Annot),unlist)
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
  clusterST(tree = tree, terms = terms, method = method,
            threshold = threshold,
            ic_method = ic_method,
            submethod = submethod,
            normType = normType,
            weights = weights,
            annotUniverse= annotUniverse, 
            useCache = useCache,
            useIgraph = useIgraph,
            verbose = verbose,
            ...)
}

#' Calculate semantic similarity between terms with optional term restriction and background
#' Calculate Weighted Semantic Similarity Between Terms
#'
#' @description
#' A unified interface for calculating semantic similarity between ontology terms using
#' various established methods with support for weighted relationships. The function supports
#' multiple similarity calculation approaches including term overlap, information content
#' based methods, and graph-based methods. It can work with both standard GO annotations
#' and custom gene sets.
#'
#' @details
#' The function implements several semantic similarity calculation methods with weighted
#' edge support:
#' \itemize{
#'   \item Term-based methods use direct comparison of gene sets
#'   \item Information Content (IC) based methods use the specificity of terms
#'   \item Graph-based methods use the ontology structure with edge weights
#'   \item Hybrid methods combine multiple approaches
#' }
#'
#' @param tree An OntologyTree object containing the ontology structure and annotations.
#' @param terms Character vector of terms between which to calculate similarity.
#' @param new_annot Either a named list mapping terms to gene sets, or a data frame 
#'        containing Term and GeneID columns from enrichment analysis.
#' @param use_background Logical indicating whether to use original annotations from
#'        the ontology tree. Default is FALSE.
#' @param method Character string specifying the similarity calculation method:
#' \describe{
#'   \item{term}{Term overlap based methods using gene set comparison}
#'   \item{ic}{Universal information content using term specificity}
#'   \item{resnik}{Resnik's information content based similarity
#'                 (\doi{10.1613/jair.514}, \doi{10.1186/1471-2105-9-S5-S4})}
#'   \item{lin}{Lin's semantic similarity (\doi{10.5555/645527.657297})}
#'   \item{faith}{FaITH semantic similarity (\doi{10.1007/978-3-642-17746-0_39})}
#'   \item{rel}{Relevance similarity (\doi{10.1186/1471-2105-7-302})}
#'   \item{simic}{SimIC similarity (\doi{10.48550/arXiv.1001.0958})}
#'   \item{gogo}{GOGO similarity (\doi{10.1038/s41598-018-33219-y})}
#'   \item{wang}{Wang's semantic similarity (\doi{10.1093/bioinformatics/btm087})}
#'   \item{anc}{Ancestor based similarity using common ancestors}
#' }
#'
#' @param ic_method Character string specifying IC calculation method for 
#'        faith, rel, and simic methods:
#' \describe{
#'   \item{universal}{Universal-based IC using term frequencies}
#'   \item{annotation}{Annotation-based IC using gene annotations}
#'   \item{wang}{Wang's IC calculation method}
#'   \item{offspring}{Offspring-based IC using term descendants}
#' }
#'
#' @param submethod For "term" method, specifies the set similarity measure:
#' \describe{
#'   \item{kappa}{Cohen's kappa coefficient for categorical agreement}
#'   \item{overlap}{Overlap coefficient: |A∩B| / min(|A|,|B|)}
#'   \item{jaccard}{Jaccard similarity: |A∩B| / |A∪B|}
#'   \item{dice}{Dice coefficient: 2|A∩B| / (|A|+|B|)}
#' }
#'
#' @param normType For Resnik similarity, specifies normalization method:
#' \describe{
#'   \item{unif}{Normalize by log of maximum annotation count}
#'   \item{max}{Normalize by maximum IC value}
#'   \item{sum}{Normalize by maximum pairwise IC}
#'   \item{none}{No normalization}
#' }
#'
#' @param weights Named numeric vector of weights for edge types in GOGO/Wang methods.
#'        Names must match relation types (e.g., "is_a", "part_of"), values should be < 1.
#' @param useIgraph Logical indicating whether to use igraph package for calculations.
#'        Recommended for large graphs. Default is FALSE.
#' @param useCache Logical indicating whether to use cached IC values when applicable.
#'        Default is TRUE.
#' @param verbose Logical indicating whether to print progress messages.
#'        Default is FALSE.
#' @param ... Additional parameters passed to internal functions.
#'
#' @return A matrix of pairwise semantic similarity scores between terms.
#'
#' @examples
#' # Load example data
#' data(goTree)
#' 
#' # Define edge weights
#' weights <- c("is_a" = 0.8, "part_of" = 0.6)
#' 
#'
#' # Example 1: Using custom gene sets
#' geneset <- list(
#'   "GO:0000086" = c("TAF2", "FOXM1", "RRM1", "ATR"),
#'   "GO:0000422" = c("LRBA", "EIF2AK1", "SPATA18"),
#'   "GO:0001678" = c("SYBU", "KCNB1", "PLA2G6", "LIN28A")
#' )
#' sim2 <- simtermW(goTree$human, names(geneset), 
#'                  new_annot = geneset,
#'                  method = "wang",
#'                  weights = weights)
#'
#' # Example 2: Using enrichment results
#' \dontrun{
#' # Assuming 'enrich_res' is from enrichment analysis
#' sim3 <- simtermW(goTree$human, 
#'                  enrich_res$Term,
#'                  new_annot = enrich_res,
#'                  method = "resnik",
#'                  normType = "max",
#'                  weights = weights)
#' }
#'
#' @references
#' \itemize{
#'   \item Resnik P (1995). Using information content to evaluate semantic similarity in a taxonomy.
#'   \item Lin D (1998). An information-theoretic definition of similarity.
#'   \item Wang JZ, et al (2007). A new method to measure the semantic similarity of GO terms.
#'   \item Yang H, et al (2012). GOGO: An improved algorithm to measure the semantic similarity.
#' }
#'
#' @seealso
#' \code{\link{simterm}} for calculating semantic similarity between terms
#' @importFrom stats cor
#' @importFrom methods is
#' @return Term clusters object
#' @author Kai Guo
#' @export
simtermW <- function(tree, terms, new_annot, use_background = FALSE, ...) {
  
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
  simterm(tree = tree, terms = terms, method = c("term", "ic", "resnik", "lin","faith", "rel", "simic", "gogo", "wang", "anc"),
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
