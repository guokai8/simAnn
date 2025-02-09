#' Calculate semantic similarity between terms
#' 
#' @description 
#' A unified interface for calculating semantic similarity between ontology terms using
#' various methods. Each method has its own characteristics and is suitable for
#' different scenarios.
#' 
#' @details
#' Available methods include:
#' 
#' ## Term Overlap Methods (method = "term")
#' 
#' * kappa: Cohen's kappa coefficient
#' 
#' * overlap: Overlap coefficient
#' 
#' * jaccard: Jaccard similarity coefficient
#' 
#' * dice: Dice coefficient
#' 
#' 
#' ## Information Content Based Methods
#' 
#' * ic: Universal IC that only depends on the structure
#' 
#' * resnik: Resnik's method with multiple normalization options
#' 
#' * lin: Lin's similarity method
#'
#' * faith: FaITH similarity for sparse annotations
#' 
#' * rel: Relevance similarity with IC-based correction
#' 
#' * simic: SimIC similarity with improved correction
#' 
#' 
#' ## Graph-based Methods
#' * gogo: GOGO method optimized for GO terms
#' 
#' * wang: Wang's method considering relation types
#' 
#' * anc: Simple ancestor-based Jaccard similarity
#' 
#' 
#' @param tree An OntologyTree object
#' @param terms A character vector of term IDs/names, or numeric indices
#' @param method Similarity method to use:
#' 
#' * "term": Term overlap based methods
#' 
#' * "ic": Universal information content
#' 
#' * "resnik": Resnik's information content 
#' Link: 
#' \doi{10.1613/jair.514}, \doi{10.1186/1471-2105-9-S5-S4}, \doi{10.1186/1471-2105-11-562}, \doi{10.1155/2013/292063}.
#' 
#' * "lin": Lin similarity
#' Link: \doi{10.5555/645527.657297}.
#' 
#' * "faith": FaITH similarity
#' Link: \doi{10.1007/978-3-642-17746-0_39}.
#' 
#' * "rel": Relevance similarity
#' Link: \doi{10.1186/1471-2105-7-302}.
#' 
#'* "simic": SimIC similarity
#' Link: \doi{10.48550/arXiv.1001.0958}.
#' 
#' * "gogo": GOGO similarity
#' Link: \doi{10.1038/s41598-018-33219-y}.
#' 
#' * "wang": Wang's method
#' Link: \doi{10.1093/bioinformatics/btm087}.
#' 
#' * "anc": Ancestor based similarity
#'
#' @param ic_method For "faith","rel",and "simic" methods:
#'   * "universal": Calculate universal-based Information Content
#'   * "annotation": Calculate annotation-based Information Content
#'   * "wang": Wang's Information Content
#'   * "offspring": Calculate offspring-based Information Content
#' @param submethod For "term" method only:
#' * "kappa": Cohen's kappa
#' 
#' Denote two sets `A` and `B` as the items annotated to term `a` and `b`. The similarity value is [the kappa coeffcient](https://en.wikipedia.org/wiki/Cohen%27s_kappa)
#' of the two sets. 
#' 
#' * "overlap": Overlap coefficient
#' 
#' Denote two sets `A` and `B` as the items annotated to terms `a` and `b`. The similarity value is the overlap coefficient
#' of the two sets, defined as `length(intersect(A, B))/min(length(A), length(B))`.
#' 
#' * "jaccard": Jaccard similarity
#' 
#' Denote two sets `A` and `B` as the items annotated to term `a` and `b`. The similarity value is the Jaccard coeffcient
#' of the two sets, defined as `length(intersect(A, B))/length(union(A, B))`.   
#' 
#' * "dice": Dice coefficient
#' 
#' Denote two sets `A` and `B` as the items annotated to term `a` and `b`. The similarity value is the Dice coeffcient
#' of the two sets, defined as `2*length(intersect(A, B))/(length(A) + length(B))`.
#'   
#' @param normType For Resnik similarity only:
#'   * "unif": By log of max annotation count
#'   
#'   * "max": By maximum IC value
#'   
#'   * "sum": By maximum pairwise IC
#'   
#'   * "none": No normalization
#'   
#' @param weights Named vector of weights for GOGO/Wang methods. 
#'    Names must match relation types (e.g., "is_a", "part_of"), values < 1.
#' @param useIgraph Use igraph-based calculation (faster for large graphs)
#' @param useCache Use cached IC values when applicable
#' @param verbose Print progress messages
#' @return A symmetric numeric matrix of similarity scores with row/column names
#' @examples
#' \dontrun{
#' # Basic usage
#' tree <-  buildGOTree(namespace = "BP", orgDb = "org.Hs.eg.db")
#' sim <- simterm(tree, c("GO:0008150", "GO:0000011"), method = "ic")
#' 
#' # Term overlap similarity
#' sim <- simterm(tree, terms, method = "term", submethod = "kappa")
#' 
#' # Resnik with normalization
#' sim <- simterm(tree, terms, method = "resnik", normType = "max")
#' 
#' # Wang method with custom weights
#' weights <- c("is_a" = 0.8, "part_of" = 0.6)
#' sim <- simterm(tree, terms, method = "wang", weights = weights)
#' 
#' # GOGO method with basic calculation
#' sim <- simterm(tree, terms, method = "gogo", useIgraph = FALSE)
#' }
#' @seealso 
#' Individual similarity functions:
#' * [calculateTermSimilarity()]
#' * [calculateUniversalSim()]
#' * [calculateResnikSim()]
#' * [calculateWangSim()]
#' @author Kai Guo
#' @export
simterm <- function(tree, terms, 
                  method = c("term", "ic", "resnik", "lin", "faith", "rel", "simic", 
                             "gogo", "wang", "anc"),
                  ic_method = NULL,
                  submethod = c("kappa", "overlap", "jaccard", "dice"),
                  normType = "max",
                  annotUniverse= NULL,
                  weights = NULL,
                  useIgraph = TRUE, 
                  useCache = TRUE,
                  verbose = TRUE) {
  # Input validation
  method <- match.arg(method)
  if(!inherits(tree, "OntologyTree")) {
    stop("'tree' must be an OntologyTree object")
  }
  if(length(terms) == 0) {
    stop("No terms provided")
  }
  
  # Validate method-specific parameters
  if(method == "term") {
    submethod <- match.arg(submethod)
  }
  
  # Validate normalization type for Resnik
  if(method == "resnik") {
    normMap <- c("unif" = "Nunif", "max" = "Nmax", 
                 "sum" = "Nunivers", "none" = "none")
    if(!normType %in% names(normMap)) {
      stop("Invalid normalization type for Resnik similarity")
    }
  }
  if(is.null(ic_method)){
    ic_method = defaultIC(tree)
  }
 # ic_method <- match.arg(ic_method)
  # Validate weights for Wang and GOGO methods
  if(method %in% c("wang", "gogo") && !is.null(weights)) {
    if(is.null(names(weights))) {
      stop("weights must be a named vector")
    }
    if(any(weights >= 1)) {
      stop("All weights must be less than 1")
    }
    if(!all(names(weights) %in% c("is_a", "part_of"))) {
      stop("Invalid relation types in weights")
    }
  }
  terms <- term2id(tree,terms)
  # Calculate similarity based on method
  tryCatch({
    sim <- switch(method,
                  "term" = calculateTermSimilarity(tree, termIds = terms, method = submethod, annotUniverse=annotUniverse),
                  "ic" = calculateUniversalSim(tree, terms, verbose),
                  "resnik" = calculateResnikSim(tree, terms, normMap[normType], verbose),
                  "faith" = calculateFaithSim(tree, terms, method= ic_method, verbose = verbose),
                  "lin" = calculateLinSim(tree, terms, method= ic_method, verbose = verbose),
                  "rel" = calculateRelevanceSim(tree, terms, method= ic_method, verbose = verbose),
                  "simic" = calculateSimICSim(tree, terms,  method= ic_method,verbose = verbose),
                  "gogo" = {
                    if(is.null(weights)) {
                      weights <- c("is_a" = 0.4, "part_of" = 0.3)
                    }
                    calculateGOGOSim(tree, terms, weights, 
                                     method = if(useIgraph) "igraph" else "basic", 
                                     verbose = verbose)
                  },
                  "wang" = {
                    if(is.null(weights)) {
                      weights <- c("is_a" = 0.8, "part_of" = 0.6)
                    }
                    calculateWangSim(tree, terms, weights, 
                                     method = if(useIgraph) "igraph" else "basic", 
                                     verbose = verbose)
                  },
                  "anc" = calculateAncestorSim(tree, terms, verbose)
    )
    
    # Validate results
    if(!is.matrix(sim)) {
      stop("Method did not return a valid similarity matrix")
    }
    if(any(is.na(sim))) {
      warning("Similarity calculation produced NA values")
    }
    
    return(sim)
    
  }, error = function(e) {
    stop("Error calculating similarity: ", e$message)
  })
}

#' Calculate similarity between ontology terms
#' 
#' @param tree OntologyTree object
#' @param termIds Vector of term IDs to calculate similarity for
#' @param annotUniverse Optional subset of annotations to consider
#' @param method Similarity method ('kappa', 'jaccard', 'dice', or 'overlap')
#' @return Matrix of similarity scores
#' @importFrom methods as
#' @author Kai Guo
#' @export
calculateTermSimilarity <- function(
    tree, 
    termIds, 
    annotUniverse = NULL, 
    method = c("kappa", "jaccard", "dice", "overlap")
) {
  # Validate inputs
  if(length(tree@annotations$list) == 0) {
    stop("Annotations must be set in the ontology tree")
  }
  if(!requireNamespace("proxyC", quietly = TRUE)) {
    stop("Package 'proxyC' is required")
  }
  
  # Process annotation universe
  annotIndices <- if(!is.null(annotUniverse)) {
    which(tree@annotations$names %in% annotUniverse)
  } else {
    seq_along(tree@annotations$names)
  }
  validTerms <- validateAnnotatedTerms(tree, termIds)
  termIds <- termIds[validTerms]
  # Get annotation matrix for selected terms
  annotMatrix <- getTermAnnotations(tree, termIds)
  annotMatrix <- annotMatrix[, annotIndices, drop = FALSE]
  annotMatrix <- as(annotMatrix, "sparseMatrix")
  
  # Calculate similarity based on method
  method <- match.arg(method)
  simMatrix <- switch(method,
                      "kappa" = calculateKappaSimilarity(annotMatrix),
                      "overlap" = calculateOverlapSimilarity(annotMatrix),
                      proxyC::simil(annotMatrix, method = method)
  )
  
  # Format result matrix
  simMatrix <- as.matrix(simMatrix)
  diag(simMatrix) <- 1
  rownames(simMatrix) <- colnames(simMatrix) <- tree@termNames[termIds]
  
  return(simMatrix)
}

#' Calculate Kappa similarity between terms
#' 
#' @param annotMatrix Annotation matrix
#' @return Similarity matrix
#' @author Kai Guo
#' @keywords internal
calculateKappaSimilarity <- function(annotMatrix) {
  # Calculate total annotations
  totalAnnots <- ncol(annotMatrix)
  
  # Calculate observed agreement
  observed <- proxyC::simil(annotMatrix, method = "simple matching")
  
  # Calculate expected agreement
  annotCounts <- Matrix::rowSums(annotMatrix)
  nonAnnotCounts <- abs(Matrix::rowSums(annotMatrix - 1))
  expected <- (cross_multiply(annotCounts) + 
                 cross_multiply(nonAnnotCounts)) / 
    (totalAnnots * totalAnnots)
  
  # Calculate kappa
  kappa <- (observed - expected) / (1 - expected)
  return(kappa)
}

#' Calculate overlap similarity between terms
#' 
#' @param annotMatrix Annotation matrix
#' @return Similarity matrix
#' @author Kai Guo
#' @keywords internal
calculateOverlapSimilarity <- function(annotMatrix) {
  # Calculate annotation counts
  annotCounts <- Matrix::rowSums(annotMatrix)
  
  # Calculate Dice similarity and adjust by annotation counts
  diceSim <- proxyC::simil(annotMatrix, method = "dice")
  overlapSim <- diceSim * cross_sum(annotCounts) / (2 * cross_min(annotCounts))
  
  return(overlapSim)
}

#' Cross multiplication helper
#' 
#' @param x Numeric vector
#' @return Cross product matrix
#' @author Kai Guo
#' @keywords internal
cross_multiply <- function(x) {
  outer(x, x)
}

#' Cross sum helper
#' 
#' @param x Numeric vector
#' @return Matrix of pairwise sums
#' @author Kai Guo
#' @keywords internal
cross_sum <- function(x) {
  outer(x, x, "+")
}

#' Validate and filter terms with zero annotations
#' 
#' @param tree An OntologyTree object containing the ontology structure
#' @param termIds A numeric vector of term indices to validate
#' @return A logical vector indicating which terms have annotations
#' @author Kai Guo
#' @keywords internal
validateAnnotatedTerms <- function(tree, termIds) {
  # Get annotation counts for terms
  annotCounts <- getAnnotations(tree)[termIds]
  
  # Find terms with zero annotations
  zeroAnnots <- annotCounts == 0
  
  # Report removed terms
  if(any(zeroAnnots)) {
    message("Removing ", sum(zeroAnnots), " terms with no annotations")
    
    if(all(zeroAnnots)) {
      stop("No annotated terms remain")
    }
  }
  
  return(!zeroAnnots)
}

#' Calculate Universal similarity between terms
#' 
#' @param tree An OntologyTree object containing the ontology structure
#' @param terms A character vector of term IDs or names, or a numeric vector of term indices
#' @param verbose Logical scalar indicating whether to print progress messages
#' @return A numeric matrix containing similarity scores between term pairs
#' @examples
#' \dontrun{
#' sim <- calculateUniversalSim(tree, c("GO:0008150", "GO:0003674"))
#' }
#' @author Kai Guo
#' @export
calculateUniversalSim <- function(tree, terms, verbose = TRUE) {
  if(verbose) {
    message("Similarity method: Universal")
  }
  
  termIds <- term2id(tree, terms, strict = FALSE)
  IC <- calculateIC(tree, "universal", verbose = FALSE)[termIds]
  micaValues <- findMICA(tree, termIds, "universal", verbose = verbose)
  
  sim <- 2 * micaValues / cross_max(IC)
  dimnames(sim) <- list(tree@termNames[termIds], tree@termNames[termIds])
  
  return(sim)
}

#' Calculate Resnik similarity between terms
#' 
#' @param tree An OntologyTree object containing the ontology structure
#' @param terms A character vector of term IDs or names, or a numeric vector of term indices
#' @param normMethod Character scalar specifying the normalization method: "Nunif", "Nmax", "Nunivers", or "none"
#' @param verbose Logical scalar indicating whether to print progress messages
#' @details The normalization methods are:
#' * Nunif: Normalize by log of maximum annotation count
#' * Nmax: Normalize by maximum IC value
#' * Nunivers: Normalize by maximum pairwise IC value
#' * none: No normalization
#' @return A numeric matrix containing similarity scores between term pairs
#' @examples
#' \dontrun{
#' sim <- calculateResnikSim(tree, c("GO:0008150", "GO:0003674"), normMethod = "Nmax")
#' }
#' @author Kai Guo
#' @export
calculateResnikSim <- function(tree, terms, normMethod = "Nmax", verbose = TRUE) {
  if(!normMethod %in% c("Nunif", "Nmax", "Nunivers", "none")) {
    stop("normMethod must be one of: 'Nunif', 'Nmax', 'Nunivers', 'none'")
  }
  
  if(verbose) {
    message("Similarity method: Resnik (1999) + ", normMethod)
  }
  
  termIds <- term2id(tree, terms, strict = FALSE)
  IC <- calculateIC(tree, "annotation", verbose = FALSE)[termIds]
  
  validTerms <- validateAnnotatedTerms(tree, termIds)
  termIds <- termIds[validTerms]
  IC <- IC[validTerms]
  
  micaValues <- findMICA(tree, termIds, "annotation", verbose = verbose)
  
  sim <- switch(normMethod,
                "Nunif" = {
                  maxN <- attr(getAnnotations(tree), "N")
                  micaValues/log(maxN)
                },
                "Nmax" = micaValues/max(IC),
                "Nunivers" = micaValues/cross_max(IC),
                "none" = micaValues
  )
  
  dimnames(sim) <- list(tree@termNames[termIds], tree@termNames[termIds])
  sim[is.na(sim)] <- 1
  
  return(sim)
}

#' Calculate Lin similarity between terms
#' 
#' @param tree An OntologyTree object containing the ontology structure
#' @param terms A character vector of term IDs or names, or a numeric vector of term indices
#' @param method Character scalar specifying the IC calculation method, e.g., "annotation"
#' @param verbose Logical scalar indicating whether to print progress messages
#' @return A numeric matrix containing similarity scores between term pairs
#' @examples
#' \dontrun{
#' sim <- calculateLinSim(tree, c("GO:0008150", "GO:0003674"))
#' }
#' @author Kai Guo
#' @export
calculateLinSim <- function(tree, terms, method = "annotation", verbose = TRUE) {
  if(verbose) {
    message("Similarity method: Lin (1998)")
  }
  
  termIds <- term2id(tree, terms, strict = FALSE)
  IC <- calculateIC(tree, method, verbose = FALSE)[termIds]
  if(method == "annotation"){
    validTerms <- validateAnnotatedTerms(tree, termIds)
    termIds <- termIds[validTerms]
    IC <- IC[validTerms]
  }
  micaValues <- findMICA(tree, termIds, method, verbose = verbose)
  
  sim <- 2*micaValues/cross_sum(IC) 
  dimnames(sim) <- list(tree@termNames[termIds], tree@termNames[termIds])
  sim[is.na(sim)] <- 1
  
  return(sim)
}

#' Calculate FaITH similarity between terms
#' 
#' @param tree An OntologyTree object containing the ontology structure
#' @param terms A character vector of term IDs or names, or a numeric vector of term indices
#' @param method Character scalar specifying the IC calculation method, e.g., "annotation"
#' @param verbose Logical scalar indicating whether to print progress messages
#' @return A numeric matrix containing similarity scores between term pairs
#' @examples
#' \dontrun{
#' sim <- calculateFaithSim(tree, c("GO:0008150", "GO:0003674"))
#' }
#' @author Kai Guo
#' @export
calculateFaithSim <- function(tree, terms, method = "annotation", verbose = TRUE) {
  if(verbose) {
    message("Similarity method: FaITH (2010)")
  }
  
  termIds <- term2id(tree, terms, strict = FALSE)
  IC <- calculateIC(tree, method, verbose = FALSE)[termIds]
  
  validTerms <- validateAnnotatedTerms(tree, termIds)
  termIds <- termIds[validTerms]
  IC <- IC[validTerms]
  
  micaValues <- findMICA(tree, termIds, method, verbose = verbose)
  
  sim <- micaValues/(cross_sum(IC) - micaValues)
  dimnames(sim) <- list(tree@termNames[termIds], tree@termNames[termIds])
  sim[is.na(sim)] <- 1
  
  return(sim)
}

#' Calculate Relevance similarity between terms
#' 
#' @param tree An OntologyTree object containing the ontology structure
#' @param terms A character vector of term IDs or names, or a numeric vector of term indices
#' @param method Character scalar specifying the IC calculation method, e.g., "annotation"
#' @param verbose Logical scalar indicating whether to print progress messages
#' @return A numeric matrix containing similarity scores between term pairs
#' @details The similarity score is modified by a correction factor of (1 - e^(-IC_MICA))
#' @examples
#' \dontrun{
#' sim <- calculateRelevanceSim(tree, c("GO:0008150", "GO:0003674"))
#' }
#' @author Kai Guo
#' @export
calculateRelevanceSim <- function(tree, terms, method = "annotation", verbose = TRUE) {
  if(verbose) {
    message("Similarity method: Relevance (2006)")
  }
  
  termIds <- term2id(tree, terms, strict = FALSE)
  IC <- calculateIC(tree, method, verbose = FALSE)[termIds]
  
  validTerms <- validateAnnotatedTerms(tree, termIds)
  termIds <- termIds[validTerms]
  IC <- IC[validTerms]
  
  micaValues <- findMICA(tree, termIds, method, verbose = verbose)
  
  sim <- 2 * micaValues/cross_sum(IC)
  sim[is.na(sim)] <- 1
  dimnames(sim) <- list(tree@termNames[termIds], tree@termNames[termIds])
  
  correctionFactor <- 1 - exp(-micaValues)
  return(correctionFactor * sim)
}

#' Calculate SimIC similarity between terms
#' 
#' @param tree An OntologyTree object containing the ontology structure
#' @param terms A character vector of term IDs or names, or a numeric vector of term indices
#' @param method Character scalar specifying the IC calculation method, e.g., "annotation"
#' @param verbose Logical scalar indicating whether to print progress messages
#' @return A numeric matrix containing similarity scores between term pairs
#' @details The similarity score is modified by a correction factor of (1 - 1/(1 + IC_MICA))
#' @examples
#' \dontrun{
#' sim <- calculateSimICSim(tree, c("GO:0008150", "GO:0003674"))
#' }
#' @author Kai Guo
#' @export
calculateSimICSim <- function(tree, terms, method = "annotation", verbose = TRUE) {
  if(verbose) {
    message("Similarity method: SimIC (2010)")
  }
  
  termIds <- term2id(tree, terms, strict = FALSE)
  IC <- calculateIC(tree, method, verbose = FALSE)[termIds]
  
  validTerms <- validateAnnotatedTerms(tree, termIds)
  termIds <- termIds[validTerms]
  IC <- IC[validTerms]
  
  micaValues <- findMICA(tree, termIds, method, verbose = verbose)
  
  sim <- 2 * micaValues/cross_sum(IC)
  sim[is.na(sim)] <- 1
  dimnames(sim) <- list(tree@termNames[termIds], tree@termNames[termIds])
  
  correctionFactor <- 1 - 1/(1 + micaValues)
  return(correctionFactor * sim)
}

#' Calculate GOGO similarity between terms
#' 
#' @param tree An OntologyTree object containing the ontology structure
#' @param terms A character vector of term IDs or names, or a numeric vector of term indices
#' @param contributionFactors Named numeric vector specifying contribution factors for each relation type
#' @param method Character scalar specifying calculation method: "igraph" or "basic"
#' @param verbose Logical scalar indicating whether to print progress messages
#' @return A numeric matrix containing similarity scores between term pairs
#' @details The contributionFactors must be named and correspond to relation types in the tree.
#'    All values must be less than 1.
#' @examples
#' \dontrun{
#' factors <- c("is_a" = 0.4, "part_of" = 0.3)
#' sim <- calculateGOGOSim(tree, c("GO:0008150", "GO:0003674"), 
#'                        contributionFactors = factors)
#' }
#' @author Kai Guo
#' @export
calculateGOGOSim <- function(tree, terms, 
                             contributionFactors = c("is_a" = 0.4, "part_of" = 0.3),
                             method = "igraph", 
                             verbose = TRUE) {
  if(verbose) {
    message("Similarity method: GOGO (2018)")
  }
  
  if(length(tree@edgeTypes) == 0) {
    stop("Relations not set in the tree")
  }
  
  relationLevels <- attr(tree@edgeTypes, "levels")
  if(is.null(names(contributionFactors))) {
    stop("contributionFactors must be a named numeric vector")
  }
  
  # Process relation names directly
  names(contributionFactors) <- gsub("isa", "is_a", names(contributionFactors))
  names(contributionFactors) <- gsub(" ", "_", names(contributionFactors))
  
  # Extend contribution factors if needed
  if(!is.null(tree@relTree)) {
    cf <- contributionFactors
    for(nm in names(contributionFactors)) {
      if(nm %in% tree@relTree@termNames) {
        offspring <- getOffspringCount(tree@relTree, nm)
        if(length(offspring)) {
          cf[offspring] <- contributionFactors[nm]
        }
      }
    }
    contributionFactors <- cf
  }
  
  if(length(setdiff(relationLevels, names(contributionFactors))) > 0) {
    stop("Contribution factors must be provided for all relations")
  }
  if(any(contributionFactors >= 1)) {
    stop("All contribution factors must be less than 1")
  }
  
  termIds <- term2id(tree, terms, strict = FALSE)
  sim <- matrix(0, nrow = length(termIds), ncol = length(termIds))
  dimnames(sim) <- list(tree@termNames[termIds], tree@termNames[termIds])
  diag(sim) <- 1
  
  if(length(termIds) <= 1) return(sim)
  
  if(method == "igraph") {
    g <- convertToIgraph(tree)
    em <- get.edgelist(g)
    c <- max(contributionFactors)/(1 - max(contributionFactors))
    nc <- getChildCount(tree)
    E(g)$weight <- contributionFactors[E(g)$relation] + 1/(c + nc[em[,1]])
    
    ancestors <- getGroupAncestors(tree, termIds, includeSelf = TRUE)
    d <- distances(g, v = ancestors, to = termIds, mode = "out", 
                   weights = -log(E(g)$weight))
    s <- exp(-d)
    sim <- wang_similarity_sv(s)
  } else {
    sim <- wang_similarity(tree, termIds, 
                           unname(contributionFactors[relationLevels]), TRUE)
  }
  
  dimnames(sim) <- list(tree@termNames[termIds], tree@termNames[termIds])
  return(sim)
}

#' Calculate Wang similarity between terms
#' 
#' @param tree An OntologyTree object containing the ontology structure
#' @param terms A character vector of term IDs or names, or a numeric vector of term indices
#' @param contributionFactors Named numeric vector specifying contribution factors for each relation type
#' @param method Character scalar specifying calculation method: "igraph" or "basic"
#' @param verbose Logical scalar indicating whether to print progress messages
#' @return A numeric matrix containing similarity scores between term pairs
#' @details The contributionFactors must be named and correspond to relation types in the tree.
#'    All values must be less than 1.
#' @examples
#' \dontrun{
#' # Calculate similarity with default contribution factors
#' sim <- calculateWangSim(tree, c("GO:0008150", "GO:0003674"))
#' 
#' # Use custom contribution factors
#' factors <- c("is_a" = 0.8, "part_of" = 0.6)
#' sim <- calculateWangSim(tree, terms, contributionFactors = factors)
#' }
#' @author Kai Guo
#' @export
calculateWangSim <- function(tree, terms, 
                             contributionFactors = c("is_a" = 0.8, "part_of" = 0.6),
                             method = "igraph", 
                             verbose = TRUE) {
  if(verbose) {
    message("Similarity method: Wang (2007)")
  }
  
  if(length(tree@edgeTypes) == 0) {
    stop("Relations not set in the tree")
  }
  
  relationLevels <- attr(tree@edgeTypes, "levels")
  if(is.null(names(contributionFactors))) {
    stop("contributionFactors must be a named numeric vector")
  }
  
  # Process relation names directly
  names(contributionFactors) <- gsub("isa", "is_a", names(contributionFactors))
  names(contributionFactors) <- gsub(" ", "_", names(contributionFactors))
  
  # Extend contribution factors if needed
  if(!is.null(tree@relTree)) {
    cf <- contributionFactors
    for(nm in names(contributionFactors)) {
      if(nm %in% tree@relTree@termNames) {
        offspring <- getOffspringCount(tree@relTree, nm)
        if(length(offspring)) {
          cf[offspring] <- contributionFactors[nm]
        }
      }
    }
    contributionFactors <- cf
  }
  
  if(length(setdiff(relationLevels, names(contributionFactors))) > 0) {
    stop("Contribution factors must be provided for all relations")
  }
  if(any(contributionFactors >= 1)) {
    stop("All contribution factors must be less than 1")
  }
  
  termIds <- term2id(tree, terms, strict = FALSE)
  sim <- matrix(0, nrow = length(termIds), ncol = length(termIds))
  dimnames(sim) <- list(tree@termNames[termIds], tree@termNames[termIds])
  diag(sim) <- 1
  
  if(length(termIds) <= 1) return(sim)
  
  if(method == "igraph") {
    g <- convertToIgraph(tree)
    E(g)$weight <- contributionFactors[E(g)$relation]
    
    ancestors <- getGroupAncestors(tree, termIds, includeSelf = TRUE)
    d <- distances(g, v = ancestors, to = termIds, mode = "out", 
                   weights = -log(E(g)$weight))
    s <- exp(-d)
    sim <- wang_similarity_sv(s)
  } else {
    sim <- wang_similarity(tree, termIds, 
                           unname(contributionFactors[relationLevels]))
  }
  
  dimnames(sim) <- list(tree@termNames[termIds], tree@termNames[termIds])
  return(sim)
}

#' Calculate ancestor-based similarity between terms
#' 
#' @description 
#' Calculates semantic similarity between terms based on their common ancestors.
#' 
#' @details 
#' For any two terms a and b, let Sa and Sb be their respective sets of ancestor terms 
#' (including the terms themselves). The similarity is calculated as:
#' 
#' similarity = |Sa ∩ Sb| / |Sa ∪ Sb|
#' 
#' where |X| denotes the size of set X.
#' 
#' @param tree An OntologyTree object
#' @param terms A character vector of term IDs or names, or a numeric vector of term indices
#' @param verbose Logical scalar indicating whether to print progress messages
#' @return A numeric matrix containing similarity scores between term pairs
#' @examples
#' \dontrun{
#' # Calculate ancestor-based similarity for specific terms
#' sim <- calculateAncestorSim(tree, c("GO:0008150", "GO:0003674"))
#' 
#' # Calculate similarity using term indices
#' sim <- calculateAncestorSim(tree, 1:10)
#' }
#' @author Kai Guo
#' @export
calculateAncestorSim <- function(tree, terms, verbose = TRUE) {
  # Convert terms to IDs
  termIds <- term2id(tree, terms, strict = FALSE)
  
  # Calculate similarity matrix
  similarity <- calculateAncestorSimilarity(tree, termIds)
  
  # Add term names
  dimnames(similarity) <- list(tree@termNames[termIds], tree@termNames[termIds])
  
  return(similarity)
}
