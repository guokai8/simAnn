#' @title Ontology Tree Class Definition
#' @description A class representing hierarchical relationships in ontology
#' @slot termNames Character vector of term names
#' @slot termCount Integer scalar of total terms
#' @slot edgeCount Integer scalar of total relations
#' @slot parentMap List of parent term indices
#' @slot childMap List of child term indices
#' @slot edgeTypes List of semantic relations between terms
#' @slot relTree Ontology tree object for relation types
#' @slot source Character scalar marking data source
#' @slot rootNode Integer scalar of root term index
#' @slot leafNodes Integer vector of leaf term indices
#' @slot annotations List containing term annotations
#' @slot n_parents Vector containing the parents number
#' @slot n_children Vector containing the children number
#' @slot metadata Additional term metadata
#' @slot termStats List for storing arbitrary term-level statistics
#' @author Kai Guo
#' @export
OntologyTree <- setClass("OntologyTree",
                         slots = c(
                           termNames = "character",
                           termCount = "integer", 
                           edgeCount = "integer",
                           parentMap = "list",
                           childMap = "list", 
                           edgeTypes = "list",
                           relTree = "ANY",
                           source = "character",
                           rootNode = "integer", 
                           leafNodes = "integer",
                           annotations = "list",
                           metadata = "ANY",
                           n_parents = "ANY",
                           n_children = "ANY",
                           termStats = "ANY"
                         )
)

#' Create an ontology tree
#' 
#' @param parentTerms Character vector of parent terms
#' @param childTerms Character vector of child terms 
#' @param edgeTypes Character vector of relation types
#' @param relTree Ontology tree object for relation hierarchy
#' @param source Data source identifier
#' @param annotations Term annotations list
#' @param removeCycles Remove cyclic paths if TRUE
#' @param removeIsolated Remove isolated rings
#' @param verbose Print progress messages
#' @return OntologyTree object
#' @author Kai Guo
#' @export
buildOntologyTree <- function(
    parentTerms,
    childTerms = NULL,
    edgeTypes = NULL,
    relTree = NULL,
    source = "Ontology",
    annotations = NULL,
    removeCycles = FALSE,
    removeIsolated = FALSE,
    verbose = TRUE
) {
  # Handle edge list input format
  if(is.null(childTerms)) {
    parsed <- parseEdgeList(parentTerms)
    parentTerms <- parsed$parents
    childTerms <- parsed$children
  }
  
  # Validate inputs
  validateInputs(parentTerms, childTerms)
  
  # Remove self-loops
  selfLoops = parentTerms == childTerms
  if(any(selfLoops)) {
    if(verbose) {
      message("removed ", sum(selfLoops), " relations where parents are the same as children.")
    }
    parentTerms = parentTerms[!selfLoops]
    childTerms = childTerms[!selfLoops]
  }
  
  if(length(parentTerms) == 0 || length(childTerms) == 0) {
    stop("There is no relation.")
  }
  
  hasRelations = !is.null(edgeTypes)
  edgeTypesList = list()
  
  if(hasRelations && any(selfLoops)) {
    edgeTypes = edgeTypes[!selfLoops]
  }
  
  # Build term mappings  
  terms = sort(unique(c(parentTerms, childTerms)))
  termCount = length(terms)
  termIndices = structure(seq_along(terms), names = terms)
  
  # Build maps
  parentMap = buildParentMap(parentTerms, childTerms, termIndices, termCount)
  childMap = buildChildMap(parentTerms, childTerms, termIndices, termCount)
  edgeTypesList = if(hasRelations) buildRelationMap(edgeTypes, parentTerms, termIndices, termCount) else list()
  
  # Find root and leaf nodes
  root = which(vapply(parentMap, length, FUN.VALUE = integer(1)) == 0)
  leaves = which(vapply(childMap, length, FUN.VALUE = integer(1)) == 0)
  
  if(length(root) == 0) {
    stop("Cannot find the root. There might exist a cycle.")
  }
  
  # Handle multiple roots if needed
  if(length(root) > 1) {
    result = addSuperRoot(root, parentMap, childMap, edgeTypesList, terms, verbose)
    parentMap = result$parentMap
    childMap = result$childMap
    edgeTypesList = result$relationMap
    terms = result$terms
    root = result$root
    termCount = length(terms)
  }
  
  if(!is.null(relTree)) {
    if(!inherits(relTree, "OntologyTree")) {
      stop("`relTree` should be constructed by `buildOntologyTree()`.")
    }
  }
  
  # Calculate node statistics
  n_parents = vapply(parentMap, length, FUN.VALUE = integer(1))
  n_children = vapply(childMap, length, FUN.VALUE = integer(1))
  
  # Initialize termStats
  termStats = list()
  
  # Create OntologyTree object
  tree = OntologyTree(
    termNames = terms,
    termCount = termCount,
    edgeCount = sum(vapply(childMap, length, FUN.VALUE = integer(1))),
    parentMap = parentMap,
    childMap = childMap,
    edgeTypes = edgeTypesList,
    relTree = relTree,
    source = source,
    rootNode = root,
    leafNodes = leaves,
    n_parents = n_parents,
    n_children = n_children,
    termStats = termStats
  )
  
  tree@annotations = list(list = vector("list", 0), names = character(0))
  tree = addOntologyTerms(tree, annotations)
  return(tree)
}

#' Add annotations to the ontology tree
#' 
#' @param tree An OntologyTree object
#' @param annotations Term annotation list
#' @return Updated OntologyTree object
#' @author Kai Guo
#' @export
addOntologyTerms = function(tree, annotations) {
  if(!is.null(annotations)) {
    annotations = lapply(annotations, as.character)
    allAnnotations = unique(unlist(annotations))
    annotationCount = length(allAnnotations)
    annotationIndices = structure(seq_along(allAnnotations), names = allAnnotations)
    annotations = lapply(annotations, function(x) unname(annotationIndices[as.character(x)]))
    
    termCount = tree@termCount
    annotationMap = rep(list(integer(0)), termCount)
    termIndices = structure(seq_along(tree@termNames), names = tree@termNames)
    
    commonTerms = intersect(tree@termNames, names(annotations))
    if(length(commonTerms) == 0) {
      stop("No overlap between annotation names and terms.")
    }
    annotationMap[termIndices[commonTerms]] = annotations[commonTerms]
    
    tree@annotations = list(
      list = annotationMap,
      names = allAnnotations
    )
  } else {
    tree@annotations = list(list = vector("list", 0), names = character(0))
  }
  return(tree)
}

#' Create ontology tree from igraph object
#' 
#' @param graph An igraph object
#' @param edgeTypes Optional vector of edge relation types
#' @param verbose Print progress messages
#' @return OntologyTree object
#' @importFrom igraph get.edgelist is.directed is.dag
#' @author Kai Guo
#' @export
fromGraphToOntology <- function(graph, edgeTypes = NULL, verbose = TRUE) {
  if(!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for this function")
  }
  
  if(!igraph::is.directed(graph)) {
    stop("Input graph must be directed")
  }
  
  if(!igraph::is.dag(graph)) {
    stop("Input graph must be a directed acyclic graph (DAG)")
  }
  
  edges = igraph::get.edgelist(graph)
  buildOntologyTree(
    parentTerms = as.character(edges[, 1]),
    childTerms = as.character(edges[, 2]),
    edgeTypes = edgeTypes,
    source = "igraph object",
    verbose = verbose
  )
}

#' Add term statistics to OntologyTree
#' 
#' @param tree OntologyTree object
#' @param name Name of the statistic
#' @param data Data to store (must have same length as number of terms)
#' @return Updated OntologyTree object
#' @author Kai Guo
#' @export
addTermStats <- function(tree, name, data) {
  # Input validation
  if(!inherits(tree, "OntologyTree")) {
    stop("First argument must be an OntologyTree object")
  }
  if(!is.character(name) || length(name) != 1) {
    stop("'name' must be a single character string")
  }
  if(length(data) != tree@termCount) {
    stop("Length of data must match the number of terms in the tree")
  }
  
  # Add the named data to termStats
  tree@termStats[[name]] <- data
  return(tree)
}

#' Get term statistics from OntologyTree
#' 
#' @param tree OntologyTree object
#' @param name Name of the statistic to retrieve
#' @return Requested statistic data
#' @author Kai Guo
#' @export
getTermStats <- function(tree, name) {
  if(!name %in% names(tree@termStats)) {
    stop(sprintf("Statistic '%s' not found in termStats", name))
  }
  return(tree@termStats[[name]])
}

#' List available term statistics
#' 
#' @param tree OntologyTree object
#' @return Character vector of available statistic names
#' @author Kai Guo
#' @export
listTermStats <- function(tree) {
  return(names(tree@termStats))
}

#' Remove term statistics from OntologyTree
#' 
#' @param tree OntologyTree object
#' @param name Name of the statistic to remove
#' @return Updated OntologyTree object
#' @author Kai Guo
#' @export
removeTermStats <- function(tree, name) {
  if(!name %in% names(tree@termStats)) {
    warning(sprintf("Statistic '%s' not found in termStats", name))
    return(tree)
  }
  tree@termStats[[name]] <- NULL
  return(tree)
}

#' Get number of parents for each term
#' 
#' @param tree OntologyTree object
#' @param terms Optional subset of terms to get counts for (character or numeric)
#' @return Integer vector of parent counts
#' @examples
#' \dontrun{
#' # Get parent counts for all terms
#' counts <- getParentCount(tree)
#' 
#' # Get counts for specific terms
#' counts <- getParentCount(tree, c("GO:0008150", "GO:0003674"))
#' 
#' # Get counts using numeric indices
#' counts <- getParentCount(tree, 1:10)
#' }
#' @author Kai Guo
#' @export
getParentCount <- function(tree, terms = NULL) {
  counts <- tree@n_parents
  names(counts) <- tree@termNames
  
  if(!is.null(terms)) {
    if(is.numeric(terms)) {
      return(counts[terms])
    } else {
      termIds <- term2id(tree, terms, addNames = TRUE)
      return(counts[termIds[terms]])
    }
  }
  
  return(counts)
}

#' Get number of children for each term
#' 
#' @param tree OntologyTree object
#' @param terms Optional subset of terms to get counts for (character or numeric)
#' @return Integer vector of child counts
#' @examples
#' \dontrun{
#' # Get child counts for all terms
#' counts <- getChildCount(tree)
#' 
#' # Get counts for specific terms
#' counts <- getChildCount(tree, c("GO:0008150", "GO:0003674"))
#' 
#' # Get counts using numeric indices
#' counts <- getChildCount(tree, 1:10)
#' }
#' @author Kai Guo
#' @export
getChildCount <- function(tree, terms = NULL) {
  counts <- tree@n_children
  names(counts) <- tree@termNames
  
  if(!is.null(terms)) {
    if(is.numeric(terms)) {
      return(counts[terms])
    } else {
      termIds <- term2id(tree, terms, addNames = TRUE)
      return(counts[termIds[terms]])
    }
  }
  
  return(counts)
}

#' Get terms in an ontology tree
#' 
#' @param tree An OntologyTree object
#' @return Various information about terms in the tree:
#'   * `terms()`: character vector of all term names
#'   * `nTerms()`: integer scalar of total term count
#'   * `nEdges()`: integer scalar of total edge count
#'   * `nLeaves()`: integer scalar of total leaf term count
#' @examples
#' \dontrun{
#' # Create a simple tree
#' parents <- c("a", "a", "b", "b", "c", "d")
#' children <- c("b", "c", "c", "d", "e", "f")
#' tree <- buildOntologyTree(parents, children)
#' 
#' # Get term information
#' terms(tree)
#' nTerms(tree)
#' }
#' 
#' @author Kai Guo
#' @export
terms <- function(tree) {
  tree@termNames
}

#' @rdname terms
#' @author Kai Guo
#' @export
nTerms <- function(tree) {
  tree@termCount
}

#' @rdname terms
#' @author Kai Guo
#' @export
nEdges <- function(tree) {
  sum(vapply(tree@childMap, length, FUN.VALUE = integer(1)))
}

#' @rdname terms
#' @author Kai Guo
#' @export
nLeaves <- function(tree) {
  length(tree@leafNodes)
}

#' Get root and leaf terms
#' 
#' @param tree An OntologyTree object
#' @param asNames Whether to return term names instead of indices
#' @return Term information:
#'   * `root()`: root term (name or index)
#'   * `leaves()`: leaf terms (names or indices)
#'   * `isLeaf()`: logical vector indicating which terms are leaves
#' @examples
#' \dontrun{
#' # Get root and leaves
#' root(tree)
#' leaves(tree)
#' 
#' # Get indices instead of names
#' root(tree, asNames = FALSE)
#' leaves(tree, asNames = FALSE)
#' 
#' # Check if terms are leaves
#' isLeaf(tree, c("GO:0008150", "GO:0003674"))
#' }
#' @author Kai Guo
#' @export
root <- function(tree, asNames = TRUE) {
  if(asNames) {
    tree@termNames[tree@rootNode]
  } else {
    tree@rootNode
  }
}

#' @rdname root
#' @author Kai Guo
#' @export
leaves <- function(tree, asNames = TRUE) {
  if(asNames) {
    tree@termNames[tree@leafNodes]
  } else {
    tree@leafNodes
  }
}

#' @param terms Character vector of term names to check
#' @rdname root
#' @author Kai Guo
#' @export
isLeaf <- function(tree, terms) {
  termIds <- term2id(tree, terms)
  isLeaf <- termIds %in% tree@leafNodes
  stats::setNames(isLeaf, tree@termNames[termIds])
}

#' Check if terms exist in the ontology tree
#' 
#' @param tree An OntologyTree object
#' @param terms A character vector of term names or IDs to check
#' @return A logical vector indicating which terms exist in the tree
#' @examples
#' \dontrun{
#' # Check if single term exists
#' exists <- hasTerms(tree, "GO:0008150")
#' 
#' # Check multiple terms
#' exists <- hasTerms(tree, c("GO:0008150", "GO:0003674"))
#' }
#' @author Kai Guo
#' @export
hasTerms <- function(tree, terms) {
  terms %in% tree@termNames
}

#' @rdname is_annotations
#' @returns `is_annotation()` returns a logical scalar.
is_annotation <- function(tree) {
  return(length(tree@annotations$list) > 0)
}

#' provide default IC method
defaultIC <- function(tree) {
  if (is_annotation(tree)) {
    return("annotation")
  } else {
    return("offspring")
  }
}

#' Get parents for a specific term
#' 
#' @param tree OntologyTree object
#' @param termId Term ID to get parents for
#' @param asNames Whether to return term names instead of indices
#' @return Vector of parent terms
#' @author Kai Guo
#' @export
getParents <- function(tree, termId, asNames = TRUE) {
  if(!termId %in% seq_len(tree@termCount)) {
    stop("Invalid term ID")
  }
  parents <- tree@parentMap[[termId]]
  if(asNames) {
    return(tree@termNames[parents])
  }
  return(parents)
}

#' Get children for a specific term
#' 
#' @param tree OntologyTree object
#' @param termId Term ID to get children for
#' @param asNames Whether to return term names instead of indices
#' @return Vector of child terms
#' @author Kai Guo
#' @export
getChildren <- function(tree, termId, asNames = TRUE) {
  if(!termId %in% seq_len(tree@termCount)) {
    stop("Invalid term ID")
  }
  children <- tree@childMap[[termId]]
  if(asNames) {
    return(tree@termNames[children])
  }
  return(children)
}

# Helper Functions

#' Parse edge list into parent and child components
#' @param edges Character vector of edges
#' @return List of parent and child vectors
#' @author Kai Guo
#' @keywords internal
parseEdgeList <- function(edges) {
  separator <- if(any(grepl("\\s+-\\s+", edges))) "\\s+-\\s+" else "\\s*-\\s*"
  parts <- strsplit(edges, separator)
  list(
    parents = sapply(parts, `[`, 1),
    children = sapply(parts, `[`, 2)
  )
}

#' Validate parent and child inputs
#' @param parents Parent term vector
#' @param children Child term vector
#' @author Kai Guo
#' @keywords internal
validateInputs <- function(parents, children) {
  if(!is.character(parents)) stop("Parents must be character vector")
  if(!is.character(children)) stop("Children must be character vector")
  if(length(parents) != length(children)) {
    stop("Parent and child vectors must have equal length")
  }
}

#' Build parent term mapping
#' @param parents Parent terms
#' @param children Child terms
#' @param termIndices Term index mapping
#' @param termCount Total term count
#' @return Parent mapping list
#' @author Kai Guo
#' @keywords internal
buildParentMap <- function(parents, children, termIndices, termCount) {
  parentList <- split(unname(termIndices[parents]), children)
  result <- rep(list(integer(0)), termCount) 
  result[termIndices[names(parentList)]] <- parentList
  result
}

#' Build child term mapping 
#' @param parents Parent terms
#' @param children Child terms
#' @param termIndices Term index mapping
#' @param termCount Total term count
#' @return Child mapping list
#' @author Kai Guo
#' @keywords internal
buildChildMap <- function(parents, children, termIndices, termCount) {
  childList <- split(unname(termIndices[children]), parents)
  result <- rep(list(integer(0)), termCount)
  result[termIndices[names(childList)]] <- childList
  result
}

#' Build relation type mapping
#' @param relations Relation types
#' @param parents Parent terms 
#' @param termIndices Term index mapping
#' @param termCount Total term count
#' @return Relation mapping list
#' @author Kai Guo
#' @keywords internal
buildRelationMap <- function(relations, parents, termIndices, termCount) {
  if(is.null(relations)) return(NULL)
  
  # First ensure relations is a factor
  if(!is.factor(relations)) {
    relations <- as.factor(relations)
  }
  
  # Convert to integer (factor levels) and split
  relationInts <- as.integer(relations)
  relationList <- split(relationInts, parents)
  
  # Create result list
  result <- rep(list(integer(0)), termCount)
  result[termIndices[names(relationList)]] <- relationList
  attr(result, "levels") <- levels(relations)
  
  return(result)
}

#' Add super root node to handle multiple roots
#' @param roots Root term indices
#' @param parentMap Parent mapping
#' @param childMap Child mapping
#' @param relationMap Relation mapping
#' @param terms Term names
#' @param verbose Print messages flag
#' @return List of updated mappings
#' @author Kai Guo
#' @keywords internal
addSuperRoot <- function(roots, parentMap, childMap, relationMap, terms, verbose) {
  if(verbose) {
    reportMultipleRoots(roots, terms)
  }
  
  superRootIndex <- length(terms) + 1L
  parentMap[[superRootIndex]] <- integer(0)
  childMap[[superRootIndex]] <- roots
  
  for(root in roots) {
    parentMap[[root]] <- superRootIndex
  }
  
  if(!is.null(relationMap)) {
    relationMap[[superRootIndex]] <- rep(1L, length(roots))
  }
  
  terms <- c(terms, "SUPER_ROOT")
  
  list(
    parentMap = parentMap,
    childMap = childMap,
    relationMap = relationMap,
    terms = terms,
    root = superRootIndex
  )
}

#' Report multiple root nodes found
#' @param roots Root indices
#' @param terms Term names
#' @author Kai Guo
#' @keywords internal
reportMultipleRoots <- function(roots, terms) {
  rootTerms <- terms[roots]
  if(length(rootTerms) > 5) {
    sample <- paste(rootTerms[1:5], collapse = ", ")
    message(sprintf("Multiple roots found: %s and %d others", 
                    sample, length(rootTerms) - 5))
  } else {
    message(sprintf("Multiple roots found: %s", 
                    paste(rootTerms, collapse = ", ")))
  }
  message("Adding super root node")
}

#' Show method for OntologyTree
#' @param object An OntologyTree object
#' @export
setMethod("show", "OntologyTree", function(object) {
  cat("OntologyTree Object\n")
  cat("-------------------\n")
  cat("Source:       ", object@source, "\n")
  cat("Total Terms:  ", object@termCount, "\n")
  cat("Total Edges:  ", object@edgeCount, "\n")
  cat("Root Term:    ", object@termNames[object@rootNode], "\n")
  cat("Leaf Terms:   ", length(object@leafNodes), "\n")
})

#' Summary method for OntologyTree
#' @param object An OntologyTree object
#' @export
setMethod("summary", "OntologyTree", function(object) {
  cat("OntologyTree Summary\n")
  cat("---------------------\n")
  cat("Source:          ", object@source, "\n")
  cat("Total Terms:     ", object@termCount, "\n")
  cat("Total Relations: ", object@edgeCount, "\n")
  cat("Root Term:       ", object@termNames[object@rootNode], "\n")
  cat("Number of Leaves:", length(object@leafNodes), "\n")
  
  # Parent & Child statistics
  cat("Max Parents:     ", max(object@n_parents), "\n")
  cat("Max Children:    ", max(object@n_children), "\n")
  
  # Relationship types
  if (!is.null(object@edgeTypes) && length(object@edgeTypes) > 0) {
    cat("Unique Relations:", length(unique(unlist(object@edgeTypes))), "\n")
  } else {
    cat("Unique Relations: None defined\n")
  }
  
  # Annotation summary
  if (!is.null(object@annotations$list) && length(object@annotations$list) > 0) {
    cat("Annotated Terms: ", sum(sapply(object@annotations$list, length) > 0), "\n")
  } else {
    cat("Annotated Terms: 0\n")
  }
})

#' Print an OntologyTree object
#' 
#' @param x An OntologyTree object
#' @param ... Additional arguments (ignored)
#' @return None (prints summary to console)
#' @export
print.OntologyTree <- function(x, ...) {
  cat("OntologyTree Object\n")
  cat("-------------------\n")
  cat("Source:       ", x@source, "\n")
  cat("Total Terms:  ", x@termCount, "\n")
  cat("Total Edges:  ", x@edgeCount, "\n")
  cat("Root Term:    ", x@termNames[x@rootNode], "\n")
  cat("Leaf Terms:   ", length(x@leafNodes), "\n")
}

#' Summary of an OntologyTree object
#' 
#' @param object An OntologyTree object
#' @param ... Additional arguments (ignored)
#' @return A structured summary of the ontology tree
#' @export
summary.OntologyTree <- function(object, ...) {
  cat("OntologyTree Summary\n")
  cat("---------------------\n")
  cat("Source:          ", object@source, "\n")
  cat("Total Terms:     ", object@termCount, "\n")
  cat("Total Relations: ", object@edgeCount, "\n")
  cat("Root Term:       ", object@termNames[object@rootNode], "\n")
  cat("Number of Leaves:", length(object@leafNodes), "\n")
  
  # Parent & Child statistics
  cat("Max Parents:     ", max(object@n_parents), "\n")
  cat("Max Children:    ", max(object@n_children), "\n")
  
  # Relationship types
  if (!is.null(object@edgeTypes) && length(object@edgeTypes) > 0) {
    cat("Unique Relations:", length(unique(unlist(object@edgeTypes))), "\n")
  } else {
    cat("Unique Relations: None defined\n")
  }
  
  # Annotation summary
  if (!is.null(object@annotations$list) && length(object@annotations$list) > 0) {
    cat("Annotated Terms: ", sum(sapply(object@annotations$list, length) > 0), "\n")
  } else {
    cat("Annotated Terms: 0\n")
  }
}
