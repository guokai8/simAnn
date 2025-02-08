#' Get number of offspring terms for each term
#' 
#' @param tree OntologyTree object
#' @param terms Optional subset of terms to get counts for
#' @param useCache Whether to use cached counts
#' @param includeSelf Whether to include the term itself in counts
#' @return Named integer vector of offspring counts
#' @author Kai Guo
#' @export
getOffspringCount <- function(tree, terms = NULL, useCache = TRUE, includeSelf = FALSE) {
  # Try to get from cache first
  if(useCache && "offspringCount" %in% names(tree@termStats)) {
    counts <- tree@termStats$offspringCount
  } else {
    # Calculate counts and cache them
    counts <- countOffspring(tree, FALSE)
    tree@termStats$offspringCount <- counts
  }
  
  # Add names to the counts
  names(counts) <- tree@termNames
  
  # Include self if requested
  if(includeSelf) {
    counts <- counts + 1
  }
  
  # Return subset if terms specified
  if(!is.null(terms)) {
    if(is.numeric(terms)) {
      # Numeric indices
      return(counts[terms])
    } else {
      # Term names
      termIds <- match(terms, tree@termNames)
      if(any(is.na(termIds))) {
        stop("Some terms not found in the tree")
      }
      return(counts[termIds])
    }
  }
  
  return(counts)
}

#' Normalize relation types
#' 
#' @param relations Character vector of relation types
#' @return Normalized relation types
#' @author Kai Guo
#' @keywords internal
normalizeRelations <- function(relations) {
  relations <- gsub("isa", "is_a", relations)
  relations <- gsub(" ", "_", relations)
  return(relations)
}

#' Convert term names to term IDs
#' 
#' @param tree OntologyTree object
#' @param termNames Character vector of term names or numeric vector of term IDs
#' @param strict Whether to error on missing terms
#' @param addNames Whether to add term names to returned IDs
#' @return Integer vector of term IDs
#' @author Kai Guo
#' @export
term2id <- function(tree, termNames, strict = TRUE, addNames = FALSE) {
  # If numeric input, return as is
  if(is.numeric(termNames)) {
    termIds <- termNames
    
    # Single term case
  } else if(length(termNames) == 1) {
    termIds <- which(tree@termNames == termNames)
    if(length(termIds) == 0) {
      stop("Cannot find term: ", termNames)
    }
    
    # Multiple terms case
  } else {
    uniqueTerms <- unique(termNames)
    foundTerms <- tree@termNames %in% uniqueTerms
    termIds <- which(foundTerms)
    
    # Handle missing terms
    if(length(termIds) == 0) {
      stop("Cannot find any of these terms.")
    }
    if(length(termIds) != length(uniqueTerms)) {
      if(strict) {
        stop("Cannot find some of the terms in the tree.")
      } else {
        message("Removed ", length(uniqueTerms) - length(termIds), 
                " terms that cannot be found in the tree.")
      }
    }
    
    # Map found terms to original order
    termIds <- unname(structure(termIds, 
                                names = tree@termNames[termIds])[intersect(termNames, 
                                                                           tree@termNames[termIds])])
  }
  
  # Add names if requested
  if(addNames) {
    structure(termIds, names = tree@termNames[termIds])
  } else {
    termIds
  }
}

#' Calculate annotation-based Information Content
#' 
#' @param tree OntologyTree object
#' @param uniquify Whether to use unique annotations
#' @param useCache Whether to use cached values
#' @param verbose Whether to print messages
#' @return Numeric vector of IC values
#' @author Kai Guo
#' @export
calculateAnnotationIC <- function(tree, uniquify = TRUE, useCache = TRUE, verbose = TRUE) {
  if(verbose) {
    message("method: annotation")
  }
  
  # Check cache
  cacheKey <- if(uniquify) "annotation_unique" else "annotation"
  
  if(!useCache || is.null(tree@termStats[[cacheKey]])) {
    # Get annotation counts
    counts <- getAnnotations(tree, uniquify = uniquify, useCache = useCache)
    
    # Calculate IC
    maxCount <- max(counts)
    ic <- ifelse(counts == 0, NA_real_, -log(counts/maxCount))
    
    # Cache results
    if(useCache) {
      tree@termStats[[cacheKey]] <- ic
    }
    
    return(ic)
  }
  
  return(tree@termStats[[cacheKey]])
}

#' Calculate Wang's Information Content
#' 
#' @param tree OntologyTree object
#' @param contributionFactor Named numeric vector of semantic contribution factors
#' @param useCache Whether to use cached values
#' @param verbose Whether to print messages
#' @return Numeric vector of IC values
#' @author Kai Guo
#' @export
calculateWangIC <- function(tree, 
                            contributionFactor = c("is_a" = 0.8, "part_of" = 0.6), 
                            useCache = TRUE, 
                            verbose = TRUE) {
  if(verbose) {
    message("method: Wang")
  }
  
  # Check cache
  if(!useCache || is.null(tree@termStats$Wang)) {
    # Validate relations
    if(length(tree@edgeTypes) == 0) {
      stop("Relations not set in the ontology tree")
    }
    
    # Get relation levels
    relationLevels <- attr(tree@edgeTypes, "levels")
    
    # Validate contribution factors
    if(is.null(names(contributionFactor))) {
      stop("contributionFactor must be a named numeric vector corresponding to relations")
    }
    
    # Check contribution factor values
    if(any(contributionFactor >= 1)) {
      stop("All contribution factors must be less than 1")
    }
    
    # Normalize relation names
    names(contributionFactor) <- gsub("isa", "is_a", names(contributionFactor))
    names(contributionFactor) <- gsub(" ", "_", names(contributionFactor))
    
    # Extend contribution factors if relTree exists
    if(!is.null(tree@relTree)) {
      extendedFactors <- contributionFactor
      for(relationName in names(contributionFactor)) {
        if(relationName %in% tree@relTree@termNames) {
          offspring <- getOffspringCount(tree@relTree, relationName)
          if(length(offspring)) {
            extendedFactors[offspring] <- contributionFactor[relationName]
          }
        }
      }
      contributionFactor <- extendedFactors
    }
    
    # Validate all relations have factors
    if(length(setdiff(relationLevels, names(contributionFactor))) > 0) {
      stop("Contribution factors must be provided for all relations")
    }
    
    # Calculate IC values
    ic <- calWangIC(tree, unname(contributionFactor[relationLevels]))
    
    # Cache results
    if(useCache) {
      tree@termStats$Wang <- ic
    }
    
    return(ic)
  }
  
  return(tree@termStats$Wang)
}

#' Calculate Information Content for ontology terms
#' 
#' @param tree OntologyTree object
#' @param method IC calculation method ("universal", "annotation","wang" or "offspring")
#' @param useCache Whether to use cached values
#' @param contributionFactor Named numeric vector of semantic contribution factors
#' @param uniquify Whether to use unique annotations (annotation method)
#' @param verbose Whether to print messages
#' @return Named numeric vector of IC values
#' @author Kai Guo
#' @export
calculateIC <- function(tree, method = c("universal","annotation","wang","offspring"), useCache = TRUE,contributionFactor = c("is_a" = 0.8, "part_of" = 0.6),uniquify =TRUE,verbose = TRUE) {
  method <- match.arg(method)
  
  ic <- if(method == "universal") {
    calculateUniversal(tree, useCache=useCache, verbose=verbose)
  } else if(method == "annotation"){
    calculateAnnotationIC(tree, useCache=useCache, uniquify=uniquify, verbose=verbose)
  } else if(method == "wang"){
    calculateWangIC(tree, useCache=useCache, contributionFactor=contributionFactor, verbose=verbose)
  }else {
    calculateOffspring(tree, useCache=useCache, verbose=verbose)
  }
  
  names(ic) <- tree@termNames
  return(ic)
}

#' Find maximum value ancestor between term pairs
#' 
#' @param tree OntologyTree object
#' @param terms Character vector or numeric vector of term IDs
#' @param values Numeric vector of values for each term
#' @param returnNames Whether to return term names instead of IDs
#' @param distanceMethod Distance method ("longest" or "shortest")
#' @param verbose Whether to print progress messages
#' @return Matrix of ancestor IDs or names
#' @author Kai Guo
#' @export
findMaxAncestor <- function(tree, terms, values, returnNames = FALSE,
                            distanceMethod = c("longest", "shortest"),
                            verbose = TRUE) {
  # Input validation
  if(length(values) != tree@termCount) {
    stop("Length of 'values' must match the number of terms in the tree")
  }
  
  # Process terms input
  termIds <- if(is.character(terms)) {
    term2id(tree, terms)
  } else {
    terms
  }
  
  # Check for duplicates
  if(any(duplicated(termIds))) {
    stop("Terms must not be duplicated")
  }
  
  # Match distance method
  distanceMethod <- match.arg(distanceMethod)
  useLongestPath <- distanceMethod == "longest"
  
  # Calculate ancestor IDs
  ancestorMatrix <- if(verbose) {
    message("Calculating maximum value ancestors...")
    findMaxAncestorID(tree, termIds, values, useLongestPath)
  } else {
    findMaxAncestorID(tree, termIds, values, useLongestPath)
  }
  
  # Add dimension names
  termNames <- tree@termNames[termIds]
  dimnames(ancestorMatrix) <- list(termNames, termNames)
  
  # Return either IDs or names
  if(returnNames) {
    structure(tree@termNames[ancestorMatrix], 
              dim = dim(ancestorMatrix),
              dimnames = dimnames(ancestorMatrix))
  } else {
    ancestorMatrix
  }
}

#' Find Most Informative Common Ancestor (MICA) between terms
#' 
#' @param tree OntologyTree object
#' @param terms Character vector or numeric vector of term IDs
#' @param method Method for calculating Information Content
#' @param useCache Whether to use cached values
#' @param returnNames Whether to return term names instead of IDs
#' @param distanceMethod Distance method ("longest" or "shortest")
#' @param verbose Whether to print progress messages
#' @return Matrix of MICA term IDs or names
#' @author Kai Guo
#' @export
findMICAterm <- function(tree, 
                     terms, 
                     method="annotation",
                     useCache=TRUE,
                     returnNames = TRUE, 
                     distanceMethod = c("longest", "shortest"), 
                     verbose = TRUE) {
  # Calculate Information Content
  IC <- calculateIC(tree, method=method, useCache=useCache, verbose = verbose)
  
  # Find ancestor with maximum IC
  findMaxAncestor(tree, 
                  terms, 
                  values = IC, 
                  returnNames = returnNames,
                  distanceMethod = match.arg(distanceMethod),
                  verbose = verbose)
}

#' Calculate Most Informative Common Ancestor (MICA) values
#' 
#' @param tree OntologyTree object
#' @param terms Character vector or numeric vector of term IDs
#' @param method Method for calculating Information Content
#' @param useCache Whether to use cached values
#' @param verbose Whether to print progress messages
#' @return Matrix of MICA values
#' @author Kai Guo
#' @export
findMICA <- function(tree, 
                           terms, 
                           method="annotation",
                           useCache=TRUE,
                           verbose = TRUE) {
  # Calculate Information Content
  IC <- calculateIC(tree, method = method, useCache = useCache, verbose = verbose)
  # Process terms input
  termIds <- if(is.character(terms)) {
    term2id(tree, terms)
  } else {
    terms
  }
  
  # Check for duplicates
  if(any(duplicated(termIds))) {
    stop("Terms must not be duplicated")
  }
  # Find ancestor with maximum IC values
  # Add dimension names
  termNames <- tree@termNames[termIds]
  ancestorMatrix <- findMaxAncestorValues(tree, termIds, values = IC)
  dimnames(ancestorMatrix) <- list(termNames, termNames)
  ancestorMatrix
}

#' Validate and filter terms with zero annotations
#' 
#' @param tree OntologyTree object
#' @param termIds Term IDs to validate
#' @return Logical vector indicating which terms to keep
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
  
  # Return logical vector of terms to keep
  return(!zeroAnnots)
}

#' Convert OntologyTree to igraph object
#' 
#' @param tree An OntologyTree object to convert
#' @return An igraph object representing the ontology structure
#' @importFrom igraph graph_from_edgelist V set_edge_attr
#' @examples
#' \dontrun{
#' g <- convertToIgraph(tree)
#' }
#' @author Kai Guo
#' @export
convertToIgraph <- function(tree) {
  if(!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required")
  }
  
  # Extract child node information
  childMap <- tree@childMap
  children <- unlist(childMap)
  parents <- rep(seq_along(tree@termNames), 
                 times = vapply(childMap, length, FUN.VALUE = integer(1)))
  
  # Process edge types if available
  edgeTypes <- NULL
  if(length(tree@edgeTypes) > 0) {
    relationLevels <- attr(tree@edgeTypes, "levels")
    edgeTypes <- unlist(tree@edgeTypes)
    edgeTypes <- relationLevels[edgeTypes]
  }
  
  # Create igraph object
  graph <- igraph::graph_from_edgelist(cbind(parents, children))
  igraph::V(graph)$name <- tree@termNames
  
  # Add edge attributes if available
  if(!is.null(edgeTypes)) {
    graph <- igraph::set_edge_attr(graph, "relation", value = edgeTypes)
  }
  
  return(graph)
}

#' extract genes based on index 
#' @author Kai Guo
#' @return vector of gene id
extract_genes_from_terms <- function(tree, term_ids) {
  unique(tree@annotations$names[unlist(tree@annotations$list[term2id(tree, term_ids)])])
}
#' All Supported methods and the paramaters
#' @author Kai Guo
#' @export
AllSimMethod<-function(){
  d<-data.frame(paramaters=c("lin","resnik","faith","rel","simic","gogo","ic","wang","term","term","term","term","anc"),
                submethods=c("","","","","","","","","kappa","jaccard","dice","overlap",""),
                methods=c("Lin_1998","Resnik_1999" , "FaITH_2010" , "Relevance_2006","SimIC_2010" ,"universal" ,"Wang_2007","GOGO_2018" ,"Kappa" ,
                             "Jaccard" ,  "Dice" , "Overlap" , "Ancestor"))
  print(d)
}

#' All supported IC methods and the paramaters
#' @author Kai Guo
#' @export
AllICMethod<-function(){
  d<-c("offspring","annotation","universal","wang")
  d
}

#' Get Database Package Name for a Given Species
#'
#' This function returns the appropriate organism database package name for the specified species.
#' The package name is used for ID conversion (e.g. via AnnotationDbi) and is determined by a lookup
#' table. Allowed species are: "anopheles", "arabidopsis", "bovine", "celegans", "canine", "fly",
#' "zebrafish", "ecoli", "ecsakai", "chicken", "human", "mouse", "rhesus", "malaria", "chipm",
#' "rat", "toxoplasma", "streptomyces", "pig", "yeast", "xenopus", and "warm". If an unsupported
#' species is provided, the function returns "unsupported".
#'
#' @param species A character string specifying the species name (e.g., "human", "mouse").
#'
#' @return A character string containing the database package name corresponding to the input species.
#'
#' @details For example, for "human" the function returns "org.Hs.eg.db". This function is intended for
#' internal use in the package.
#'
#' @examples
#' .getdb("human")
#' .getdb("mouse")
#'
#' @keywords internal
.getdb <- function(species = "human") {
  species <- tryCatch(match.arg(species, c("anopheles", "arabidopsis", "bovine", "celegans", "canine", "fly", 
                                           "zebrafish", "ecoli", "ecsakai", "chicken", "human", "mouse",
                                           "rhesus", "malaria", "chipm", "rat",
                                           "toxoplasma", "streptomyces", "pig", "yeast", "xenopus", "warm")),
                      error = function(cond) { return("unsupported") })
  
  db_map <- c(
    anopheles    = "org.Ag.eg.db",
    arabidopsis  = "org.At.tair.db",
    bovine       = "org.Bt.eg.db",
    celegans     = "org.Ce.eg.db",
    canine       = "org.Cf.eg.db",
    fly          = "org.Dm.eg.db",
    zebrafish    = "org.Dr.eg.db",
    ecoli        = "org.EcK12.eg.db",
    ecsakai      = "org.EcSakai.eg.db",
    chicken      = "org.Gg.eg.db",
    human        = "org.Hs.eg.db",
    mouse        = "org.Mm.eg.db",
    rhesus       = "org.Mmu.eg.db",
    malaria      = "org.Pf.plasmo.db",
    chipm        = "org.Pt.eg.db",
    rat          = "org.Rn.eg.db",
    toxoplasma   = "org.Tgondii.eg.db",
    streptomyces = "org.Sco.eg.db",
    pig          = "org.Ss.eg.db",
    yeast        = "org.Sc.sgd.db",
    xenopus      = "org.Xl.eg.db",
    warm         = NA_character_
  )
  
  return(db_map[species])
}

#' Get Abbreviated Species Code
#'
#' This function returns the abbreviated species code for the specified species.
#' The abbreviated code is used in KEGG queries and internal references.
#'
#' @param species A character string specifying the species name (e.g., "human", "arabidopsis").
#'
#' @return A character string containing the abbreviated species code.
#'
#' @details For example, "human" is mapped to "hsa", "mouse" to "mmu", and "fly" to "dme". This
#' function uses a lookup table and is intended for internal use in the package.
#'
#' @examples
#' .getspeices("human")
#' .getspeices("arabidopsis")
#'
#' @keywords internal
.getspeices <- function(species = "human") {
  species <- tryCatch(match.arg(species, c("anopheles", "arabidopsis", "bovine", "celegans", "canine", "fly", 
                                           "zebrafish", "ecoli", "ecsakai", "chicken", "human", "mouse", "rhesus",
                                           "malaria", "chipm", "rat", "toxoplasma", "sco", "pig", "yeast",
                                           "xenopus", "warm")),
                      error = function(cond){ return("unsupported") })
  
  sp_map <- c(
    anopheles    = "aga",
    arabidopsis  = "ath",
    bovine       = "bta",
    celegans     = "cel",
    canine       = "cfa",
    fly          = "dme",
    zebrafish    = "dre",
    ecoli        = "eco",
    ecsakai      = "ecs",
    chicken      = "gga",
    human        = "hsa",
    mouse        = "mmu",
    rhesus       = "mcc",
    malaria      = "pfa",
    chipm        = "ptr",
    rat          = "rno",
    toxoplasma   = "tgondii",
    sco          = "sco",
    pig          = "ssc",
    yeast        = "sce",
    xenopus      = "xla",
    warm         = "warm"
  )
  
  return(sp_map[species])
}


