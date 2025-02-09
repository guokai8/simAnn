#' Build Gene Ontology Tree
#' 
#' @param species A character string specifying the species name (e.g., "human", "mouse").
#' @param namespace GO namespace ('BP', 'CC', or 'MF')
#' @param relations Character vector of relation types to include
#' @param builtin Logical. If \code{TRUE} and \code{species} is one of the built-in types (\code{"human"},
#'   \code{"mouse"}, \code{"rat"}, or \code{"ko"}), the function loads a prebuilt GO tree from
#'   the built-in data. Default is \code{TRUE}.
#' @param evidenceCode Vector of evidence codes to filter
#' @param includeAlt Include alternative terms
#' @param verbose Print progress messages
#' @return OntologyTree object
#' @examples
#' # Basic BP ontology
#' go_bp <- buildGOTree("BP")
#' 
#' # MF ontology with specific relations
#' go_mf <- buildGOTree("MF", relations = c("part_of", "regulates"))
#' 
#' # With organism annotations
#' \dontrun{
#' tree <- buildGOTree(species="human",namespace="BP", keytype="SYMBOL",builtin=TRUE)
#' }
#' @import GO.db AnnotationDbi
#' @author Kai Guo
#' @export
buildGOTree <- function(
    species = "human",
    namespace = c("BP", "CC", "MF"),
    keytype = "SYMBOL",
    relations = "part_of",
    builtin=TRUE,
    evidenceCode = NULL,
    includeAlt = FALSE,
    verbose = TRUE
) {
  # Input validation
  if(species %in% c("human", "mouse", "rat") && isTRUE(builtin)) {
    data(goTree)  # Assumes built-in data object 'keggTree' is available.
    return(goTree[[species]])
  }
  orgDb <- .getdb(species = species)
  namespace <- match.arg(namespace)
  validateInputsGO(namespace, relations, orgDb, evidenceCode)
  
  # Get relationships
  df <- getGORelationships(namespace, verbose)
  
  # Process relations
  df <- processRelations(df, relations, verbose)
  
  # Get annotations
  annotations <- if(!is.null(orgDb)) {
    getAnnotationsdb(orgDb, namespace=namespace, evidenceCode=evidenceCode,keytype=keytype, verbose=verbose)
  }
  
  # Build the tree structure
  tree <- buildTreeStructure(
    df = df,
    annotations = annotations,
    namespace = namespace,
    includeAlt = includeAlt,
    verbose = verbose
  )
  
  # Add metadata and return
  addGOMetadata(tree)
}

#' Get GO relationships from appropriate database
#' @importFrom GO.db GOBPCHILDREN GOCCCHILDREN GOMFCHILDREN
#' @keywords internal
getGORelationships <- function(namespace, verbose) {
  if(verbose) message("Fetching GO ", namespace, " relationships...")
  
  tryCatch({
    df <- switch(namespace,
                 BP = AnnotationDbi::toTable(GO.db::GOBPCHILDREN),
                 CC = AnnotationDbi::toTable(GO.db::GOCCCHILDREN),
                 MF = AnnotationDbi::toTable(GO.db::GOMFCHILDREN)
    )
    
    # Standard column names
    colnames(df) <- c("child", "parent", "relation")
    return(df)
  }, 
  error = function(e) {
    stop("Failed to fetch GO relationships: ", e$message)
  })
}

#' Process and filter relations
#' @keywords internal
processRelations <- function(df, relations, verbose) {
  # Basic standardization of relation types in data
  l = df$relation == "isa"
  df$relation[l] = "is_a"
  df$relation = gsub(" ", "_", df$relation)
  
  # Always include is_a
  relations <- unique(c("is_a", relations))
  
  # Handle regulates relationships
  if("regulates" %in% relations) {
    relations <- c(relations, "negatively_regulates", "positively_regulates")
  }
  
  if(verbose) {
    message("Including relations: ", paste(relations, collapse = ", "))
  }
  
  # Filter by relations
  df[df$relation %in% relations, , drop = FALSE]
}

#' Get and process annotations from organism database
#' 
#' @param orgDb Organism database object or name
#' @param namespace GO namespace ("BP", "MF", "CC") 
#' @param evidenceCode Vector of evidence codes to include
#' @param keytype Type of gene IDs ("SYMBOL", "ENTREZID", "REFSEQ", "ENSEMBL")
#' @param verbose Print progress messages
#' @keywords internal
getAnnotationsdb <- function(orgDb, namespace = "BP", evidenceCode = NULL, 
                           keytype = "SYMBOL", verbose = TRUE) {
  if(verbose) message("Processing annotations...")
  
  # Validate keytype
  valid_keytypes <- c("SYMBOL", "ENTREZID", "REFSEQ", "ENSEMBL") 
  keytype <- match.arg(keytype, valid_keytypes)
  
  # Load database
  db <- loadOrgDb(orgDb)
  
  tryCatch({
    # Get annotations with specified keytype
    tb <- suppressMessages(AnnotationDbi::select(db,
                                                 keys = AnnotationDbi::keys(db, keytype = keytype),
                                                 columns = c(keytype, "GO", "EVIDENCE", "ONTOLOGY"),
                                                 keytype = keytype))
    
    # Filter namespace
    tb <- tb[tb$ONTOLOGY == namespace, , drop = FALSE]
    
    # Filter evidence codes
    if(!is.null(evidenceCode)) {
      tb <- filterByEvidence(tb, evidenceCode)
    }
    
    # Process using keytype column
    processAnnotationTable(tb, id_col = keytype)
    
  }, error = function(e) {
    stop("Failed to process annotations: ", e$message)
  })
}

#' Process annotation table to list format
#' @keywords internal  
processAnnotationTable <- function(tb, id_col = "SYMBOL") {
  split(tb[[id_col]], tb$GO)
}

#' Load organism database
#' @author Kai Guo
#' @keywords internal
loadOrgDb <- function(orgDb) {
  if(is.character(orgDb)) {
    if(!requireNamespace(orgDb, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required", orgDb))
    }
    return(getFromNamespace(orgDb, ns = orgDb))
  }
  orgDb
}

#' Filter annotation table by evidence codes
#' @author Kai Guo
#' @keywords internal
filterByEvidence <- function(tb, evidenceCode) {
  filtered <- tb[tb$EVIDENCE %in% evidenceCode, , drop = FALSE]
  if(nrow(filtered) == 0) {
    stop("No annotations remain after evidence code filtering")
  }
  filtered
}

#' Process annotation table to list format 
#' @param tb Annotation table with gene IDs and GO terms
#' @param id_col Column name for gene IDs 
#' @author Kai Guo
#' @keywords internal
processAnnotationTable <- function(tb, id_col = "SYMBOL") {
  tb <- tb[, c(which(colnames(tb) == id_col), 
               which(colnames(tb) == "GO")), drop = FALSE]
  
  annotations <- split(tb[,1], tb[,2])  
  lapply(annotations, function(x) unique(as.character(x)))
}

#' Build the core tree structure
#' @keywords internal
#' @author Kai Guo
buildTreeStructure <- function(df, annotations, namespace, includeAlt, verbose) {
  # Build relation tree
  relTree <- buildOntologyTree(
    parentTerms = c("regulates", "regulates"),
    childTerms = c("negatively_regulates", "positively_regulates")
  )
  
  # Get GO.db version for source tracking
  goVersion <- packageVersion("GO.db")
  
  # Create main ontology tree
  buildOntologyTree(
    parentTerms = df$parent,
    childTerms = df$child,
    edgeTypes = df$relation,
    relTree = relTree,
    annotations = annotations,
    source = sprintf("GO %s / GO.db %s", namespace, goVersion)  )
}

#' Add GO term metadata
#' @keywords internal
#' @author Kai Guo
addGOMetadata <- function(tree) {
  # Fetch term information
  goTerms <- GO.db::GOTERM[tree@termNames]
  
  # Build metadata dataframe
  metadata <- data.frame(
    TermID = AnnotationDbi::GOID(goTerms),
    TermAnnot = AnnotationDbi::Term(goTerms),
    definition = AnnotationDbi::Definition(goTerms),
    namespace = AnnotationDbi::Ontology(goTerms),
    row.names = tree@termNames,
    stringsAsFactors = FALSE
  )
  
  # Add metadata to tree
  tree@metadata <- metadata
  tree
}

#' Validate function inputs
#' @keywords internal
#' @author Kai Guo
validateInputsGO <- function(namespace, relations, orgDb, evidenceCode) {
  if(!requireNamespace("GO.db", quietly = TRUE)) {
    stop("Package 'GO.db' is required")
  }
  
  if(!is.null(evidenceCode) && is.null(orgDb)) {
    stop("Evidence codes can only be used when orgDb is specified")
  }
  
  if(!is.null(relations)) {
    if(!is.character(relations)) {
      stop("Relations must be a character vector")
    }
    if(any(is.na(relations))) {
      warning("NA values in relations will be removed")
    }
  }
}
#' count the annotation numbers
#' @author Kai Guo
getAnnotations <- function(tree, terms = NULL, uniquify = TRUE, useCache = TRUE) {
  # Get from cache or calculate
  if(uniquify) {
    if(!useCache || is.null(tree@termStats$n_annotations_unique)) {
      if(length(tree@annotations$list) == 0) stop("No annotations found in tree")
      counts <- countAnnotations(tree, TRUE)
      if(useCache) tree@termStats$n_annotations_unique <- counts
    } else {
      counts <- tree@termStats$n_annotations_unique
    }
  } else {
    if(!useCache || is.null(tree@termStats$n_annotations)) {
      if(length(tree@annotations$list) == 0) stop("No annotations found in tree")
      counts <- countAnnotations(tree, FALSE)
      if(useCache) tree@termStats$n_annotations <- counts
    } else {
      counts <- tree@termStats$n_annotations
    }
  }
  
  # Add term names
  names(counts) <- tree@termNames
  
  # Return subset if terms provided
  if(!is.null(terms)) {
    if(!is.numeric(terms)) {
      terms <- term2id(tree, terms)
    }
    return(counts[terms])
  }
  
  return(counts)
}