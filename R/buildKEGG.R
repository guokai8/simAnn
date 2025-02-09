#' Convert IDs using an organism database
#'
#' This helper function converts a vector of keys from one type to another using AnnotationDbi::mapIds.
#'
#' @param species A string specifying the species code (e.g., "human","mouse","rat").
#' @param keys A vector of IDs to convert.
#' @param fkeytype The key type of the input IDs.
#' @param tkeytype The target key type (e.g., "SYMBOL").
#' @return A vector of converted IDs.
#' @keywords internal
.idc <- function(species, keys, fkeytype, tkeytype) {
  dbname <- .getdb(species = species)
  orgDb <- loadOrgDb(dbname)
  unlist(mapIds(orgDb, keys = as.vector(keys),
                column = tkeytype,
                keytype = fkeytype,
                multiVals = "first"))
}

#' Get Gene–Pathway Mapping from KEGG
#'
#' This function retrieves pathway information from KEGG by calling
#' keggLink("pathway", species) and returns a data frame with two columns:
#'   - GeneID: The gene (or enzyme) identifier.
#'   - PathwayID: The associated KEGG pathway ID.
#'
#' Optionally, the gene IDs can be converted using the internal helper function .idc.
#' @importFrom KEGGREST keggLink
#' @param species A string specifying the species code (e.g., "human").
#'                Default is "hsa".
#' @param keytype A character string specifying the target key type (e.g., "SYMBOL").
#'                Default is "SYMBOL".
#'
#' @return A data frame with two columns: GeneID and Level3
#'
#' @examples
#' \dontrun{
#'   # For human, converting KEGG gene IDs to SYMBOL:
#'   mapping <- getGenePathwayMapping("human", keytype =  "SYMBOL")
#'   head(mapping)
#' }
#'
#' @export
getGenePathwayMapping <- function(species = "human", keytype = "ENTREZID") {
  
  if(species=="ko"){
    spe="ko"
    keytype = "ENTREZID"
  }else{
    spe = .getspeices(species = species)
  }
  tmp <- keggLink("pathway", spe)
  tmp <- substr(tmp,9,13)
  names(tmp) <- sub('.*:','',names(tmp))
  # Construct a data frame with two columns.
  df <- data.frame(GeneID = names(tmp),
                   Level3 = as.vector(tmp),
                   stringsAsFactors = FALSE)
  
  if(keytype != "ENTREZID") {
    df$GeneID <- .idc(species, df$GeneID, fkeytype = "ENTREZID", tkeytype = keytype)
  }
  df <- na.omit(df)
  return(df)
}


#' Build KEGG Ontology Tree
#'
#' This function retrieves KEGG BRITE data using \code{KEGGREST::keggGet} with a query of the form
#' \code{"br:<species>00001"}. The default species is \code{"ko"}. When the \code{builtin} parameter is
#' set to \code{TRUE} and the supplied \code{species} is one of \code{"human"}, \code{"mouse"},
#' \code{"macaca"}, \code{"rat"}, or \code{"ko"}, a prebuilt KEGG ontology tree is loaded from the
#' built-in data. Otherwise, the function downloads the KEGG BRITE data, parses its hierarchical structure,
#' and constructs the ontology tree.
#'
#' The hierarchical structure is parsed into a data frame that includes the following columns:
#' \itemize{
#'   \item \code{GeneID} - the identifier of the leaf node (an enzyme K-code or gene ID),
#'   \item \code{Anno} - the annotation for the leaf node,
#'   \item \code{Level3} - the Level 3 code,
#'   \item \code{Level3Anno} - the Level 3 annotation,
#'   \item \code{Level2} - the Level 2 code,
#'   \item \code{Level2Anno} - the Level 2 annotation,
#'   \item \code{Level1} - the Level 1 code,
#'   \item \code{Level1Anno} - the Level 1 annotation.
#' }
#'
#' Please note the following adjustments made for internal consistency:
#' \itemize{
#'   \item The KEGG code for Drug Development is designated as \code{09200}.
#'   \item Temporary codes were assigned to Level2 (e.g., \code{09113} and \code{09114} for temporary use).
#'   \item To resolve duplicate entries, the following modifications were applied:
#'     \item \code{541X} represents Biosynthesis of various nucleotide sugars,
#'     \item \code{541} represents O-Antigen nucleotide sugar biosynthesis,
#'     \item \code{710X} represents Carbon fixation by Calvin cycle,
#'     \item \code{710} represents Carbon fixation in photosynthetic organisms,
#'     \item \code{720X} represents Other carbon fixation pathways, and
#'     \code{720} represents Carbon fixation pathways in prokaryotes.
#' }
#'
#' The function then calls the existing \code{buildOntologyTree()} function (which expects a data frame
#' with columns \code{"GeneID"}, \code{"Level3"}, \code{"Level2"}, and \code{"Level1"}) to construct the
#' final \code{OntologyTree} object.
#'
#' Metadata is generated solely from the level information (Levels A, B, and C) and is stored as a data frame
#' with columns \code{id}, \code{name}, and \code{namespace}. (Leaf node data, including the K-codes or gene IDs,
#' do not appear in the metadata; their annotations are maintained in the \code{tree@annotations} slot.)
#'
#' Finally, the species information is appended as an attribute to the metadata.
#'
#' @param species A string specifying the species prefix (e.g., \code{"ko"} or \code{"hsa"}). Default is \code{"ko"}.
#' @param keytype For species \code{"hsa"}, a character vector specifying the gene identifier type. Allowed values
#'   are \code{"SYMBOL"} and \code{"ENTREZID"}. The default is \code{"SYMBOL"}.
#' @param builtin Logical. If \code{TRUE} and \code{species} is one of the built-in types (\code{"human"},
#'   \code{"mouse"}, \code{"rat"}, or \code{"ko"}), the function loads a prebuilt KEGG tree from
#'   the built-in data. Default is \code{TRUE}.
#' @param verbose Logical. When \code{TRUE}, the function prints progress messages. Default is \code{TRUE}.
#'
#' @return An \code{OntologyTree} object. The object contains:
#'   \itemize{
#'     \item A \code{metadata} slot holding a data frame with columns \code{id}, \code{name}, and \code{namespace}.
#'     \item Leaf node annotation information stored in \code{tree@annotations}.
#'     \item An attribute on the metadata specifying the species (in the form \code{"<species>00001"}).
#'   }
#'
#' @importFrom KEGGREST keggLink
#' @export
buildKEGGTree <- function(species = "ko", keytype = c("SYMBOL","ENTREZID"), builtin=TRUE, verbose = TRUE) {
  if (verbose) message("Fetching data using ", species, "...")
  if(species %in% c("human", "mouse", "rat") && isTRUE(builtin)) {
    data(keggTree)  # Assumes built-in data object 'keggTree' is available.
    return(keggTree[[species]])
  }
  pathway <- getGenePathwayMapping(species = species,keytype = "SYMBOL")
  ### load the pathway level information
  data("leaf_df")
  leaf <- subset(leaf_df,Level3%in%unique(pathway$Level3))
  dd <- merge(pathway,leaf,by.x="Level3",by.y="Level3",all.x=T)
  # Build annotation information for leaf nodes.
  leaf_meta <- unique(dd[, c("GeneID", c("Level3","Level2","Level1"))])
  annotations <- c(split(leaf_meta$GeneID,leaf_meta$Level1),split(leaf_meta$GeneID,leaf_meta$Level2),split(leaf_meta$GeneID,leaf_meta$Level3))
  tree <- buildOntologyTree(parentTerms = c(leaf$Level1,leaf$Level2),childTerms = c(leaf$Level2,leaf$Level3),edgeTypes=rep("is_a",nrow(leaf)*2),annotations = annotations,source = "KEGG", verbose = verbose)
  # Build metadata from levels A, B, and C.
  level1_meta <- leaf[, c("Level1", "Level1Anno")]
  level1_df <- if(nrow(level1_meta) > 0) {
    data.frame(id = level1_meta$Level1, name = level1_meta$Level1Anno,
               namespace = "Level1", stringsAsFactors = FALSE)
  } else { data.frame() }
  level2_meta <- leaf[, c("Level2", "Level2Anno")]
  level2_df <- if(nrow(level2_meta) > 0) {
    data.frame(id = level2_meta$Level2, name = level2_meta$Level2Anno,
               namespace = "Level2", stringsAsFactors = FALSE)
  } else { data.frame() }
  level3_meta <- leaf[, c("Level3", "Level3Anno")]
  level3_df <- if(nrow(level3_meta) > 0) {
    data.frame(id = level3_meta$Level3, name = level3_meta$Level3Anno,
               namespace = "Level3", stringsAsFactors = FALSE)
  } else { data.frame() }
  metadata_df <- rbind(level1_df, level2_df, level3_df)
  
  all_terms <- tree@termNames
  final_metadata <- data.frame(
    id = all_terms,
    name = NA_character_,
    namespace = NA_character_,
    row.names = all_terms,
    stringsAsFactors = FALSE
  )
  for(i in 1:nrow(metadata_df)) {
    term_id <- metadata_df$id[i]
    if(term_id %in% all_terms) {
      final_metadata[term_id, "id"] <- term_id
      final_metadata[term_id, "name"] <- metadata_df$name[i]
      final_metadata[term_id, "namespace"] <- metadata_df$namespace[i]
    }
  }
  tree@metadata <- final_metadata
  # Build metadata from levels A, B, and C (exclude leaf nodes).
  level1_meta <- leaf[, c("Level1", "Level1Anno")]
  if (nrow(level1_meta) > 0) {
    level1_df <- data.frame(id = level1_meta$Level1,
                            name = level1_meta$Level1Anno,
                            namespace = "Level1",
                            stringsAsFactors = FALSE)
  } else {
    level1_df <- data.frame()
  }
  level2_meta <- leaf[, c("Level2", "Level2Anno")]
  if (nrow(level2_meta) > 0) {
    level2_df <- data.frame(id = level2_meta$Level2,
                            name = level2_meta$Level2Anno,
                            namespace = "Level2",
                            stringsAsFactors = FALSE)
  } else {
    level2_df <- data.frame()
  }
  level3_meta <- leaf[, c("Level3", "Level3Anno")]
  if (nrow(level3_meta) > 0) {
    level3_df <- data.frame(id = level3_meta$Level3,
                            name = level3_meta$Level3Anno,
                            namespace = "Level3",
                            stringsAsFactors = FALSE)
  } else {
    level3_df <- data.frame()
  }
  metadata_df <- rbind(level1_df, level2_df, level3_df)
  
  all_terms <- tree@termNames
  final_metadata <- data.frame(
    TermID = all_terms,
    TermAnnot = NA_character_,
    namespace = NA_character_,
    row.names = all_terms,
    stringsAsFactors = FALSE
  )
  for(i in seq_len(nrow(metadata_df))) {
    term_id <- metadata_df$id[i]
    if(term_id %in% all_terms) {
      final_metadata[term_id, "TermID"] <- term_id
      final_metadata[term_id, "TermAnnot"] <- metadata_df$name[i]
      final_metadata[term_id, "namespace"] <- metadata_df$namespace[i]
    }
  }
  tree@metadata <- final_metadata
  # Attach species information as an attribute to the metadata.
  attr(tree@metadata, "species") <- paste0(species, "00001")
  
  if (verbose) message("Ontology tree built successfully.")
  return(tree)
}


.getspeices<-function(species="human"){
  species=tryCatch(match.arg(species,c("anopheles","arabidopsis","bovine","celegans","canine","fly","zebrafish",
                                       "ecoli","ecsakai","chicken","human","mouse","rhesus","malaria","chipm","rat",
                                       "toxoplasma","sco","pig","yeast","xenopus","warm")),
                   error=function(cond){return("unsupported")})
  if (species == "anopheles") {
    species <- "aga"
  } else if (species == "arabidopsis") {
    species <- "ath"
  } else if (species == "bovine") {
    species <- "bta"
  } else if (species == "canine") {
    species <- "cfa"
  } else if (species == "chicken") {
    species <- "gga"
  } else if (species == "chipm") {
    species <- "ptr"
  } else if (species == "ecolik12") {
    species <- "eco"
  } else if (species == "ecsakai") {
    species <- "ecs"
  } else if (species == "fly") {
    species <- "dme"
  } else if (species == "human") {
    species <- "hsa"
  } else if (species == "malaria") {
    species <- "pfa"
  } else if (species == "mouse") {
    species <- "mmu"
  } else if (species == "pig") {
    species <- "ssc"
  } else if (species == "rat") {
    species <- "rno"
  } else if (species == "rhesus") {
    species <- "mcc"
  } else if (species == "worm" || species == "celegans") {
    species <- "cel"
  } else if (species == "xenopus") {
    species <- "xla"
  } else if (species == "yeast") {
    species <- "sce"
  } else if (species =="streptomyces"){
    species <- "sco"
  } else if (species == "zebrafish") {
    species <- "dre"
  } else {
    species <- "ko"
  }
  return(species)
}
