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
#' @description
#' Constructs a hierarchical ontology tree from KEGG BRITE data, either by loading
#' prebuilt data for supported species or by retrieving and parsing KEGG data
#' dynamically. The function handles the complex hierarchical structure of KEGG
#' pathways and creates a standardized ontology representation.
#'
#' @details
#' The function processes KEGG data in the following steps:
#' 
#' 1. Data Source:
#'    * Uses prebuilt data for supported species if builtin=TRUE
#'    * Retrieves data via KEGGREST::keggGet using "br:<species>00001" format
#' 
#' 2. Hierarchical Structure:
#'    * Parses KEGG BRITE hierarchy into a standardized format
#'    * Creates a data frame with consistent level information
#'    * Handles special cases and duplicates in KEGG pathways
#' 
#' 3. Internal Adjustments:
#' \describe{
#'   \item{Drug Development}{Assigned code 09200}
#'   \item{Level2 Temporary Codes}{Uses codes like 09113, 09114}
#'   \item{Duplicate Handling}{
#'     * 541X: Biosynthesis of various nucleotide sugars
#'     * 541: O-Antigen nucleotide sugar biosynthesis
#'     * 710X: Carbon fixation by Calvin cycle
#'     * 710: Carbon fixation in photosynthetic organisms
#'     * 720X: Other carbon fixation pathways
#'     * 720: Carbon fixation pathways in prokaryotes
#'   }
#' }
#'
#' @param species Character string specifying the species:
#' \describe{
#'   \item{ko}{KEGG Orthology (default)}
#'   \item{human}{Human}
#'   \item{mouse}{Mouse}
#'   \item{rat}{Rat}
#' }
#' @param keytype Character string specifying the gene identifier type for human data:
#' \describe{
#'   \item{SYMBOL}{Gene symbols (default)}
#'   \item{ENTREZID}{Entrez gene IDs}
#' }
#' @param builtin Logical indicating whether to use prebuilt data for supported
#'        species. Default: TRUE.
#' @param verbose Logical indicating whether to print progress messages.
#'        Default: TRUE.
#'
#' @return An OntologyTree object containing:
#' \describe{
#'   \item{metadata}{Data frame with columns:
#'     \itemize{
#'       \item id: Unique identifier for each node
#'       \item name: Descriptive name of the node
#'       \item namespace: Ontology namespace information
#'     }
#'   }
#'   \item{annotations}{Leaf node annotations including:
#'     \itemize{
#'       \item GeneID: Enzyme K-code or gene ID
#'       \item Anno: Node annotation
#'       \item Level3: Level 3 pathway code
#'       \item Level3Anno: Level 3 pathway annotation
#'       \item Level2: Level 2 pathway code
#'       \item Level2Anno: Level 2 pathway annotation
#'       \item Level1: Level 1 pathway code
#'       \item Level1Anno: Level 1 pathway annotation
#'     }
#'   }
#'   \item{attributes}{Species information in format "<species>00001"}
#' }
#'
#' @examples
#' \dontrun{
#' # Build KEGG tree for general orthology
#' ko_tree <- buildKEGGTree()
#'
#' # Build human KEGG tree with gene symbols
#' human_tree <- buildKEGGTree(
#'   species = "human",
#'   keytype = "SYMBOL",
#'   builtin = TRUE
#' )
#'
#' # Build mouse KEGG tree with verbose output
#' mouse_tree <- buildKEGGTree(
#'   species = "mouse",
#'   builtin = TRUE,
#'   verbose = TRUE
#' )
#'
#' # Build rat KEGG tree with downloaded data
#' rat_tree <- buildKEGGTree(
#'   species = "rat",
#'   builtin = FALSE,
#'   verbose = TRUE
#' )
#' }
#'
#' @seealso
#' \code{\link[KEGGREST]{keggGet}} for KEGG data retrieval
#' \code{\link{buildOntologyTree}} for ontology tree construction
#'
#' @references
#' Kanehisa, M., & Goto, S. (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes.
#' Nucleic Acids Research, 28(1), 27-30.
#'
#' @importFrom KEGGREST keggLink
#' @return OntologyTree object
#' @author Kai Guo
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
