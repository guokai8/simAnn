#' Helper function for data aggregation
#' @noRd
aggregateColumns <- function(df, col_name, method) {
  if (!col_name %in% names(df)) return(NA)
  values <- df[[col_name]]
  if (length(values) == 0) return(NA)
  
  switch(method,
         "mean" = mean(values, na.rm = TRUE),
         "max" = max(values, na.rm = TRUE),
         "min" = min(values, na.rm = TRUE),
         "sum" = sum(values, na.rm = TRUE),
         "concat" = paste(unique(values), collapse = "; "),
         stop("Unsupported aggregation method: ", method))
}

#' Helper function for fast string splitting
#' @noRd
fastSplit <- function(x, split) {
  if (is.factor(x)) x <- as.character(x)
  strsplit(x, split, fixed = TRUE)
}


#' Cluster Complex Analysis for Protein Complexes
#'
#' @description
#' Integrates enrichment analysis results with protein complex annotations to identify
#' significantly enriched functional groups within protein complexes. The function supports
#' both simple GO ID lists and detailed enrichment analysis results as input.
#'
#' @param enrich_data Either a vector of GO IDs or a data frame containing enrichment
#'        analysis results. If a data frame, it must contain GO ID and gene identifier columns.
#' @param species A character string specifying the species ("human" or "mouse").
#' @param ontology_type A character string specifying the GO ontology type ("BP", "CC", or "MF").
#'        BP: Biological Process, CC: Cellular Component, MF: Molecular Function.
#' @param gene_delim Character string for splitting gene lists. Default: ",".
#' @param enrich_go_col Column name for GO IDs in enrich_data. Default: "Annot".
#' @param enrich_gene_col Column name for genes in enrich_data. Default: "GeneID".
#' @param merged_go_col Column name for GO IDs in complex data. Default: "functions_go_id".
#' @param merged_gene_col Column name for genes in complex data. Default: "subunits_gene_name".
#' @param merged_complex_col Column name for complex names. Default: "ComplexName".
#' @param column_config List defining column aggregation rules. Each element should contain:
#'        \itemize{
#'          \item col: Column name to aggregate
#'          \item agg: Aggregation method ("mean", "max", "min", "sum", or "concat")
#'        }
#' @param join_method Join type for combining enrichment and complex data ("inner" or "left").
#'        Default: "inner".
#' @param fisher_alternative Alternative hypothesis for Fisher's exact test
#'        ("two.sided", "greater", or "less"). Default: "two.sided".
#' @param include_unmatched Include entries without complex matches. Default: TRUE.
#' @param adjust_p Perform multiple testing correction. Default: TRUE.
#' @param p_adjust_method Method for p-value adjustment (see ?p.adjust). Default: "BH".
#' @param verbose Print progress messages. Default: FALSE.
#'
#' @return A list containing two data frames:
#' \describe{
#'   \item{summary}{Aggregated results per complex including:
#'     \itemize{
#'       \item ComplexName: Name of the protein complex
#'       \item GO_ids: Matched GO term IDs
#'       \item Matches: Number of matching GO terms
#'       \item Background: Total GO terms for the complex
#'       \item Ratio: Matches/Background ratio
#'       \item P_Com: Fisher's exact test p-value
#'       \item adj_P_Com: Adjusted p-value (if adjust_p=TRUE)
#'       \item Additional columns based on column_config
#'     }
#'   }
#'   \item{matchedAnnot}{Detailed matching records between enrichment and complex data}
#' }
#'
#' @examples
#' \dontrun{
#' # Vector mode - using GO IDs only
#' go_vector <- c("GO:0006915", "GO:0012501")
#' result1 <- clusterComplex(go_vector, 
#'                          species = "human",
#'                          ontology_type = "BP")
#'
#' # Table mode - using enrichment results
#' enrich_df <- data.frame(
#'   Annot = c("GO:0006915", "GO:0012501"),
#'   GeneID = c("CASP3,BCL2", "CASP3,BAX"),
#'   Pvalue = c(0.001, 0.002)
#' )
#' result2 <- clusterComplex(enrich_df,
#'                          species = "mouse",
#'                          ontology_type = "BP",
#'                          verbose = TRUE)
#' }
#'
#' @note
#' - The function automatically loads appropriate complex data based on species
#' - For table mode, gene lists in enrich_gene_col are automatically split
#' - P-values are calculated using Fisher's exact test
#' - Multiple testing correction is applied across all complexes when adjust_p=TRUE
#'
#' @seealso
#' \code{\link[stats]{fisher.test}} for Fisher's exact test details
#' \code{\link[stats]{p.adjust}} for p-value adjustment methods
#'
#' @importFrom stats fisher.test p.adjust
#' @return a list include summary and matchedAnnot data.frame
#' @export
clusterComplex <- function(enrich_data,
                           species = c("human", "mouse"),
                           ontology_type = c("BP", "CC", "MF"),
                           # Core parameters
                           gene_delim = ",",
                           enrich_go_col = "Annot",
                           enrich_gene_col = "GeneID",
                           merged_go_col = "functions_go_id",
                           merged_gene_col = "subunits_gene_name",
                           merged_complex_col = "ComplexName",
                           # Column configurations
                           column_config = list(
                             Term = list(col = "Term", agg = "concat"),
                             Annotated = list(col = "Annotated", agg = "mean"),
                             Significant = list(col = "Significant", agg = "mean"),
                             RichFactor = list(col = "RichFactor", agg = "max"),
                             FoldEnrichment = list(col = "FoldEnrichment", agg = "max"),
                             zscore = list(col = "zscore", agg = "max"),
                             Pvalue = list(col = "Pvalue", agg = "min"),
                             Padj = list(col = "Padj", agg = "min")
                           ),
                           # Options
                           join_method = c("inner", "left"),
                           fisher_alternative = c("two.sided", "greater", "less"),
                           include_unmatched = TRUE,
                           adjust_p = TRUE,
                           p_adjust_method = "BH",
                           verbose = FALSE) {
  
  # Input validation
  species <- match.arg(species)
  ontology_type <- match.arg(ontology_type)
  join_method <- match.arg(join_method)
  fisher_alternative <- match.arg(fisher_alternative)
  
  # Helper function for logging
  logMsg <- function(msg) {
    if (verbose) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg))
  }
  
  # Load and validate data
  tryCatch({
    logMsg("Loading complex data")
    if (species == "human") {
      data("hComplex", envir = environment())
      merged_data <- get("hComplex", envir = environment())
    } else {
      data("mComplex", envir = environment())
      merged_data <- get("mComplex", envir = environment())
    }
  }, error = function(e) {
    stop("Failed to load complex data for species: ", species, "\nError: ", e$message)
  })
  
  # Filter by ontology type
  merged_data <- subset(merged_data, Ontology_type == ontology_type)
  if (nrow(merged_data) == 0) {
    stop("No data available for specified ontology type: ", ontology_type)
  }
  
  # Determine operation mode
  isVectorMode <- !is.data.frame(enrich_data) || ncol(enrich_data) == 1
  logMsg(sprintf("Running in %s mode", if(isVectorMode) "vector" else "table"))
  
  # Process input data
  if (isVectorMode) {
    go_ids <- if (is.data.frame(enrich_data)) {
      if (enrich_go_col %in% names(enrich_data)) 
        enrich_data[[enrich_go_col]] 
      else 
        enrich_data[[1]]
    } else {
      enrich_data
    }
    
    # Create standardized data frame
    enrich_df <- data.frame(
      Annot = go_ids,
      stringsAsFactors = FALSE
    )
    
    # Perform joining operation
    logMsg("Joining data")
    joined <- merge(
      enrich_df, 
      merged_data,
      by.x = "Annot",
      by.y = merged_go_col,
      all = (join_method == "left")
    )
    
  } else {
    # Validate required columns
    required_cols <- c(enrich_go_col, enrich_gene_col)
    missing_cols <- setdiff(required_cols, names(enrich_data))
    if (length(missing_cols) > 0) {
      stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }
    
    # Process table input
    logMsg("Processing enrichment table")
    
    # Split gene lists
    genes_list <- fastSplit(enrich_data[[enrich_gene_col]], gene_delim)
    enrich_sep <- data.frame(
      row = rep(seq_len(nrow(enrich_data)), lengths(genes_list)),
      gene = unlist(genes_list, use.names = FALSE),
      stringsAsFactors = FALSE
    )
    
    # Add other columns
    other_cols <- setdiff(names(enrich_data), enrich_gene_col)
    enrich_sep <- cbind(
      enrich_sep,
      enrich_data[enrich_sep$row, other_cols, drop = FALSE]
    )
    
    # Perform joining operation
    logMsg("Joining data")
    joined <- merge(
      enrich_sep,
      merged_data,
      by.x = c(enrich_go_col, "gene"),
      by.y = c(merged_go_col, merged_gene_col),
      all = (join_method == "left")
    )
  }
  
  # Filter unmatched if needed
  if (!include_unmatched) {
    joined <- joined[!is.na(joined[[merged_complex_col]]), ]
  }
  
  # Check for empty results
  if (nrow(joined) == 0) {
    warning("No matches found between enrichment data and complex annotations")
    return(list(
      summary = data.frame(),
      matchedAnnot = data.frame()
    ))
  }
  
  # Compute complex summaries
  logMsg("Computing complex summaries")
  complex_groups <- split(joined, joined[[merged_complex_col]])
  result_list <- lapply(complex_groups, function(df) {
    go_col <- if(isVectorMode) "Annot" else enrich_go_col
    
    # Initialize basic result
    result <- data.frame(
      ComplexName = unique(df[[merged_complex_col]]),
      GO = paste(unique(df[[go_col]]), collapse = ","),
      Genes = if(!isVectorMode) {
        paste(unique(df$gene), collapse = ",")
      } else {
        paste(unique(df[[merged_gene_col]]), collapse = ",")
      },
      Matches = length(unique(df[[go_col]])),
      stringsAsFactors = FALSE
    )
    
    # Add configured column aggregations
    for (col_name in names(column_config)) {
      config <- column_config[[col_name]]
      if (!is.null(config$col) && config$col %in% names(df)) {
        result[[col_name]] <- aggregateColumns(df, config$col, config$agg)
      }
    }
    
    result
  })
  result <- do.call(rbind, result_list)
  
  # Compute background statistics
  logMsg("Computing background statistics")
  bg_stats <- tapply(
    merged_data[[merged_go_col]],
    merged_data[[merged_complex_col]],
    function(x) length(unique(x))
  )
  bg <- data.frame(
    ComplexName = names(bg_stats),
    Background = as.numeric(bg_stats),
    stringsAsFactors = FALSE
  )
  
  # Merge results with background
  result <- merge(result, bg, by = "ComplexName", all.x = TRUE)
  result$Ratio <- result$Matches / result$Background
  
  # Compute Fisher's exact test statistics
  logMsg("Computing Fisher's exact test statistics")
  input_total <- if(isVectorMode) {
    length(unique(enrich_df[["Annot"]]))
  } else {
    length(unique(enrich_data[[enrich_go_col]]))
  }
  total_background <- length(unique(merged_data[[merged_go_col]]))
  
  result$P_Com <- mapply(function(matches, background) {
    tryCatch({
      contingency <- matrix(c(
        matches,                           # a: matches
        input_total - matches,             # b: input not matched
        background - matches,              # c: background matched but not in input
        total_background - background - (input_total - matches)  # d: remainder
      ), nrow = 2)
      
      stats::fisher.test(contingency, alternative = fisher_alternative)$p.value
    }, error = function(e) {
      warning("Fisher test failed for a complex: ", e$message)
      NA
    })
  }, result$Matches, result$Background)
  
  # Adjust p-values if requested
  if (adjust_p) {
    result$adj_P_Com <- stats::p.adjust(result$P_Com, method = p_adjust_method)
  }
  
  # Add metadata if available
  logMsg("Adding metadata")
  meta_cols <- c("Functional_Complex_Group", "Root", "FCG_assoc_GO_ID", "ComplexID")
  available_meta <- intersect(meta_cols, names(merged_data))
  
  if (length(available_meta) > 0) {
    meta_data <- unique(merged_data[c(merged_complex_col, available_meta)])
    result <- merge(result, meta_data, by = merged_complex_col, all.x = TRUE)
  }
  
  # Sort results by p-value
  result <- result[order(result$P_Com), ]
  rownames(result)<-NULL
  
  logMsg("Analysis complete")
  return(list(
    summary = result,
    matchedAnnot = joined
  ))
}
