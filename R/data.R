#' Built-in Gene Ontology Tree Dataset
#'
#' This built-in dataset contains Gene Ontology Trees for several species, represented as a list of
#' four \code{OntologyTree} objects (class \code{"OntologyTree"} from the \code{simAnn} package). The dataset includes:
#'
#' \describe{
#'   \item{human}{An \code{OntologyTree} object containing the GO hierarchy for human (species code \code{"hsa"}).}
#'   \item{mouse}{An \code{OntologyTree} object containing the GO hierarchy for mouse (species code \code{"mmu"}).}
#'   \item{rat}{An \code{OntologyTree} object containing the GO hierarchy for rat (species code \code{"rno"}).}
#' }
#'
#' Each \code{OntologyTree} object in this list has 15 slots that store various aspects of the hierarchical
#' ontology (e.g., term names, parent/child mappings, edge relations, annotations, metadata, and term-level statistics).
#'
#' @format A list with 4 elements. Each element is an \code{OntologyTree} object with the following slots:
#' \describe{
#'   \item{termNames}{A character vector of term names.}
#'   \item{termCount}{An integer indicating the total number of terms.}
#'   \item{edgeCount}{An integer indicating the total number of relations (edges).}
#'   \item{parentMap}{A list of parent term indices.}
#'   \item{childMap}{A list of child term indices.}
#'   \item{edgeTypes}{A list of semantic relation types between terms.}
#'   \item{relTree}{An ontology tree object for relation types.}
#'   \item{source}{A character scalar indicating the data source.}
#'   \item{rootNode}{An integer scalar indicating the index of the root term.}
#'   \item{leafNodes}{A vector of indices corresponding to the leaf terms.}
#'   \item{annotations}{A list containing term annotation information.}
#'   \item{n_parents}{A vector containing the number of parents for each term.}
#'   \item{n_children}{A vector containing the number of children for each term.}
#'   \item{termStats}{A list for storing arbitrary term-level statistics.}
#'   \item{metadata}{A data frame containing additional term metadata.}
#' }
#'
#' @source GO hierarchy database.
#'
#' @examples
#' data("goTree")
#' names(goTree)
#'
"goTree"

#' Built-in KEGG Ontology Tree Dataset
#'
#' This built-in dataset contains KEGG Ontology Trees for several species, represented as a list of
#' four \code{OntologyTree} objects (class \code{"OntologyTree"} from the \code{simAnn} package). The dataset includes:
#'
#' \describe{
#'   \item{human}{An \code{OntologyTree} object containing the KEGG BRITE hierarchy for human (species code \code{"hsa"}).}
#'   \item{mouse}{An \code{OntologyTree} object containing the KEGG BRITE hierarchy for mouse (species code \code{"mmu"}).}
#'   \item{rat}{An \code{OntologyTree} object containing the KEGG BRITE hierarchy for rat (species code \code{"rno"}).}
#' }
#'
#' Each \code{OntologyTree} object in this list has 15 slots that store various aspects of the hierarchical
#' ontology (e.g., term names, parent/child mappings, edge relations, annotations, metadata, and term-level statistics).
#'
#' @format A list with 4 elements. Each element is an \code{OntologyTree} object with the following slots:
#' \describe{
#'   \item{termNames}{A character vector of term names.}
#'   \item{termCount}{An integer indicating the total number of terms.}
#'   \item{edgeCount}{An integer indicating the total number of relations (edges).}
#'   \item{parentMap}{A list of parent term indices.}
#'   \item{childMap}{A list of child term indices.}
#'   \item{edgeTypes}{A list of semantic relation types between terms.}
#'   \item{relTree}{An ontology tree object for relation types.}
#'   \item{source}{A character scalar indicating the data source.}
#'   \item{rootNode}{An integer scalar indicating the index of the root term.}
#'   \item{leafNodes}{A vector of indices corresponding to the leaf terms.}
#'   \item{annotations}{A list containing term annotation information.}
#'   \item{n_parents}{A vector containing the number of parents for each term.}
#'   \item{n_children}{A vector containing the number of children for each term.}
#'   \item{termStats}{A list for storing arbitrary term-level statistics.}
#'   \item{metadata}{A data frame containing additional term metadata.}
#' }
#'
#' @source KEGG BRITE database.
#'
#' @examples
#' data("keggTree")
#' names(keggTree)
#'
"keggTree"


#' KEGG Leaf Hierarchy Data
#'
#' A data frame containing leaf-level hierarchical information extracted from the KEGG BRITE database.
#'
#' This dataset provides the hierarchical classification of KEGG BRITE terms at the leaf level,
#' including the codes and annotations for three levels:
#' \describe{
#'   \item{Level1}{A character vector of Level1 KEGG codes.}
#'   \item{Level1Anno}{A character vector of Level1 annotations.}
#'   \item{Level2}{A character vector of Level2 KEGG codes.}
#'   \item{Level2Anno}{A character vector of Level2 annotations.}
#'   \item{Level3}{A character vector of Level3 KEGG codes.}
#'   \item{Level3Anno}{A character vector of Level3 annotations.}
#' }
#'
#' @format A data frame with 363 observations of 6 variables.
#'
#' @source KEGG BRITE database.
#'
#' @examples
#' data("leaf_df")
#' str(leaf_df, 2)
"leaf_df"

#' Mouse Protein Complex Dataset
#'
#' This dataset contains information about protein complexes in mouse, including their functional groups,
#' associated Gene Ontology (GO) terms, and component proteins. Each complex is described with multiple
#' attributes including complex name, functional grouping, GO annotations, and subunit information.
#'
#' @format A data frame with multiple rows (one per subunit-complex combination) and 11 columns:
#' \describe{
#'   \item{ComplexName}{Character. Name of the protein complex.}
#'   \item{Functional_Complex_Group}{Character. The functional group classification of the complex.}
#'   \item{Root}{Character. Root classification identifier.}
#'   \item{FCG_assoc_GO_ID}{Character. GO identifiers associated with the functional complex group.}
#'   \item{ComplexID}{Character. Unique identifier for the complex.}
#'   \item{subunits_gene_name}{Character. Gene names of complex subunits.}
#'   \item{organism}{Character. Source organism (Mouse).}
#'   \item{functions_go_id}{Character. GO term identifiers describing complex functions.}
#'   \item{functions_go_name}{Character. Names of the GO terms.}
#'   \item{functions_go_ontology}{Character. GO domain (e.g., cellular_component).}
#'   \item{Ontology_type}{Character. Type of ontology (e.g., CC for Cellular Component).}
#' }
#'
#' This dataset provides a comprehensive view of protein complexes and their functional annotations,
#' particularly useful for studying protein complex organization and function in mouse.
#'
#' @source Compiled from protein complex databases and GO annotations.
#'
#' @examples
#' data("mComplex")
#' head(mComplex)
#' 
#' # Get unique complex names
#' unique(mComplex$ComplexName)
#'
#' # Filter complexes by GO term
#' er_complexes <- mComplex[mComplex$functions_go_id == "GO:0005783", ]
#'
"mComplex"

#' Human Protein Complex Dataset
#'
#' This dataset contains information about protein complexes in human, including their functional
#' associations, component proteins, and Gene Ontology (GO) annotations. Each row represents
#' a subunit-complex-function combination, providing detailed molecular function annotations
#' for protein complexes.
#'
#' @format A data frame with multiple rows (one per subunit-complex-function combination) and 10 columns:
#' \describe{
#'   \item{ComplexName}{Character. Name of the protein complex.}
#'   \item{Functional_Complex_Group}{Character. The functional group classification of the complex.}
#'   \item{Root}{Character. Root classification identifier.}
#'   \item{FCG_assoc_GO_ID}{Character. GO identifiers associated with the functional complex group.}
#'   \item{ComplexID}{Character. Unique identifier for the complex.}
#'   \item{subunits_gene_name}{Character. Gene names of complex subunits.}
#'   \item{functions_go_id}{Character. GO term identifiers describing complex functions.}
#'   \item{functions_go_name}{Character. Names of the GO terms.}
#'   \item{functions_go_ontology}{Character. GO domain (e.g., molecular_function).}
#'   \item{Ontology_type}{Character. Type of ontology (e.g., MF for Molecular Function).}
#' }
#'
#' This dataset is particularly useful for studying protein complex functions and their molecular
#' activities in human cells. It can be used for analyzing protein-protein interactions,
#' functional annotations, and complex organization.
#'
#' @source Compiled from protein complex databases and GO annotations.
#'
#' @examples
#' data("hComplex")
#' head(hComplex)
#' 
#' # Get unique complexes
#' unique(hComplex$ComplexName)
#'
#' # Filter by molecular function
#' kinase_binding <- hComplex[hComplex$functions_go_id == "GO:0051018", ]
#'
#' # Count unique complexes with specific function
#' length(unique(hComplex[hComplex$functions_go_name == "protein kinase A binding", "ComplexName"]))
#'
"hComplex"
