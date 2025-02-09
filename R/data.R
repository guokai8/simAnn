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
