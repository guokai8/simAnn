% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/buildKEGG.R
\name{getGenePathwayMapping}
\alias{getGenePathwayMapping}
\title{Get Gene–Pathway Mapping from KEGG}
\usage{
getGenePathwayMapping(species = "human", keytype = "ENTREZID")
}
\arguments{
\item{species}{A string specifying the species code (e.g., "human").
Default is "hsa".}

\item{keytype}{A character string specifying the target key type (e.g., "SYMBOL").
Default is "SYMBOL".}
}
\value{
A data frame with two columns: GeneID and Level3
}
\description{
This function retrieves pathway information from KEGG by calling
keggLink("pathway", species) and returns a data frame with two columns:
  - GeneID: The gene (or enzyme) identifier.
  - PathwayID: The associated KEGG pathway ID.
}
\details{
Optionally, the gene IDs can be converted using the internal helper function .idc.
}
\examples{
\dontrun{
  # For human, converting KEGG gene IDs to SYMBOL:
  mapping <- getGenePathwayMapping("human", keytype =  "SYMBOL")
  head(mapping)
}

}
