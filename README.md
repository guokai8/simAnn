# simAnn: Semantic Similarity Analysis in R

## Overview

**simAnn** is an R package designed for semantic similarity analysis in biological ontologies, particularly Gene Ontology (GO) and KEGG pathways. It provides comprehensive methods for ontology traversal, clustering, and similarity calculations using various approaches, including information content (IC) and graph-based methods. The package facilitates ontology-based clustering and visualization, supporting both standard ontological analyses and protein complex-based clustering.

## Features

- **Multiple Ontology Support**:
  - Gene Ontology (GO) analysis with species-specific annotations
  - KEGG pathway analysis with hierarchical relationships
  - Custom ontology tree construction
  - Support for different GO namespaces (BP, CC, MF)
  
- **Semantic Similarity Analysis**:
  - Term-based methods (kappa, overlap, jaccard, dice)
  - Information content methods (universal, annotation-based)
  - Graph-based methods (Wang's method, GOGO)
  - Customizable edge weights for relationship types
  - Support for new annotations and enrichment results

- **Clustering Capabilities**:
  - Standard term clustering based on similarity scores
  - Protein complex-based clustering for GO terms
  - Customizable threshold selection
  - Multiple clustering methods and parameters
  - Support for weighted clustering with gene sets

- **Visualization Tools**:
  - Network visualization of term relationships
  - Hierarchical and circular tree layouts
  - Customizable heatmaps with multiple color schemes
  - Interactive term focusing and highlighting
  - Various annotation styles and color schemes

## Installation

```r
# Install from GitHub
devtools::install_github("guokai8/simAnn")
```

## Usage

### 1. Load Required Packages

```r
library(simAnn)
# For visualization
library(ggplot2)
library(igraph)
```

### 2. Build Ontology Trees

```r
# Build GO tree with specific namespace and species
tree <- buildGOTree(
    species = "human", 
    namespace = "BP",    # Options: "BP", "CC", "MF"
    keytype = "SYMBOL",  # Options: "SYMBOL", "ENTREZID"
    builtin = TRUE,      # Use built-in data
    verbose = TRUE       # Show progress messages
)

# Build KEGG tree
kegg_tree <- buildKEGGTree(
    species = "human", 
    keytype = "SYMBOL",
    builtin = TRUE
)

# Build custom ontology tree
parents <- c("a", "a", "b", "b", "b", "c", "d")
children <- c("b", "c", "c", "d", "e", "e", "f")
annotation <- list(
    "a" = 1:3,
    "b" = 3:4,
    "c" = 5,
    "d" = 7,
    "e" = 4:7,
    "f" = 8
)
custom_tree <- buildOntologyTree(
    parentTerms = parents,
    childTerms = children,
    annotations = annotation
)
```

### 3. Calculate Semantic Similarity

#### Basic Similarity Calculation
```r
# Different similarity methods
sim_mat <- simterm(tree, 
                     terms = c("GO:0008150", "GO:0000011"),
                     method = "term")

sim_mat_resnik <- simterm(tree,
                         terms = 1:30,
                         method = "resnik",
                         normType = "max")
```

### Cluster Similar Terms with Custom Method and Parameters

```r
weights <- c("is_a" = 0.8, "part_of" = 0.6)
sim_mat<-simterm(tree,1:50, method = "wang",weights = weights)
threshold <- pickThreshold(sim_mat,thresholds = seq(0.5,0.9,0.01),cluster_method = "components",quality_metric = "modularity")
clu <- clusterST(tree, 1:50, 
                 method = "wang", 
                 weights = weights, 
                 threshold$best_threshold)
```

#### Weighted Similarity with Gene Sets
```r

# Define example gene sets
geneset <- list(
    "GO:0000086" = c("TAF2", "FOXM1", "RRM1", "ATR"),
    "GO:0000422" = c("LRBA", "EIF2AK1", "SPATA18"),
    "GO:0001678" = c("SYBU", "KCNB1", "PLA2G6", "LIN28A"),
    "GO:0001704" = c("SNAI1", "GPI", "AHDC1"),
    "GO:0001707" = c("SNAI1", "GPI", "AHDC1"),
    "GO:0002573" = c("IL25", "CAMK4", "EIF2AK1", "TNFRSF11B", "IL17A"),
    "GO:0002761" = c("CAMK4", "TNFRSF11B", "IL17A"),
    "GO:0006282" = c("TAF2", "FOXM1", "BRD8", "ATR", "PARPBP"),
    "GO:0006284" = c("POLL", "APEX2"),
    "GO:0006302" = c("FOXM1", "BRD8", "POLL", "ATR", "PARPBP"),
    "GO:0006909" = c("IL2RB", "EIF2AK1", "PLA2G6", "TUB"),
    "GO:0007498" = c("SNAI1", "GPI", "AHDC1")
)

# Define weights for different relationship types
weights <- c("is_a" = 0.8, "part_of" = 0.6)

# Calculate weighted similarity
sim_mat_weighted <- simtermW(tree,
                           terms = names(geneset),
                           new_annot = geneset,
                           method = "wang",
                           weights = weights)
```

### 4. Cluster Similar Terms with custom gene sets

####  Clustering
```r
# Find optimal threshold
threshold <- pickThreshold(
    sim_mat_weighted,
    thresholds = seq(0.1, 0.9, 0.01),
    cluster_method = "components",
    quality_metric = "modularity"
)

# Cluster with gene sets
clu <- clusterSTW(tree, 
                  terms = names(geneset),
                  new_annot = geneset,
                  method = "wang", 
                  weights = weights, 
                  threshold = threshold$best_threshold)
```

#### Complex-based Clustering
```r
# Example enrichment data
enrichment_df <- data.frame(
    Annot = c("GO:0006915", "GO:0012501"),
    GeneID = c("CASP3,BCL2", "CASP3,BAX"),
    Pvalue = c(0.001, 0.002)
)

# Cluster using protein complex information
complex_clusters <- clusterComplex(
    enrichment_df,
    species = "human",
    ontology_type = "BP",
    gene_delim = ",",
    verbose = TRUE
)

# With enrichment results from richR
# clu <- clusterComplex(
#     as.data.frame(res),
#     species = "human",
#     ontology_type = "BP",
#     verbose = TRUE
# )
```

### 5. Visualize Results

#### Tree Visualization
```r
# Basic tree layout
plotOntologyTree(custom_tree, layout = "tree")

# Focus on specific terms
plotOntologyTree(tree, 
                focus_terms = names(geneset)[1:5],
                layout = "tree")
```

#### Network Visualization
```r
# Plot cluster network
plotClusterNetwork(clu,
                  show_names = TRUE,
                  annotation_style = "combined",
                  layout = "fr",
                  vertex_size = 12)
```

#### Heatmap Visualization
```r
# Plot similarity heatmap
plotClusterHeatmap(clu,
                  show_names = TRUE,
                  color_scheme = "brewer",
                  annotation_style = "separate")
```

### 6. Available Methods

```r
# List all similarity methods
AllSimMethod()

# List all IC calculation methods
AllICMethod()
```

## Documentation

For detailed function documentation and examples:

```r
help(package = "simAnn")
?simterm      # Basic similarity calculation
?simtermW     # Weighted similarity calculation
?clusterST    # Basic clustering
?clusterSTW   # Weighted clustering
?clusterComplex   # Complex-based clustering
```

## Author

- **Kai Guo** ([guokai8@gmail.com](mailto:guokai8@gmail.com))

## License

This package is released under the GPL_v3 License.

## Citation

If you use **simAnn** in your research, please contact guokai8@gmail.com

## Acknowledgments

This project includes some modified code from the _simona_ package.
