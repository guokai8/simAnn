# simAnn: Semantic Similarity Analysis in R

## Overview

**simAnn** is an R package designed for semantic similarity analysis in biological ontologies, including Gene Ontology (GO) and KEGG pathways. It provides methods for ontology traversal, clustering, and similarity calculations using various approaches, including information content and graph-based methods. The package facilitates ontology-based clustering and visualization, supporting the analysis of hierarchical relationships among biological terms and genes.

## Features

- **Ontology Tree Construction**: Builds hierarchical relationships between terms based on ontological data.
- **Semantic Similarity Calculation**: Implements term-to-term and group-wise similarity calculations.
- **Gene Ontology (GO) Analysis**: Processes gene annotations and extracts meaningful relationships.
- **KEGG Pathway Analysis**: Supports pathway-based semantic similarity calculations.
- **Flexible Clustering Methods**: Groups terms based on similarity scores and hierarchical structures.
- **Custom Thresholding and Weighting**: Allows fine-tuned similarity computations with adjustable parameters.
- **Network and Heatmap Visualization**: Generates term similarity networks and heatmaps.
- **Support for Multiple Ontologies**: Works with GO, KEGG, and other biological ontology systems.

## Installation

```r
# Install from GitHub
devtools::install_github("guokai8/simAnn")
```

## Usage

### 1. Load the Package

```r
library(simAnn)
```

### 2. Build Ontology Tree

```r
tree<-buildGOTree(species="human",namespace = "BP",keytype = "SYMBOL")
# or try some examples
parents  = c("a", "a", "b", "b", "c", "d")
children = c("b", "c", "c", "d", "e", "f")
annotation = list(
    "a" = 1:3,
    "b" = 3:4,
    "c" = 5,
    "d" = 7,
    "e" = 4:7,
    "f" = 8
)
tre<-buildOntologyTree(parentTerms = parents,childTerms = children,annotations = annotation)
trek <- buildKEGGTree(species = "human", keytype = "SYMBOL")
```

### 3. Compute Term Similarity

```r
sim_mat <- simterm(tree, terms = c("GO:0008150", "GO:0003674"))
# or
sim_mat <- simterm(tre, terms = 1:6, method ="lin")

```

### 4. Compute Similarity Get Clusters

```r
##choose the thresholds
sim_mat<-simterm(tre,1:6)
threshold <- pickThreshold(sim_mat,thresholds = seq(0.1,0.9,0.01),cluster_method = "components","modularity")
clu <- clusterST(tre, terms = 1:6, threshold = threshold$best_threshold)
```

### 5. Cluster Similar Terms with Custom Method and Parameters

```r
weights <- c("is_a" = 0.8, "part_of" = 0.6)
sim_mat<-simterm(tre,1:50, method = "wang",weights = weights)
threshold <- pickThreshold(sim_mat,thresholds = seq(0.5,0.9,0.01),cluster_method = "components","modularity")
clu <- clusterST(tree, 1:50, 
                 method = "wang", 
                 weights = weights, 
                 threshold$best_threshold)
```
### 6. Cluster Similar Terms with Custom Gene Sets (a list with term name and gene id)

```r
clu <- clusterSTW(tree, names(geneset), geneset,
                 method = "wang", 
                 weights = weights, 
                 threshold = 0.2)
```

### 7. Visualize the Network

```r
plotClusterNetwork(cluster_result = clu)
```

### 8. Visualize Heatmap of Similarity Clusters

```r
plotClusterHeatmap(cluster_result = clu)
```
### 8. All Available methods

```r
AllSimMethod()
AllICMethod()
```

## Documentation

To see a full list of functions and their descriptions, use:

```r
help(package = "simAnn")
```

## Author

- **Kai Guo** ([guokai8@gmail.com](mailto:guokai8@gmail.com))

## License

This package is released under the GPL_v3 License.

## Citation

If you use **simAnn** in your research, please contact guokai8@gmail.com

## Acknowledge
This project includes some modified code from _simona_ package. 



