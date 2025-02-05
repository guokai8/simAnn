# simAnn: Semantic Similarity Analysis in R

## Overview

**simAnn** is an R package designed for semantic similarity analysis in biological ontologies, including Gene Ontology (GO) and KEGG pathways. It provides methods for ontology traversal, clustering, and similarity calculations using information content-based approaches. The package supports analyzing hierarchical relationships among biological terms and genes, facilitating ontology-based clustering and visualization.

## Features

- **Ontology Tree Construction**: Builds hierarchical relationships between terms based on ontological data.
- **Semantic Similarity Calculation**: Implements term-to-term and group-wise similarity calculations based on shared ancestors.
- **Gene Ontology (GO) Analysis**: Processes gene annotations and extracts meaningful relationships.
- **KEGG Pathway Analysis**: Supports pathway-based semantic similarity calculations.
- **Clustering of Similar Terms**: Groups terms based on similarity scores and hierarchical structures.
- **Network Visualization**: Generates parent-child networks and term similarity networks.
- **Support for Multiple Ontologies**: Works with GO, KEGG, and other biological ontology systems.

## Installation

```r
# Install from source
devtools::install_github("guokai8/simAnn")
```

## Usage

### 1. Load the Package

```r
library(simAnn)
```

### 2. Build Ontology Tree

```r
tree <- buildGOTree(namespace = "BP",org)
```

### 3. Compute Term Similarity

```r
sim_matrix <- clusterST(tree, terms = c("GO:0008150", "GO:0003674"))
```

### 4. Cluster Similar Terms

```r
clusters <- clusterST(terms = c("GO:0008150", "GO:0003674"), method = "ic")
```

### 5. Visualize the Network

```r
plotClusterNetwork(cluster_result = clusters)
```

## Author

- **Kai Guo** ([kaiguo@guokai8.com](mailto:guokai8@gmail.com))

## License

This package is released under the MIT License.

## Citation

If you use **simAnn** in your research, please contact guokai8@gmail.com

