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
tree <- buildOntologyTree(namespace = "BP")
```

### 3. Compute Term Similarity

```r
sim_matrix <- clusterST(tree, terms = c("GO:0008150", "GO:0003674"))
```

### 4. Compute KEGG Pathway Similarity

```r
sim_kegg <- clusterST(tree, terms = c("hsa04110", "hsa04060"))
```

### 5. Cluster Similar Terms with Custom Method and Parameters

```r
clu <- clusterST(tree, 1:30, 
                 method = "wang", 
                 weights = weights, 
                 threshold = 0.2)
```

### 6. Visualize the Network

```r
plotClusterNetwork(cluster_result = clu)
```

### 7. Visualize Heatmap of Similarity Clusters

```r
plotClusterHeatmap(cluster_result = clu)
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



