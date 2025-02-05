#include <Rcpp.h>
using namespace Rcpp;

// Constants for traversal direction
const int TRAVERSE_TO_PARENTS = -1;
const int TRAVERSE_TO_CHILDREN = 1;

// Constants for distance calculation method
const bool USE_LONGEST_PATH = true;
const bool USE_SHORTEST_PATH = false;

/**
 * Perform breadth-first traversal on the ontology tree
 * 
 * @param tree Ontology tree object
 * @param startNodes Vector of starting node indices
 * @param useLongestPath Whether to use longest path for distance calculation
 * @param backgroundNodes Background nodes to consider (optional)
 * @param direction Traversal direction (TO_CHILDREN or TO_PARENTS)
 * @return Vector of distances for all nodes
 */
IntegerVector traverseTreeBFS(
    S4 tree, 
    IntegerVector startNodes = IntegerVector(0), 
    bool useLongestPath = USE_LONGEST_PATH,
    LogicalVector backgroundNodes = LogicalVector(0), 
    int direction = TRAVERSE_TO_CHILDREN
) {
  // Determine which node map to use based on direction
  String nodeMapName = (direction == TRAVERSE_TO_CHILDREN) ? "childMap" : "parentMap";
  List nodeMap = tree.slot(nodeMapName);
  
  // If no start nodes provided, use default nodes based on direction
  if(startNodes.size() == 0) {
    startNodes = (direction == TRAVERSE_TO_CHILDREN) ? 
    tree.slot("rootNode") : tree.slot("leafNodes");
  }
  
  // Convert to 0-based indices
  IntegerVector startIndices = startNodes - 1;
  
  // Initialize distance vector and visited flags
  int nodeCount = nodeMap.size();
  IntegerVector distances(nodeCount, -1);
  bool hasBackground = (backgroundNodes.size() > 0);
  
  // Initialize current level nodes
  LogicalVector currentLevelNodes(nodeCount, false);
  for(int i = 0; i < startIndices.size(); i++) {
    int nodeIndex = startIndices[i];
    currentLevelNodes[nodeIndex] = true;
    distances[nodeIndex] = 0;
  }
  
  // Perform BFS traversal
  int nodesInCurrentLevel = sum(currentLevelNodes);
  while(nodesInCurrentLevel > 0) {
    for(int i = 0; i < nodeCount; i++) {
      if(!currentLevelNodes[i]) continue;
      
      // Process current node
      currentLevelNodes[i] = false;
      int currentDistance = distances[i];
      
      // Process neighbors
      IntegerVector neighbors = nodeMap[i];
      for(int j = 0; j < neighbors.size(); j++) {
        int neighborIndex = neighbors[j] - 1;
        
        // Skip if not in background (when background is specified)
        if(hasBackground && !backgroundNodes[neighborIndex]) continue;
        
        int newDistance = currentDistance + 1;
        
        // Update distance based on strategy
        if(distances[neighborIndex] == -1) {
          // First visit to this node
          distances[neighborIndex] = newDistance;
        } else if(useLongestPath) {
          // Update if new path is longer
          distances[neighborIndex] = std::max(distances[neighborIndex], newDistance);
        } else {
          // Update if new path is shorter
          distances[neighborIndex] = std::min(distances[neighborIndex], newDistance);
        }
        
        currentLevelNodes[neighborIndex] = true;
      }
    }
    nodesInCurrentLevel = sum(currentLevelNodes);
  }
  
  return distances;
}

/**
 * Calculate depths of all nodes in the tree (distance from root)
 * 
 * @param tree Ontology tree object
 * @return Vector of node depths
 */
// [[Rcpp::export]]
IntegerVector calculateNodeDepths(S4 tree) {
  IntegerVector rootNode(1);
  rootNode[0] = tree.slot("rootNode");
  return traverseTreeBFS(
    tree, 
    rootNode, 
    USE_LONGEST_PATH, 
    LogicalVector(0), 
    TRAVERSE_TO_CHILDREN
  );
}


// [[Rcpp::export]]
NumericVector calculateUniversalIC(S4 tree, bool verbose = true) {
  if(verbose) {
    Rcout << "method: Universal\n";
  }
  
  // Get basic tree information
  List parentMap = tree.slot("parentMap");
  IntegerVector childCounts = tree.slot("n_children");
  int termCount = as<int>(tree.slot("termCount"));
  int rootNode = as<int>(tree.slot("rootNode")) - 1;  // Convert to 0-based index
  
  // Get node depths
  IntegerVector depths = calculateNodeDepths(tree);
  
  // Initialize IC values
  NumericVector ic(termCount, NA_REAL);
  ic[rootNode] = 0;
  
  // Process nodes level by level
  int currentDepth = 0;
  while(true) {
    currentDepth++;
    
    // Find nodes at current depth
    LogicalVector currentLevelNodes(termCount);
    int nodesAtCurrentDepth = 0;
    for(int i = 0; i < termCount; i++) {
      if(depths[i] == currentDepth) {
        currentLevelNodes[i] = true;
        nodesAtCurrentDepth++;
      }
    }
    
    // Break if no nodes at current depth
    if(nodesAtCurrentDepth == 0) break;
    
    // Calculate IC for nodes at current depth
    for(int i = 0; i < termCount; i++) {
      if(!currentLevelNodes[i]) continue;
      
      // Get parents of current node
      IntegerVector parents = parentMap[i];
      double sumIC = 0;
      
      // Sum IC values from parents
      for(int j = 0; j < parents.size(); j++) {
        int parentIndex = parents[j] - 1;  // Convert to 0-based index
        sumIC += ic[parentIndex] + log(childCounts[parentIndex]);
      }
      
      ic[i] = sumIC;
    }
  }
  
  return ic;
}


// _find_ancestors
inline void traverse_ancestors(const List& parentMap, int node, LogicalVector& ancestors) {
  IntegerVector parents = parentMap[node];
  for(int i = 0; i < parents.size(); i++) {
    int parent = parents[i] - 1;
    if(!ancestors[parent]) {
      ancestors[parent] = true;
      traverse_ancestors(parentMap, parent, ancestors);
    }
  }
}


inline void traverse_offspring(const List& childMap, const int node, LogicalVector& offspring, std::vector<int>& changedIndices) {
  const IntegerVector& children = childMap[node];
  const int n = children.size();
  for(int i = 0; i < n; i++) {
    const int child = children[i] - 1;
    if(!offspring[child]) {
      offspring[child] = true;
      changedIndices.push_back(child);
      traverse_offspring(childMap, child, offspring, changedIndices);
    }
  }
}

// [[Rcpp::export]]
IntegerMatrix getTermAnnotations(S4 otr, IntegerVector nodes) {
    List childMap = otr.slot("childMap");
    List annotation = otr.slot("annotations");
    List lt_annotation = annotation["list"];
    CharacterVector anno_names = annotation["names"];
    
    const int n = childMap.size();
    const int m = nodes.size();
    const int n_anno = anno_names.size();
    
    IntegerMatrix mat(m, n_anno);
    LogicalVector offspring(n);
    std::vector<int> changedIndices;
    changedIndices.reserve(n);
    
    for(int i = 0; i < m; i++) {
        // Reset for new term
        changedIndices.clear();
        offspring[nodes[i]-1] = true;
        changedIndices.push_back(nodes[i]-1);
        
        // Get all offspring
        traverse_offspring(childMap, nodes[i]-1, offspring, changedIndices);
        
        // Process annotations for all marked nodes
        for(size_t j = 0; j < changedIndices.size(); j++) {
            IntegerVector anno = lt_annotation[changedIndices[j]];
            for(int k = 0; k < anno.size(); k++) {
                mat(i, anno[k]-1) = 1;
            }
        }
        
        // Reset marked nodes
        for(size_t j = 0; j < changedIndices.size(); j++) {
            offspring[changedIndices[j]] = false;
        }
    }
    
    return mat;
}


// [[Rcpp::export]]
NumericMatrix cross_min(const NumericVector& x) {
  int n = x.length();
  NumericMatrix res(n, n);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      res(i,j) = std::min(x[i], x[j]);
    }
  }
  return res;
}


// [[Rcpp::export]]
IntegerVector countOffspring(S4 tree, bool includeSelf = false) {
  const List& childMap = tree.slot("childMap");
  const int n = childMap.size();
  IntegerVector num(n);
  LogicalVector offspring(n);
  std::vector<int> changedIndices;
  changedIndices.reserve(n);
  
  for(int i = 0; i < n; i++) {
    if(changedIndices.capacity() < n) {
      changedIndices.reserve(n);
    }
    changedIndices.clear();
    
    // 
    traverse_offspring(childMap, i, offspring, changedIndices);
    
    // 
    num[i] = changedIndices.size();
    if(includeSelf) num[i] += 1;
    
    // 
    for(size_t j = 0; j < changedIndices.size(); j++) {
      offspring[changedIndices[j]] = false;
    }
  }
  
  return num;
}


// [[Rcpp::export]]
NumericVector calculateOffspringIC(S4 tree, bool verbose = true) {
  if(verbose) {
    Rcout << "method: Offspring\n";
  }
  
  // Get offspring counts including self
  IntegerVector offspringCounts = countOffspring(tree, true);
  
  // Calculate IC values
  double maxCount = max(offspringCounts);
  int n = offspringCounts.length();
  NumericVector ic(n);
  
  for(int i = 0; i < n; i++) {
    double p = offspringCounts[i] / maxCount;
    ic[i] = -log(p);
  }
  
  return ic;
}

// [[Rcpp::export]]
NumericVector calculateIC(S4 tree, std::string method = "universal", bool useCache = true, bool verbose = true) {
  List termStats = tree.slot("termStats");
  std::string statName = method == "universal" ? "Universal" : "Offspring";
  
  // Try to get from cache
  if(useCache && termStats.containsElementNamed(statName.c_str())) {
    if(verbose) Rcout << "Using cached " << method << " IC values\n";
    return as<NumericVector>(termStats[statName]);
  }
  
  // Calculate IC values based on method
  NumericVector ic = method == "universal" ? 
  calculateUniversalIC(tree, verbose) : 
    calculateOffspringIC(tree, verbose);
  
  // Cache results
  if(useCache) {
    List updatedStats = clone(termStats);
    updatedStats[statName] = ic;
    tree.slot("termStats") = updatedStats;
  }
  
  return ic;
}

// [[Rcpp::export]]
NumericVector calculateUniversal(S4 tree, bool useCache = true, bool verbose = true) {
  return calculateIC(tree, "universal", useCache, verbose);
}

// [[Rcpp::export]]
NumericVector calculateOffspring(S4 tree, bool useCache = true, bool verbose = true) {
  return calculateIC(tree, "offspring", useCache, verbose);
}


// [[Rcpp::export]]
IntegerVector countAnnotations(S4 tree, bool unify = true) {
  List childMap = tree.slot("childMap");
  List annotations = tree.slot("annotations");
  List annotList = annotations["list"];
  CharacterVector annotNames = annotations["names"];
  
  int totalAnnots = annotNames.size();
  int termCount = childMap.size();
  
  IntegerVector annotCounts(termCount, 0);
  IntegerVector annotSizes(termCount);
  
  if(!unify) {
    for(int i = 0; i < termCount; i++) {
      IntegerVector currentAnnots = annotList[i];
      annotSizes[i] = currentAnnots.size();
    }
  }
  
  LogicalVector offspring(termCount, false);
  std::vector<int> changedIndices;
  changedIndices.reserve(termCount);
  
  for(int i = 0; i < termCount; i++) {
    changedIndices.clear();
    
    // Add self
    offspring[i] = true;
    changedIndices.push_back(i);
    
    traverse_offspring(childMap, i, offspring, changedIndices);
    
    if(unify) {
      LogicalVector uniqueAnnots(totalAnnots, false);
      
      // Use changedIndices instead of checking all nodes
      for(size_t j = 0; j < changedIndices.size(); j++) {
        IntegerVector currentAnnots = annotList[changedIndices[j]];
        for(int k = 0; k < currentAnnots.size(); k++) {
          uniqueAnnots[currentAnnots[k]-1] = true;
        }
      }
      annotCounts[i] = sum(uniqueAnnots);
      
    } else {
      int count = 0;
      for(size_t j = 0; j < changedIndices.size(); j++) {
        count += annotSizes[changedIndices[j]];
      }
      annotCounts[i] = count;
    }
    
    // Reset only changed nodes
    for(size_t j = 0; j < changedIndices.size(); j++) {
      offspring[changedIndices[j]] = false;
    }
  }
  
  return annotCounts;
}

// [[Rcpp::export]]
IntegerVector getLongestDistToOffspring(S4 tree, IntegerVector fromNodes, LogicalVector backgroundNodes = LogicalVector(0)) {
  return traverseTreeBFS(tree, fromNodes, USE_LONGEST_PATH, backgroundNodes, TRAVERSE_TO_CHILDREN);
}

// Overload for single node
IntegerVector getLongestDistToOffspring(S4 tree, int fromNode, LogicalVector backgroundNodes = LogicalVector(0)) {
  IntegerVector nodes(1);
  nodes[0] = fromNode;
  return traverseTreeBFS(tree, nodes, USE_LONGEST_PATH, backgroundNodes, TRAVERSE_TO_CHILDREN);
}

// [[Rcpp::export]]
IntegerVector getShortestDistToOffspring(S4 tree, IntegerVector fromNodes, LogicalVector backgroundNodes = LogicalVector(0)) {
  return traverseTreeBFS(tree, fromNodes, USE_SHORTEST_PATH, backgroundNodes, TRAVERSE_TO_CHILDREN);
}

// Overload for single node
IntegerVector getShortestDistToOffspring(S4 tree, int fromNode, LogicalVector backgroundNodes = LogicalVector(0)) {
  IntegerVector nodes(1);
  nodes[0] = fromNode;
  return traverseTreeBFS(tree, nodes, USE_SHORTEST_PATH, backgroundNodes, TRAVERSE_TO_CHILDREN);
}


#define SET_UNION 1
#define SET_INTERSECT 2

// [[Rcpp::export]]
IntegerVector getGroupAncestors(S4 tree, IntegerVector nodes, int type = 1, bool includeSelf = false) {
  List parentMap = tree.slot("parentMap");
  int n = parentMap.size();
  
  if(type == 1) {  // union
    LogicalVector ancestors(n, false);
    for(int i = 0; i < nodes.length(); i++) {
      int currentNode = nodes[i] - 1;  // Convert to 0-based index
      traverse_ancestors(parentMap, currentNode, ancestors);
      if(includeSelf) ancestors[currentNode] = true;
    }
    
    // Convert to indices
    int resultSize = sum(ancestors);
    IntegerVector result(resultSize);
    int idx = 0;
    for(int i = 0; i < n; i++) {
      if(ancestors[i]) {
        result[idx++] = i + 1;  // Convert back to 1-based index
      }
    }
    return result;
    
  } else {  // intersection
    LogicalVector temp(n, false);
    LogicalVector result(n, true);
    
    for(int i = 0; i < nodes.length(); i++) {
      int currentNode = nodes[i] - 1;
      traverse_ancestors(parentMap, currentNode, temp);
      if(includeSelf) temp[currentNode] = true;
      
      result = result & temp;
      std::fill(temp.begin(), temp.end(), false);
    }
    
    // Convert to indices
    int resultSize = sum(result);
    IntegerVector finalResult(resultSize);
    int idx = 0;
    for(int i = 0; i < n; i++) {
      if(result[i]) {
        finalResult[idx++] = i + 1;
      }
    }
    return finalResult;
  }
}


inline void addChildrenInBackground(const List& childMap, 
                                    int node,
                                    LogicalVector& offspring, 
                                    const LogicalVector& background) {
  // Only proceed if current node is in background
  if(background[node]) {
    IntegerVector children = childMap[node];
    
    // Process each child
    for(int i = 0; i < children.size(); i++) {
      int childNode = children[i] - 1;  // Convert to 0-based index
      
      // Add child if it's in background and not already processed
      if(background[childNode] && !offspring[childNode]) {
        offspring[childNode] = true;
        addChildrenInBackground(childMap, childNode, offspring, background);
      }
    }
  }
}

inline void getOffspringInBackground(const List& childMap,
                                     int node,
                                     LogicalVector& offspring,
                                     const LogicalVector& background,
                                     bool includeSelf = false) {
  // Mark self if required
  if(includeSelf) {
    offspring[node] = true;
  }
  
  // Get all descendants within background
  addChildrenInBackground(childMap, node, offspring, background);
}

// Convert integer vector to logical vector
inline LogicalVector intToLogical(const IntegerVector& indices, int size) {
  LogicalVector result(size, false);
  for(int i = 0; i < indices.length(); i++) {
    result[indices[i]] = true;
  }
  return result;
}

// [[Rcpp::export]]
IntegerMatrix findMaxAncestorID(S4 tree, IntegerVector nodes, NumericVector values, 
                                bool useLongestPath = true) {
  List childMap = tree.slot("childMap");
  int rootNode = tree.slot("rootNode");
  
  int termCount = childMap.size();
  int nodeCount = nodes.size();
  
  // Initialize result matrices
  NumericMatrix scores(nodeCount, nodeCount);  // Store values
  IntegerMatrix ancestorIDs(nodeCount, nodeCount);  // Store ancestor IDs
  IntegerMatrix distances(nodeCount, nodeCount);    // Store distances
  
  // Create node mapping
  IntegerVector nodeMapping(termCount, -1);
  for(int i = 0; i < nodeCount; i++) {
    nodeMapping[nodes[i]-1] = i;
    ancestorIDs(i, i) = nodes[i];
  }
  
  if(nodeCount <= 1) return ancestorIDs;
  
  // Get all ancestors
  IntegerVector allAncestors = getGroupAncestors(tree, nodes, 1, true);
  LogicalVector offspring(termCount);
  LogicalVector allAncestorsFlags = intToLogical(allAncestors - 1, termCount);
  
  // Process each ancestor
  for(int k = 0; k < allAncestors.size(); k++) {
    int currentAncestor = allAncestors[k] - 1;
    
    // Get offspring within ancestors
    if(allAncestors[k] == rootNode) {
      offspring = clone(allAncestorsFlags);
    } else {
      getOffspringInBackground(childMap, currentAncestor, offspring, 
                               allAncestorsFlags, true);
    }
    
    // Convert offspring to indices and process
    IntegerVector offspringNodes(sum(offspring));
    int idx = 0;
    for(int i = 0; i < termCount; i++) {
      if(offspring[i]) {
        offspringNodes[idx++] = i + 1;
      }
    }
    std::fill(offspring.begin(), offspring.end(), false);
    
    if(offspringNodes.size() <= 1) continue;
    
    // Calculate distances
    IntegerVector depths = useLongestPath ? 
    getLongestDistToOffspring(tree, allAncestors[k], allAncestorsFlags) :
      getShortestDistToOffspring(tree, allAncestors[k], allAncestorsFlags);
    
    // Update matrices
    for(int i = 0; i < offspringNodes.size() - 1; i++) {
      int id1 = nodeMapping[offspringNodes[i]-1];
      if(id1 >= 0) {
        for(int j = i+1; j < offspringNodes.size(); j++) {
          int id2 = nodeMapping[offspringNodes[j]-1];
          if(id2 >= 0) {
            double currentValue = values[currentAncestor];
            int dist = depths[offspringNodes[i]-1] + depths[offspringNodes[j]-1];
            
            // Update if better value found
            if(ancestorIDs(id1, id2) == 0 || currentValue > scores(id1, id2)) {
              scores(id1, id2) = currentValue;
              ancestorIDs(id1, id2) = currentAncestor + 1;
              ancestorIDs(id2, id1) = ancestorIDs(id1, id2);
              distances(id1, id2) = dist;
              distances(id2, id1) = dist;
            }
            // Update if same value but better distance
            else if(std::abs(currentValue - scores(id1, id2)) < 1e-10) {
              bool updateDistance = useLongestPath ? 
              (dist < distances(id1, id2)) : 
              (dist > distances(id1, id2));
              
              if(updateDistance) {
                scores(id1, id2) = currentValue;
                ancestorIDs(id1, id2) = currentAncestor + 1;
                ancestorIDs(id2, id1) = ancestorIDs(id1, id2);
                distances(id1, id2) = dist;
                distances(id2, id1) = dist;
              }
            }
          }
        }
      }
    }
  }
  
  return ancestorIDs;
}



// [[Rcpp::export]]
NumericMatrix findMaxAncestorValues(S4 tree, IntegerVector nodes, NumericVector values) {
  // Get basic tree information
  List childMap = tree.slot("childMap");
  int rootNode = tree.slot("rootNode");
  
  int termCount = childMap.size();
  int nodeCount = nodes.size();
  
  // Initialize result matrix and node mapping
  NumericMatrix scores(nodeCount, nodeCount);
  IntegerVector nodeMapping(termCount, -1);
  
  // Set diagonal values and create node mapping
  for(int i = 0; i < nodeCount; i++) {
    nodeMapping[nodes[i]-1] = i;
    scores(i, i) = values[nodes[i]-1];
  }
  
  if(nodeCount <= 1) return scores;
  
  Rcout << "Collecting all ancestors of input terms..." << std::endl;
  
  // Get all ancestors and convert to logical vector
  IntegerVector allAncestors = getGroupAncestors(tree, nodes, 1, true);
  LogicalVector offspring(termCount);
  LogicalVector allAncestorsFlags = intToLogical(allAncestors - 1, termCount);
  
  // Process each ancestor
  for(int k = 0; k < allAncestors.size(); k++) {
    int currentAncestor = allAncestors[k] - 1;
    
    // Get offspring within ancestors
    if(allAncestors[k] == rootNode) {
      offspring = clone(allAncestorsFlags);
    } else {
      getOffspringInBackground(childMap, currentAncestor, offspring, 
                               allAncestorsFlags, true);
    }
    
    // Convert offspring to indices
    IntegerVector offspringNodes(sum(offspring));
    int idx = 0;
    for(int i = 0; i < termCount; i++) {
      if(offspring[i]) {
        offspringNodes[idx++] = i + 1;
      }
    }
    std::fill(offspring.begin(), offspring.end(), false);
    
    if(offspringNodes.size() <= 1) continue;
    
    // Update scores for pairs of offspring
    double currentValue = values[currentAncestor];
    for(int i = 0; i < offspringNodes.size() - 1; i++) {
      int id1 = nodeMapping[offspringNodes[i]-1];
      if(id1 >= 0) {
        for(int j = i+1; j < offspringNodes.size(); j++) {
          int id2 = nodeMapping[offspringNodes[j]-1];
          if(id2 >= 0) {
            // Update max value if necessary
            if(scores(id1, id2) < currentValue) {
              scores(id1, id2) = currentValue;
              scores(id2, id1) = currentValue;
            }
          }
        }
      }
    }
  }
  
  return scores;
}

// [[Rcpp::export]]
NumericMatrix calculateAncestorSimilarity(S4 tree, IntegerVector termIds) {
  List childMap = tree.slot("childMap");
  const int termCount = childMap.size();
  const int querySize = termIds.size();
  
  // Initialize result matrices
  NumericMatrix similarity(querySize, querySize);
  NumericVector termOccurrence(querySize);
  NumericMatrix intersectionCounts(querySize, querySize);
  
  // Create term mapping
  IntegerVector termMapping(termCount, -1);
  for(int i = 0; i < querySize; i++) {
    termMapping[termIds[i] - 1] = i;
  }
  
  // Handle special cases
  if(querySize <= 1) {
    if(querySize == 1) {
      similarity(0, 0) = 1;
    }
    return similarity;
  }
  
  // Get ancestors and prepare background
  IntegerVector ancestorTerms = getGroupAncestors(tree, termIds, 1, true);
  LogicalVector offspring(termCount);
  LogicalVector ancestorBackground = intToLogical(ancestorTerms - 1, termCount);
  
  // Process each ancestor
  for(int k = 0; k < ancestorTerms.size(); k++) {
    getOffspringInBackground(childMap, ancestorTerms[k] - 1, offspring, 
                             ancestorBackground, true);
    
    // Get offspring terms
    IntegerVector currentOffspring(sum(offspring));
    int idx = 0;
    for(int i = 0; i < termCount; i++) {
      if(offspring[i]) {
        currentOffspring[idx++] = i + 1;
      }
    }
    std::fill(offspring.begin(), offspring.end(), false);
    
    int offspringCount = currentOffspring.size();
    if(offspringCount == 0) continue;
    
    // Count occurrences and intersections
    for(int i = 0; i < offspringCount; i++) {
      int id1 = termMapping[currentOffspring[i] - 1];
      if(id1 >= 0) {
        termOccurrence[id1]++;
        
        // Count intersections
        for(int j = i + 1; j < offspringCount; j++) {
          int id2 = termMapping[currentOffspring[j] - 1];
          if(id2 >= 0) {
            intersectionCounts(id1, id2)++;
            intersectionCounts(id2, id1) = intersectionCounts(id1, id2);
          }
        }
      }
    }
  }
  
  // Calculate final similarity scores
  for(int i = 0; i < querySize; i++) {
    similarity(i, i) = 1;
    for(int j = i + 1; j < querySize; j++) {
      double intersect = intersectionCounts(i, j);
      similarity(i, j) = intersect / 
        (termOccurrence[i] + termOccurrence[j] - intersect);
      similarity(j, i) = similarity(i, j);
    }
  }
  
  return similarity;
}