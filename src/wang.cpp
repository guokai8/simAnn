#include <Rcpp.h>
using namespace Rcpp;
inline void add_parents(const List& parentMap, const int node, LogicalVector& ancestors) {
  const IntegerVector& parents = parentMap[node];
  for(int i = 0; i < parents.size(); i++) {
    const int parent = parents[i] - 1;
    if(!ancestors[parent]) {
      ancestors[parent] = true;
      add_parents(parentMap, parent, ancestors);
    }
  }
}

inline double wangs(const List& childMap, const List& edgeTypes, 
                    const NumericVector& contribution, const int node, const int end, 
                    const LogicalVector& background) {
  if(node == end) return 1.0;
  
  const IntegerVector& children = childMap[node];
  const IntegerVector& relations = edgeTypes[node];
  double maxScore = 0.0;
  
  for(int i = 0; i < children.size(); i++) {
    const int child = children[i] - 1;
    if(background[child]) {
      const double score = wangs(childMap, edgeTypes, contribution, child, end, background) * 
        contribution[relations[i] - 1];
      if(score > maxScore) maxScore = score;
    }
  }
  return maxScore;
}

// [[Rcpp::export]]
NumericVector calculateWangIC(S4 tree, NumericVector contribution) {
  const List parentMap = tree.slot("parentMap");
  const List childMap = tree.slot("childMap");
  const List edgeTypes = tree.slot("edgeTypes");
  const int n = parentMap.size();
  
  NumericVector ic(n);
  LogicalVector ancestors(n);
  
  for(int i = 0; i < n; i++) {
    std::fill(ancestors.begin(), ancestors.end(), false);
    ancestors[i] = true;
    add_parents(parentMap, i, ancestors);
    
    double sum = 0.0;
    for(int j = 0; j < n; j++) {
      if(ancestors[j]) {
        sum += wangs(childMap, edgeTypes, contribution, j, i, ancestors);
      }
    }
    ic[i] = sum;
  }
  
  return ic;
}