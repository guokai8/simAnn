#include <Rcpp.h>
using namespace Rcpp;

inline LogicalVector to_logical(IntegerVector i, int n) {
  LogicalVector l(n);
  for(int k = 0; k < i.size(); k++) l[i[k] - 1] = true;
  return l;
}

inline IntegerVector find_indices(LogicalVector l) {
  int n2 = sum(l);
  if(n2 == 0) return IntegerVector(0);
  IntegerVector ind(n2);
  int i2 = 0;
  for(int i = 0; i < l.size(); i++) 
    if(l[i]) ind[i2++] = i;
    return ind;
}

inline void traverse_otr(List list, int node, LogicalVector& result, const LogicalVector& background, bool is_parent) {
  IntegerVector nodes = list[node];
  if(nodes.size() > 0) {
    for(int i = 0; i < nodes.size(); i++) {
      int next = nodes[i] - 1;
      if(is_parent) {
        if(!result[next]) {
          result[next] = true;
          traverse_otr(list, next, result, background, is_parent);
        }
      } else if(background[node]) {
        if(background[next] && !result[next]) {
          result[next] = true;
          traverse_otr(list, next, result, background, is_parent);
        }
      }
    }
  }
}

inline double calc_value(const List& childMap, const List& edgeTypes,
                         const NumericVector& contribution, int node, int end,
                         const LogicalVector& background, bool correct, double c) {
  if(node == end) return 1.0;
  
  IntegerVector children = childMap[node];
  IntegerVector relations = edgeTypes[node];
  int nc = 0;
  double max_s = 0.0;
  
  for(int i = 0; i < children.size(); i++) {
    if(background[children[i] - 1]) {
      nc++;
      double child_s = calc_value(childMap, edgeTypes, contribution,
                                  children[i] - 1, end, background, correct, c);
      double s = child_s * (correct ? (1/(c + nc) + contribution[relations[i] - 1]) 
                              : contribution[relations[i] - 1]);
      if(s > max_s) max_s = s;
    }
  }
  return max_s;
}

// [[Rcpp::export]]
NumericMatrix wang_similarity(S4 otr, IntegerVector nodes, NumericVector contribution, bool correct = false) {
  List childMap = otr.slot("childMap");
  List edgeTypes = otr.slot("edgeTypes");
  int n = childMap.size();
  int m = nodes.size();
  NumericMatrix sim(m, m);
  
  if(m <= 1) {
    if(m == 1) sim(0,0) = 1.0;
    return sim;
  }
  
  IntegerVector node_map(n, -1);
  for(int i = 0; i < m; i++) node_map[nodes[i] - 1] = i;
  
  LogicalVector ancestors(n);
  for(int i = 0; i < nodes.size(); i++) {
    ancestors[nodes[i] - 1] = true;
    traverse_otr(otr.slot("parentMap"), nodes[i] - 1, ancestors, LogicalVector(0), true);
  }
  
  LogicalVector offspring(n);
  NumericVector ic(m);
  double c = correct ? max(contribution)/(1 - max(contribution)) : 0.0;
  
  IntegerVector ancestor_indices = find_indices(ancestors);
  for(int k = 0; k < ancestor_indices.size(); k++) {
    std::fill(offspring.begin(), offspring.end(), false);
    int ak = ancestor_indices[k];
    offspring[ak] = true;
    traverse_otr(childMap, ak, offspring, ancestors, false);
    
    IntegerVector cur_offspring = find_indices(offspring);
    for(int i = 0; i < cur_offspring.size(); i++) {
      int id1 = node_map[cur_offspring[i]];
      if(id1 >= 0) {
        ic[id1] += calc_value(childMap, edgeTypes, contribution,
                              ak, cur_offspring[i], ancestors, correct, c);
        
        for(int j = i + 1; j < cur_offspring.size(); j++) {
          int id2 = node_map[cur_offspring[j]];
          if(id2 >= 0) {
            sim(id1, id2) += calc_value(childMap, edgeTypes, contribution,
                ak, cur_offspring[i], ancestors, correct, c) +
                  calc_value(childMap, edgeTypes, contribution,
                             ak, cur_offspring[j], ancestors, correct, c);
          }
        }
      }
    }
  }
  
  for(int i = 0; i < m - 1; i++) {
    sim(i,i) = 1.0;
    for(int j = i + 1; j < m; j++) {
      sim(i,j) /= (ic[i] + ic[j]);
      sim(j,i) = sim(i,j);
    }
  }
  sim(m-1,m-1) = 1.0;
  
  return sim;
}

// [[Rcpp::export]]
NumericMatrix wang_similarity_sv(const NumericMatrix& sv) {
  int na = sv.nrow();
  int n = sv.ncol();
  NumericMatrix sim(n, n);
  
  if(n <= 1) {
    if(n == 1) sim(0,0) = 1.0;
    return sim;
  }
  
  // Pre-compute information content
  NumericVector ic(n);
  for(int i = 0; i < n; i++) {
    double sum = 0.0;
    for(int j = 0; j < na; j++) sum += sv(j, i);
    ic[i] = sum;
  }
  
  // Fill diagonal and compute similarity
  for(int i = 0; i < n - 1; i++) {
    sim(i,i) = 1.0;
    const double* sv_i = &sv(0,i);
    
    for(int j = i + 1; j < n; j++) {
      double sum = 0.0;
      const double* sv_j = &sv(0,j);
      
      for(int k = 0; k < na; k++) {
        if(std::abs(sv_i[k]) > 1e-10 && std::abs(sv_j[k]) > 1e-10) {
          sum += sv_i[k] + sv_j[k];
        }
      }
      
      sim(i,j) = sum/(ic[i] + ic[j]);
      sim(j,i) = sim(i,j);
    }
  }
  sim(n-1,n-1) = 1.0;
  
  return sim;
}
