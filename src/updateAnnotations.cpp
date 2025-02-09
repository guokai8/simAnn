#include <Rcpp.h>
#include <sstream>
#include <unordered_set>
#include <map>
#include <vector>
#include <algorithm>
#include <cctype>
#include <cfloat>
using namespace Rcpp;
// [[Rcpp::export]]
std::vector<std::string> splitString(const std::string &s, char delim) {
  std::vector<std::string> tokens;
  std::istringstream iss(s);
  std::string token;
  while (std::getline(iss, token, delim)) {
    if (!token.empty()) {
      tokens.push_back(token);
    }
  }
  return tokens;
}

// -----------------------------------------------------------------------------
// 2) Helper function: union of two CharacterVectors
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
CharacterVector unionChar(const CharacterVector &v1, const CharacterVector &v2) {
  std::unordered_set<std::string> s;
  for (auto str : v1) {
    s.insert(as<std::string>(str));
  }
  for (auto str : v2) {
    s.insert(as<std::string>(str));
  }
  CharacterVector res(s.size());
  int idx = 0;
  for (auto &str : s) {
    res[idx++] = str;
  }
  return res;
}

// [[Rcpp::export]]
S4 updateAnnotations(S4 tree, SEXP new_annot, bool use_background = false) {
  // Clone the original tree to avoid modifying it
  S4 new_tree = clone(tree);
  
  // Check that 'tree' is an OntologyTree
  if (!tree.inherits("OntologyTree")) {
    stop("'tree' must be an OntologyTree object");
  }
  
  // Convert new_annot from data.frame to list if needed
  List annot_list;
  if (Rf_inherits(new_annot, "data.frame")) {
    DataFrame df(new_annot);
    CharacterVector annotCol = df["Annot"];
    CharacterVector geneIDCol = df["GeneID"];
    int n = geneIDCol.size();
    annot_list = List(n);
    CharacterVector listNames(n);
    
    for (int i = 0; i < n; i++) {
      std::string term = as<std::string>(annotCol[i]);
      listNames[i] = term;
      std::string genes_str = as<std::string>(geneIDCol[i]);
      std::vector<std::string> tokens = splitString(genes_str, ',');
      CharacterVector cv(tokens.size());
      for (size_t j = 0; j < tokens.size(); j++) {
        cv[j] = tokens[j];
      }
      annot_list[i] = cv;
    }
    annot_list.attr("names") = listNames;
  } else {
    annot_list = as<List>(new_annot);
  }
  
  // Keep only valid terms that exist in new_tree@termNames
  CharacterVector treeTermNames = new_tree.slot("termNames");
  CharacterVector annotNames = annot_list.attr("names");
  std::vector<std::string> valid_terms;
  
  for (int i = 0; i < annotNames.size(); i++) {
    std::string term = as<std::string>(annotNames[i]);
    for (int j = 0; j < treeTermNames.size(); j++) {
      if (term == as<std::string>(treeTermNames[j])) {
        valid_terms.push_back(term);
        break;
      }
    }
  }
  if (valid_terms.empty()) {
    stop("No valid terms found");
  }
  
  // Subset annot_list to those valid_terms only
  List newAnnotList(valid_terms.size());
  for (size_t i = 0; i < valid_terms.size(); i++) {
    newAnnotList[valid_terms[i]] = annot_list[valid_terms[i]];
  }
  
  // Use R's term2id function to get 1-based term indices for valid_terms
  Function term2id("term2id");
  IntegerVector term_indices = term2id(new_tree, wrap(valid_terms));
  
  // Access new_tree@annotations
  List annotSlot = new_tree.slot("annotations");
  List annot_map = annotSlot["list"];     // integer vectors
  CharacterVector treeAnnNames = annotSlot["names"]; // gene ID strings
  
  // Create a map from integer index to gene IDs
  List annot_map_working = clone(annot_map);
  
  // For each valid term, either merge or replace existing gene IDs
  if (use_background) {
    // Merge the old gene IDs with the new gene IDs
    for (int i = 0; i < term_indices.size(); i++) {
      int tid = term_indices[i] - 1; // Convert to 0-based
      IntegerVector old_indices = annot_map[tid];
      CharacterVector old_ids(old_indices.size());
      
      // Convert indices to gene IDs
      for(int j = 0; j < old_indices.size(); j++) {
        if(old_indices[j] > 0 && old_indices[j] <= treeAnnNames.size()) {
          old_ids[j] = treeAnnNames[old_indices[j] - 1];
        }
      }
      
      CharacterVector new_ids = as<CharacterVector>(newAnnotList[valid_terms[i]]);
      CharacterVector merged = unionChar(old_ids, new_ids);
      annot_map_working[tid] = merged;
    }
  } else {
    // Replace
    for (int i = 0; i < term_indices.size(); i++) {
      int tid = term_indices[i] - 1; // Convert to 0-based
      annot_map_working[tid] = as<CharacterVector>(newAnnotList[valid_terms[i]]);
    }
  }
  
  // Now unify the updated annotation map
  std::vector<std::string> all_genes;
  for (int i = 0; i < annot_map_working.size(); i++) {
    CharacterVector genes = as<CharacterVector>(annot_map_working[i]);
    for (int j = 0; j < genes.size(); j++) {
      std::string gene = as<std::string>(genes[j]);
      if (!gene.empty()) {
        all_genes.push_back(gene);
      }
    }
  }
  std::sort(all_genes.begin(), all_genes.end());
  all_genes.erase(std::unique(all_genes.begin(), all_genes.end()), all_genes.end());
  
  // Create a map from gene ID -> 1-based index
  std::map<std::string, int> gene_index_map;
  for (size_t i = 0; i < all_genes.size(); i++) {
    gene_index_map[all_genes[i]] = i + 1;
  }
  
  // Convert back to integer indices
  List new_annot_map(annot_map.size());
  for (int i = 0; i < annot_map_working.size(); i++) {
    CharacterVector genes = as<CharacterVector>(annot_map_working[i]);
    IntegerVector indices(genes.size());
    for (int j = 0; j < genes.size(); j++) {
      std::string gene = as<std::string>(genes[j]);
      if (!gene.empty()) {
        indices[j] = gene_index_map[gene];
      }
    }
    new_annot_map[i] = indices;
  }
  
  // Update the new_tree's annotations slots
  annotSlot["list"] = new_annot_map;
  annotSlot["names"] = wrap(all_genes);
  new_tree.slot("annotations") = annotSlot;
  
  return new_tree;
}