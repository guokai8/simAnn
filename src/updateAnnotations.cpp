#include <Rcpp.h>
#include <sstream>
#include <unordered_set>
#include <map>
#include <vector>
#include <algorithm>
#include <cctype>
#include <cfloat>

using namespace Rcpp;

// -----------------------------------------------------------------------------
// 1) Helper function: trim whitespace
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
std::string trim(const std::string &s) {
  size_t start = s.find_first_not_of(" \t\n\r");
  size_t end = s.find_last_not_of(" \t\n\r");
  return (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
}

// -----------------------------------------------------------------------------
// 2) Helper function: split a string by a fixed delimiter
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
std::vector<std::string> splitString(const std::string &s, char delim) {
  std::vector<std::string> tokens;
  std::istringstream iss(s);
  std::string token;
  while (std::getline(iss, token, delim)) {
    tokens.push_back(trim(token));
  }
  return tokens;
}

// -----------------------------------------------------------------------------
// 3) Helper function: union of two CharacterVectors
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
CharacterVector unionChar(const CharacterVector &v1, const CharacterVector &v2) {
  std::unordered_set<std::string> s;
  for (int i = 0; i < v1.size(); i++) {
    s.insert(as<std::string>(v1[i]));
  }
  for (int i = 0; i < v2.size(); i++) {
    s.insert(as<std::string>(v2[i]));
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
  
  // Keep only valid terms that exist in tree@termNames
  CharacterVector treeTermNames = tree.slot("termNames");
  CharacterVector annotNames = annot_list.attr("names");
  std::vector<std::string> valid_terms;
  
  for (int i = 0; i < annotNames.size(); i++) {
    std::string term = as<std::string>(annotNames[i]);
    // Check if 'term' is in treeTermNames
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
  IntegerVector term_indices = term2id(tree, wrap(valid_terms));

  // Access tree@annotations (S4), get the 'list' and 'names' slots
  List annotSlot = tree.slot("annotations");

  List annot_map = annotSlot["list"];     // integer vectors
  CharacterVector treeAnnNames = annotSlot["names"]; // gene ID strings
  
  // If the 'names' slot is NULL, we set it to empty
  if (Rf_isNull(treeAnnNames)) {
    treeAnnNames = CharacterVector();
  }
  
  // Convert each element of annot_map from integer indices (1-based) to CharacterVector of gene IDs
  List annot_map_char(annot_map.size());
  for (int i = 0; i < annot_map.size(); i++) {
    IntegerVector idxVec = annot_map[i];
    int m = idxVec.size();
    CharacterVector genes(m);
    for (int j = 0; j < m; j++) {
      int idx = idxVec[j]; // 1-based index
      if (idx < 1 || idx > treeAnnNames.size()) {
        genes[j] = "";
      } else {
        genes[j] = treeAnnNames[idx - 1]; // convert 1-based -> 0-based
      }
    }
    annot_map_char[i] = genes;
  }
  
  // For each valid term, either merge or replace existing gene IDs
  if (use_background) {
    // Merge the old gene IDs with the new gene IDs
    for (int i = 0; i < term_indices.size(); i++) {
      int tid = term_indices[i]; // 1-based
      CharacterVector old_ids = as<CharacterVector>(annot_map_char[tid]);
      CharacterVector new_ids = as<CharacterVector>(newAnnotList[valid_terms[i]]);
      CharacterVector merged = union_(old_ids, new_ids);
      annot_map_char[tid] = merged;
    }
  } else {
    // Replace
    for (int i = 0; i < term_indices.size(); i++) {
      int tid = term_indices[i];
      annot_map_char[tid] = as<CharacterVector>(newAnnotList[valid_terms[i]]);
    }
  }
  
  // Now unify the updated annotation map. Each element is a CharacterVector of gene IDs.
  // Build a single vector all_genes, sort/unique to get a unified set.
  std::vector<std::string> all_genes;
  for (int i = 0; i < annot_map_char.size(); i++) {
    CharacterVector genes = annot_map_char[i];
    for (int j = 0; j < genes.size(); j++) {
      all_genes.push_back(as<std::string>(genes[j]));
    }
  }
  std::sort(all_genes.begin(), all_genes.end());
  all_genes.erase(std::unique(all_genes.begin(), all_genes.end()), all_genes.end());
  
  // Create a map from gene ID -> 1-based index
  std::map<std::string,int> gene_index_map;
  for (size_t i = 0; i < all_genes.size(); i++) {
    gene_index_map[ all_genes[i] ] = i + 1;
  }
  
  // Convert each element in annot_map_char from gene IDs to indices
  List new_annot_map(annot_map_char.size());
  for (int i = 0; i < annot_map_char.size(); i++) {
    CharacterVector genes = annot_map_char[i];
    int len = genes.size();
    IntegerVector idxVec(len);
    for (int j = 0; j < len; j++) {
      std::string g = as<std::string>(genes[j]);
      auto it = gene_index_map.find(g);
      if (it != gene_index_map.end()) {
        idxVec[j] = it->second;
      } else {
        idxVec[j] = 0; // or some sentinel if not found
      }
    }
    new_annot_map[i] = idxVec;
  }
  
  // Finally, update the tree@annotations slot
  annotSlot["list"] = new_annot_map;
  annotSlot["names"] = wrap(all_genes);
  tree.slot("annotations") = annotSlot;
  
  return tree;
}


