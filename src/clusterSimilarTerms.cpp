#include <Rcpp.h>
#include <sstream>
#include <unordered_set>
#include <map>
#include <vector>
#include <algorithm>
#include <cfloat>  // for DBL_MAX
using namespace Rcpp;

std::vector<std::string> getGeneIDsForTerms(const std::vector<int>& termIds,
                                            const List &annotationsList,
                                            const CharacterVector &annotationsNames) {
  std::unordered_set<std::string> geneSet;
  for (size_t i = 0; i < termIds.size(); i++) {
    int tid = termIds[i]; // 1-based
    if (tid > 0 && tid <= annotationsList.size()) {
      IntegerVector geneIdx = annotationsList[tid - 1];  // Convert to 0-based
      for (int j = 0; j < geneIdx.size(); j++) {
        int gidx = geneIdx[j]; // Assume gene index is also 1-based
        if (gidx > 0 && gidx <= annotationsNames.size()) {
          geneSet.insert(as<std::string>(annotationsNames[gidx - 1]));
        }
      }
    }
  }
  std::vector<std::string> genes;
  for (auto &g : geneSet)
    genes.push_back(g);
  return genes;
}

// [[Rcpp::export]]
List clusterSimilarTermsCPP(
    SEXP sim_matrix,                  // sim_matrix is a data.frame (or matrix) with row and column names (this code does not use dimnames to obtain names)
    const IntegerVector &terms_vec,   // 1-based term ID vector (e.g., obtained from term2id)
    const List &parentMap,            // Parent map (1-based, index: termID-1)
    const CharacterVector &termNames, // Character vector of descriptive term names (1-based, index: termID-1), whose names attribute stores the GO id
    const NumericVector &icScores,    // IC values (1-based, index: termID-1)
    const List &annotationsList,      // tree@annotations$list, with 1-based gene indices
    const CharacterVector &annotationsNames, // tree@annotations$names
    const double threshold = 0.8
) {
  // Convert sim_matrix to NumericMatrix
  NumericMatrix simMat = as<NumericMatrix>(sim_matrix);
  
  // Try to obtain the names attribute of termNames to use as the GO id vector
  CharacterVector goIds = termNames.attr("names");
  
  // ---------------------------
  // 1) Initialize containers for cluster results
  // ---------------------------
  std::map<std::string, List> parentClusters;
  std::map<std::string, List> singletonClusters;
  
  // Record processed term IDs (1-based)
  std::unordered_set<int> processedTerms;
  processedTerms.reserve(terms_vec.size());
  
  // ---------------------------
  // 2) Iterate through terms_vec (1-based)
  // ---------------------------
  for (int i = 0; i < terms_vec.size(); i++) {
    int currentTermId = terms_vec[i]; // 1-based
    if (processedTerms.find(currentTermId) != processedTerms.end())
      continue;
    
    // Current term descriptive name: obtained directly from termNames
    std::string curTermName = "";
    if (currentTermId > 0 && currentTermId <= termNames.size())
      curTermName = as<std::string>(termNames[currentTermId - 1]);
    
    // Current term GO id: obtained from goIds (if available)
    std::string curGoId = "";
    if (!Rf_isNull(goIds) && goIds.size() >= (size_t)currentTermId)
      curGoId = as<std::string>(goIds[currentTermId - 1]);
    else
      curGoId = curTermName;
    
    // If the term is already in parentClusters, skip it
    if (parentClusters.find(curTermName) != parentClusters.end())
      continue;
    
    // Get the similarities from the i-th row of simMat (assumes simMat row order corresponds to terms_vec)
    std::vector<double> similarities(terms_vec.size());
    for (int col = 0; col < terms_vec.size(); col++) {
      similarities[col] = simMat(i, col);
    }
    
    // Find indices with similarity >= threshold (excluding itself)
    std::vector<int> similarIndices;
    similarIndices.reserve(terms_vec.size());
    for (int col = 0; col < (int)similarities.size(); col++) {
      if (col == i) continue;
      if (similarities[col] >= threshold)
        similarIndices.push_back(col);
    }
    
    // ---------------------------
    // 2.1 If similar terms are found, build a cluster
    // ---------------------------
    if (!similarIndices.empty()) {
      // Merge the current term with similar terms; store in allTermIds (all 1-based)
      std::vector<int> allTermIds;
      allTermIds.push_back(currentTermId);
      for (auto idx : similarIndices)
        allTermIds.push_back(terms_vec[idx]);
      
      // Extract the gene IDs corresponding to this cluster
      std::vector<std::string> clusterGeneIDs = getGeneIDsForTerms(allTermIds, annotationsList, annotationsNames);
      
      // Determine if a parent-child relationship exists: iterate through allTermIds to check if a term appears in another term's parent list
      bool parentFound = false;
      int parentTermId = -1;
      for (int tid1 : allTermIds) {
        if (tid1 > 0 && tid1 <= parentMap.size()) {
          IntegerVector pvec = parentMap[tid1 - 1];
          std::unordered_set<int> pset;
          for (int k = 0; k < pvec.size(); k++)
            pset.insert(pvec[k]);
          for (int tid2 : allTermIds) {
            if (tid1 == tid2) continue;
            if (pset.find(tid2) != pset.end()) {
              parentTermId = tid2;
              parentFound = true;
              break;
            }
          }
          if (parentFound)
            break;
        }
      }
      // If no parent is found, choose the term with the smallest IC value as the parent
      if (!parentFound) {
        double minIC = DBL_MAX;
        int bestId = -1;
        for (int tid : allTermIds) {
          if (tid > 0 && tid <= icScores.size()) {
            double ic = icScores[tid - 1];
            if (ic < minIC) {
              minIC = ic;
              bestId = tid;
            }
          }
        }
        parentTermId = bestId;
      }
      
      // Child terms: allTermIds excluding the parent
      std::vector<int> childTermIds;
      for (int tid : allTermIds)
        if (tid != parentTermId)
          childTermIds.push_back(tid);
        
        // Compute the similarity between each child and the parent
        std::vector<double> childSimilarities;
        childSimilarities.reserve(childTermIds.size());
        int pidx = -1;
        for (int idx = 0; idx < terms_vec.size(); idx++) {
          if (terms_vec[idx] == parentTermId) { pidx = idx; break; }
        }
        for (int cid : childTermIds) {
          int cidx = -1;
          for (int idx = 0; idx < terms_vec.size(); idx++) {
            if (terms_vec[idx] == cid) { cidx = idx; break; }
          }
          if (pidx >= 0 && cidx >= 0)
            childSimilarities.push_back(simMat(pidx, cidx));
          else
            childSimilarities.push_back(0.0);
        }
        
        // Get parent node information:
        // Parent GO id: obtained from goIds
        std::string parentGoId = "";
        if (parentTermId > 0 && parentTermId <= termNames.size()) {
          if (!Rf_isNull(goIds) && goIds.size() >= (size_t)parentTermId)
            parentGoId = as<std::string>(goIds[parentTermId - 1]);
          else
            parentGoId = as<std::string>(termNames[parentTermId - 1]);
        }
        // Parent descriptive name: obtained from termNames
        std::string parentTermName = "";
        if (parentTermId > 0 && parentTermId <= termNames.size())
          parentTermName = as<std::string>(termNames[parentTermId - 1]);
        
        // Get all child node information: both GO id and descriptive name from termNames
        std::vector<std::string> childGoIds;
        std::vector<std::string> childTermNames;
        childGoIds.reserve(childTermIds.size());
        childTermNames.reserve(childTermIds.size());
        for (int cid : childTermIds) {
          std::string goid = "";
          std::string tname = "";
          if (cid > 0 && cid <= termNames.size()) {
            if (!Rf_isNull(goIds) && goIds.size() >= (size_t)cid)
              goid = as<std::string>(goIds[cid - 1]);
            else
              goid = as<std::string>(termNames[cid - 1]);
            tname = as<std::string>(termNames[cid - 1]);
          }
          childGoIds.push_back(goid);
          childTermNames.push_back(tname);
        }
        
        // Calculate the mean and maximum similarity of the child nodes
        double sumSim = 0.0, maxSim = -DBL_MAX;
        for (double v : childSimilarities) {
          sumSim += v;
          if (v > maxSim)
            maxSim = v;
        }
        double meanSim = childSimilarities.empty() ? 0.0 : (sumSim / childSimilarities.size());
        
        // When merging with an existing cluster, handle new child nodes as follows:
        // According to the R logic: new_similarities <- sim_matrix[match(parent_term_id, terms_vec), match(new_children, termNames)]
        auto it = parentClusters.find(parentTermName);
        if (it != parentClusters.end()) {
          List existingCluster = it->second;
          CharacterVector existingChildren = existingCluster["children"]; // Existing child GO ids
          NumericVector existingSim = existingCluster["child_similarities"];
          CharacterVector existingGeneIDs = existingCluster["GeneID"];
          CharacterVector existingChildNames = existingCluster["child_names"]; // Descriptive names
          
          std::vector<std::string> mergedChildGoIds = as< std::vector<std::string> >(existingChildren);
          std::vector<std::string> mergedChildNames = as< std::vector<std::string> >(existingChildNames);
          std::unordered_set<std::string> existingChildrenSet(mergedChildGoIds.begin(), mergedChildGoIds.end());
          std::vector<double> newSimilaritiesCombined = as< std::vector<double> >(existingSim);
          
          // Get the index of the parent in terms_vec
          int pIndex = -1;
          for (int idx = 0; idx < terms_vec.size(); idx++) {
            if (terms_vec[idx] == parentTermId) { pIndex = idx; break; }
          }
          // Iterate over the current child nodes (childTermIds corresponding to childGoIds) and add those not already present
          for (size_t k = 0; k < childTermIds.size(); k++) {
            // Get the GO id of the child
            std::string childGoid = "";
            if (childTermIds[k] > 0 && childTermIds[k] <= termNames.size()) {
              if (!Rf_isNull(goIds) && goIds.size() >= (size_t)childTermIds[k])
                childGoid = as<std::string>(goIds[childTermIds[k] - 1]);
              else
                childGoid = as<std::string>(termNames[childTermIds[k] - 1]);
            }
            if (existingChildrenSet.find(childGoid) == existingChildrenSet.end()) {
              mergedChildGoIds.push_back(childGoid);
              // Get the descriptive name directly from termNames based on childTermIds
              std::string childDesc = "";
              if (childTermIds[k] > 0 && childTermIds[k] <= termNames.size())
                childDesc = as<std::string>(termNames[childTermIds[k] - 1]);
              mergedChildNames.push_back(childDesc);
              // Compute the similarity between the new child and the parent using terms_vec indices
              int cIndex = -1;
              for (int idx = 0; idx < terms_vec.size(); idx++) {
                if (terms_vec[idx] == childTermIds[k]) { cIndex = idx; break; }
              }
              double newSim = 0.0;
              if (pIndex != -1 && cIndex != -1)
                newSim = simMat(pIndex, cIndex);
              newSimilaritiesCombined.push_back(newSim);
            }
          }
          
          double sumSim2 = 0.0, maxSim2 = -DBL_MAX;
          for (double v : newSimilaritiesCombined) {
            sumSim2 += v;
            if (v > maxSim2)
              maxSim2 = v;
          }
          double meanSim2 = newSimilaritiesCombined.empty() ? 0.0 : (sumSim2 / newSimilaritiesCombined.size());
          
          // Merge gene IDs: take the union of existingGeneIDs and clusterGeneIDs
          std::vector<std::string> existingGeneVec = as< std::vector<std::string> >(existingGeneIDs);
          std::vector<std::string> clusterGeneIDsNew = getGeneIDsForTerms(allTermIds, annotationsList, annotationsNames);
          std::unordered_set<std::string> geneSet(existingGeneVec.begin(), existingGeneVec.end());
          for (auto &g : clusterGeneIDsNew)
            geneSet.insert(g);
          std::vector<std::string> mergedGeneIDs;
          for (auto &g : geneSet)
            mergedGeneIDs.push_back(g);
          
          std::string oldRelation = as<std::string>(existingCluster["relation_type"]);
          std::string newRelation = parentFound ? "parent-child" : oldRelation;
          
          List updatedCluster = List::create(
            _["cluster_term"] = parentGoId,         // Parent GO id
            _["term_name"]    = parentTermName,       // Parent descriptive name
            _["mean_similarity"] = meanSim2,
            _["max_similarity"]  = maxSim2,
            _["relation_type"] = newRelation,
            _["children"] = wrap(mergedChildGoIds),    // Child GO ids
            _["child_similarities"] = newSimilaritiesCombined,
            _["child_names"] = wrap(mergedChildNames), // Child descriptive names
            _["GeneID"] = wrap(mergedGeneIDs)
          );
          parentClusters[parentTermName] = updatedCluster;
        } else {
          std::string relType = parentFound ? "parent-child" : "similarity-based";
          CharacterVector childGoIdsR = wrap(childGoIds);   // Child GO ids
          NumericVector childSimsR(childSimilarities.begin(), childSimilarities.end());
          CharacterVector childNamesR = wrap(childTermNames); // Child descriptive names
          
          std::vector<std::string> clusterGeneIDsNew = getGeneIDsForTerms(allTermIds, annotationsList, annotationsNames);
          
          List newCluster = List::create(
            _["cluster_term"] = parentGoId,      // Parent GO id
            _["term_name"]    = parentTermName,    // Parent descriptive name
            _["mean_similarity"] = meanSim,
            _["max_similarity"]  = maxSim,
            _["relation_type"] = relType,
            _["children"] = childGoIdsR,         // Child GO ids
            _["child_similarities"] = childSimsR,
            _["child_names"] = childNamesR,        // Child descriptive names
            _["GeneID"] = wrap(clusterGeneIDsNew)
          );
          parentClusters[parentTermName] = newCluster;
        }
        
        // Mark all processed terms (1-based)
        for (int tid : allTermIds)
          processedTerms.insert(tid);
        
    } else {
      // ---------------------------
      // 2.2 If no similar term is found, then treat it as a singleton cluster
      // ---------------------------
      std::vector<int> singleTermIds;
      singleTermIds.push_back(currentTermId);
      std::vector<std::string> singleGeneIDs = getGeneIDsForTerms(singleTermIds, annotationsList, annotationsNames);
      
      List singletonCluster = List::create(
        _["cluster_term"] = curGoId,          // GO id
        _["term_name"]    = curTermName,        // Descriptive name
        _["mean_similarity"] = 0.0,
        _["max_similarity"] = 0.0,
        _["relation_type"] = "singleton",
        _["children"] = CharacterVector::create(curGoId),    // Child is the same GO id
        _["child_similarities"] = NumericVector::create(0.0),
        _["child_names"] = CharacterVector::create(curTermName), // Descriptive name
        _["GeneID"] = wrap(singleGeneIDs)
      );
      singletonClusters[curTermName] = singletonCluster;
      processedTerms.insert(currentTermId);
    }
  }
  
  // ---------------------------
  // 3) Merge parentClusters and singletonClusters
  // ---------------------------
  std::vector<List> clusters;
  clusters.reserve(parentClusters.size() + singletonClusters.size());
  for (auto &kv : parentClusters)
    clusters.push_back(kv.second);
  for (auto &kv : singletonClusters)
    clusters.push_back(kv.second);
  
  // ---------------------------
  // 4) Sort clusters (by relation_type and mean_similarity)
  // ---------------------------
  auto getRelationTypeValue = [&](const std::string &rt) {
    if (rt == "singleton")
      return 2;
    else if (rt == "parent-child")
      return 0;
    else
      return 1;
  };
  std::sort(clusters.begin(), clusters.end(), [&](const List &A, const List &B) {
    std::string ra = as<std::string>(A["relation_type"]);
    std::string rb = as<std::string>(B["relation_type"]);
    double ma = as<double>(A["mean_similarity"]);
    double mb = as<double>(B["mean_similarity"]);
    int va = getRelationTypeValue(ra);
    int vb = getRelationTypeValue(rb);
    if (va != vb)
      return va < vb;
    else {
      if (ma != mb)
        return ma > mb;
      else
        return false;
    }
  });
  
  int n = clusters.size();
  CharacterVector summary_cluster_term(n);
  NumericVector summary_mean_similarity(n);
  NumericVector summary_max_similarity(n);
  CharacterVector summary_relation_type(n);
  // Combine the children vector into a single string
  CharacterVector summary_children_str(n);
  CharacterVector summary_geneID(n);
  CharacterVector summary_term_name(n);
  
  for (int i = 0; i < n; i++) {
    List cl = clusters[i];
    summary_cluster_term[i] = as<std::string>(cl["cluster_term"]);
    summary_mean_similarity[i] = as<double>(cl["mean_similarity"]);
    summary_max_similarity[i] = as<double>(cl["max_similarity"]);
    summary_relation_type[i] = as<std::string>(cl["relation_type"]);
    
    CharacterVector children = cl["children"];
    std::stringstream ss2;
    for (int k = 0; k < children.size(); k++) {
      if (k > 0)
        ss2 << ",";
      ss2 << as<std::string>(children[k]);
    }
    summary_children_str[i] = ss2.str();
    
    CharacterVector genes = cl["GeneID"];
    std::stringstream ss;
    for (int k = 0; k < genes.size(); k++) {
      if (k > 0)
        ss << ",";
      ss << as<std::string>(genes[k]);
    }
    summary_geneID[i] = ss.str();
    
    summary_term_name[i] = as<std::string>(cl["term_name"]);
  }
  
  DataFrame summaryDF = DataFrame::create(
    _["cluster_term"]   = summary_cluster_term,
    _["mean_similarity"] = summary_mean_similarity,
    _["max_similarity"] = summary_max_similarity,
    _["relation_type"]  = summary_relation_type,
    _["children"]       = summary_children_str,
    _["GeneID"]         = summary_geneID,
    _["term_name"]      = summary_term_name
  );
  summaryDF.attr("stringsAsFactors") = false;
  
  // Construct cluster_matrix: record parent/child/similarity information row by row
  std::vector<std::string> v_parent, v_parent_name, v_child, v_child_name, v_relation_type, v_geneID;
  std::vector<double> v_similarity;
  
  for (int iCl = 0; iCl < n; iCl++) {
    List cl = clusters[iCl];
    std::string parent = as<std::string>(cl["cluster_term"]);
    std::string parentName = as<std::string>(cl["term_name"]);
    std::string relation = as<std::string>(cl["relation_type"]);
    
    CharacterVector cchildren = cl["children"];
    CharacterVector cchildNames = cl["child_names"];
    NumericVector cchildSims = cl["child_similarities"];
    CharacterVector cgenes = cl["GeneID"];
    
    std::stringstream ssg;
    for (int k = 0; k < cgenes.size(); k++) {
      if (k > 0)
        ssg << ",";
      ssg << as<std::string>(cgenes[k]);
    }
    std::string genesAll = ssg.str();
    
    if (relation == "singleton") {
      v_parent.push_back(parent);
      v_parent_name.push_back(parentName);
      v_child.push_back(as<std::string>(cchildren[0]));
      v_child_name.push_back(as<std::string>(cchildNames[0]));
      v_similarity.push_back(0.0);
      v_relation_type.push_back(relation);
      v_geneID.push_back(genesAll);
    } else {
      for (int k = 0; k < cchildren.size(); k++) {
        v_parent.push_back(parent);
        v_parent_name.push_back(parentName);
        v_child.push_back(as<std::string>(cchildren[k]));
        v_child_name.push_back(as<std::string>(cchildNames[k]));
        v_similarity.push_back(cchildSims[k]);
        v_relation_type.push_back(relation);
        v_geneID.push_back(genesAll);
      }
    }
  }
  
  DataFrame clusterMat = DataFrame::create(
    _["parent"] = v_parent,
    _["parent_name"] = v_parent_name,
    _["child"] = v_child,
    _["child_name"] = v_child_name,
    _["similarity"] = v_similarity,
    _["relation_type"] = v_relation_type,
    _["GeneID"] = v_geneID
  );
  
  int total_clusters = n;
  int parent_child_clusters = 0;
  int similarity_clusters = 0;
  int singleton_terms = 0;
  double sumMeanSim = 0.0;
  
  for (int iCl = 0; iCl < n; iCl++) {
    List cl = clusters[iCl];
    std::string relation = as<std::string>(cl["relation_type"]);
    double ms = as<double>(cl["mean_similarity"]);
    sumMeanSim += ms;
    if (relation == "parent-child")
      parent_child_clusters++;
    else if (relation == "similarity-based")
      similarity_clusters++;
    else if (relation == "singleton")
      singleton_terms++;
  }
  double meanSimilarityAll = (n == 0) ? 0.0 : (sumMeanSim / n);
  
  List stats = List::create(
    _["total_clusters"] = total_clusters,
    _["parent_child_clusters"] = parent_child_clusters,
    _["similarity_clusters"] = similarity_clusters,
    _["singleton_terms"] = singleton_terms,
    _["mean_similarity"] = NumericVector::create(meanSimilarityAll)
  );
  
  List result = List::create(
    _["clusters"] = clusters,
    _["summary"] = summaryDF,
    _["cluster_matrix"] = clusterMat,
    _["similarity_matrix"] = sim_matrix, // Returned as is
    _["stats"] = stats
  );
  
  result.attr("method") = "convertedCpp";
  result.attr("threshold") = threshold;
  result.attr("timestamp") = Datetime(::time(nullptr));
  
  return result;
}
