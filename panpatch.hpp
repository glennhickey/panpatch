#pragma once

#include <unordered_map>
#include "handlegraph/path_handle_graph.hpp"

using namespace std;
using namespace handlegraph;

// compute overlap identity between set of ref_paths (treated colletively) with each path in the
// other_paths 
// the identity is computed in windows of size w. for a given window the identity is
// #exact matchces / w
// the overall identity  is average over all windows, excluding empty windows
unordered_map<path_handle_t, double> compute_overlap_identity(const PathHandleGraph* graph,
                                                              const vector<path_handle_t>& tgt_paths,
                                                              const vector<path_handle_t>& other_paths,
                                                              int64_t w=1000);

// process the above coverage map to select which paths for each sample best cover the input
// path.  in other words, it tries to find the best haplotype for each other sample
unordered_map<string, vector<path_handle_t>> select_sample_covers(const PathHandleGraph* graph,
                                                                  const unordered_map<path_handle_t, double>& coverage_map);

// sort contigs against a given reference genome
// it's expected that the other_paths are all part of the same sample and haplotype
// (see above function to parse them out)
multimap<pair<int64_t, int64_t>, path_handle_t> sort_overlapping_paths(const PathHandleGraph* graph,
                                                                       const path_handle_t& ref_path,
                                                                       const vector<path_handle_t>& other_paths);
// 
void greedy_patch(const PathHandleGraph* graph,
                  const path_handle_t& ref_path,
                  const vector<path_handle_t>& tgt_paths,
                  const vector<string>& sample_names,
                  const unordered_map<string, vector<path_handle_t>>& sample_covers);
                  
                  
