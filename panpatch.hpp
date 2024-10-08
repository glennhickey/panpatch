#pragma once

#include <unordered_map>
#include "handlegraph/path_handle_graph.hpp"

using namespace std;
using namespace handlegraph;

// compute overlap identity between given ref_path with all other paths that overlap it
// the identity is computed in windows of size w. for a given window the identity is
// #exact matchces / w
// the overall identity  is average over all windows, excluding empty windows
unordered_map<path_handle_t, double> compute_overlap_identity(const PathHandleGraph* graph,
                                                              const path_handle_t& ref_path,
                                                              int64_t w=1000);
                                                              
