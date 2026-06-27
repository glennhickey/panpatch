#pragma once

#include <unordered_map>
#include <map>
#include <string>
#include <fstream>
#include "handlegraph/path_handle_graph.hpp"

using namespace std;
using namespace handlegraph;

// BED regions parsed from file: contig_name -> vector of [start, end) intervals
using BedRegions = unordered_map<string, vector<pair<int64_t, int64_t>>>;

// Sorted, merged intervals in reference coordinates representing excluded regions
using ExcludedRefRegions = vector<pair<int64_t, int64_t>>;

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

// quick and dirty telomere checker!!
// returns position of left and right telomere (wrt to respective ends)
pair<int64_t, int64_t> find_telomeres(const PathHandleGraph* graph,
                                      const path_handle_t path,
                                      double threshold=0.95);

// parse a BED file into regions keyed by contig name
// skips comment/header lines (starting with #, track, browser)
BedRegions parse_bed_file(const string& bed_filename);

// convert BED regions (in target assembly coordinates) to reference coordinate ranges
// using the anchor map to map between coordinate systems
ExcludedRefRegions bed_to_ref_regions(const PathHandleGraph* graph,
                                      const vector<path_handle_t>& tgt_paths,
                                      const unordered_map<int64_t, int64_t>& ref_anchors,
                                      const BedRegions& bed_regions);

// find a series of anchors along the reference path
// anchors have this property:
// 1) they overlap other path(s)
// 2) they are not redundant (ie there are never more than two consecutive anchors)
//    that overlap the same paths
// these anchors will be the basis of the cover search
// maps reference position to anchor
unordered_map<int64_t, int64_t> find_anchors(const PathHandleGraph* graph,
                                             const path_handle_t& ref_path,
                                             const vector<path_handle_t>& tgt_paths,
                                             const unordered_set<path_handle_t>& path_set);

// search along a path to the next anchor
// the bool flag returned is true if search went backwards on path
// direction: -1 backward, 0 both, 1 forward
pair<step_handle_t, bool> find_next_anchor_on_path(const PathHandleGraph* graph,
                                                   const unordered_map<int64_t, int64_t>& anchors,
                                                   step_handle_t step,
                                                   int64_t pos,
                                                   int64_t direction,
                                                   bool cross_Ns,
                                                   int64_t bound=10000);

// thread the anchors using the paths
// return a set of intervals that represent the patched genome
// the bool flag is true if the interval is reversed on the path
// algorithm:
// for every interval, scan along the input paths (in order of priority) until one is found that hits a
// subsequent anchor.  add the interval found and repeat. 
vector<tuple<step_handle_t, step_handle_t, bool>> thread_intervals(const PathHandleGraph* graph,
                                                                   const path_handle_t& ref_path,
                                                                   const unordered_map<int64_t, int64_t>& ref_anchors,
                                                                   const vector<path_handle_t>& tgt_paths,
                                                                   const vector<path_handle_t>& other_paths,
                                                                   const ExcludedRefRegions& excluded_regions = {});

// smooth out the threaded intervals
vector<tuple<step_handle_t, step_handle_t, bool>> smooth_intervals(const PathHandleGraph* graph,
                                                                   const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals);

// extend the intervals to the telomeres where appropriate
vector<tuple<step_handle_t, step_handle_t, bool>> extend_intervals(const PathHandleGraph* graph,
                                                                   const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals);

// patch a missing telomere onto a terminal end by handing off to a foreign cover that has one.
// for each terminal end lacking a telomere, walk inward from the tip to the nearest node also
// visited by a foreign cover, and if that cover runs outward to a real telomere, splice it on
// (replacing the target's capless tail beyond the shared handoff node).
vector<tuple<step_handle_t, step_handle_t, bool>> extend_to_telomeres(const PathHandleGraph* graph,
                                                                      const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals,
                                                                      const unordered_map<string, vector<path_handle_t>>& sample_covers,
                                                                      const vector<string>& sample_names,
                                                                      double telo_threshold=0.8,
                                                                      int64_t max_handoff=500000,
                                                                      bool verbose=false);

// print the intervals in a bed-like format
void print_intervals(const PathHandleGraph* graph,
                     const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals,
                     std::ostream& out = std::cout);

// get the dna sequence of a list of intervals (ie the patched assembly for a contig)
string intervals_to_sequence(const PathHandleGraph* graph,
                             const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals);

// run the patching for a single chromosome (ref_path) for a given haplotype (tgt_paths)
// returns the intervals forming the patched assembly
vector<tuple<step_handle_t, step_handle_t, bool>> greedy_patch(const PathHandleGraph* graph,
                                                               const path_handle_t& ref_path,
                                                               const vector<path_handle_t>& tgt_paths,
                                                               const vector<string>& sample_names,
                                                               const unordered_map<string, vector<path_handle_t>>& sample_covers,
                                                               const BedRegions& bed_regions = {},
                                                               bool patch_ends = false,
                                                               double telo_threshold = 0.8,
                                                               int64_t max_telomere_patch = 500000,
                                                               bool verbose = false);

// partial-patch cleanup: excise repeat-region-misjoin foreign interior grafts in place (merging the
// flanking target pieces), keeping the rest of the patch.  reverts a graft on low k-mer recovery
// (content) or low flank anchoring (locus).  run before revert_bad_patch.
void excise_bad_interior_grafts(const PathHandleGraph* graph,
                                vector<tuple<step_handle_t, step_handle_t, bool>>& intervals,
                                const string& target_sample,
                                double min_recovery,
                                int64_t min_replaced,
                                double min_flank,
                                int64_t flank_window,
                                std::unordered_map<path_handle_t, int64_t>& excised_nonN);

// return the input intervals unmodified if it failed to find a reasonable patch
bool revert_bad_patch(const PathHandleGraph* graph,
                      const path_handle_t& ref_path,
                      const vector<path_handle_t>& tgt_paths,
                      const vector<string>& sample_names,
                      const vector<tuple<step_handle_t, step_handle_t, bool>>& in_intervals,                      
                      vector<tuple<step_handle_t, step_handle_t, bool>>& out_intervals,
                      string default_sample,
                      double threshold,
                      double graft_recovery = 50.0,
                      int64_t graft_min_bp = 10000,
                      double telo_threshold = 0.8,
                      const std::unordered_map<path_handle_t, int64_t>& excised_nonN = {});

// make sure all intervals have the correct orientation (and assert fail if not)
void check_intervals(const PathHandleGraph* graph,
                     const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals);

// validate that intervals have telomeres at both ends and no internal telomeres
// returns true if valid, false otherwise
// if verbose is true, prints diagnostic messages
bool validate_telomeres(const PathHandleGraph* graph,
                        const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals,
                        double threshold=0.8,
                        bool verbose=false);

// log telomere information for each contig in the intervals
// outputs to cout with # prefix showing which contigs have telomeres and where
void log_contig_telomeres(const PathHandleGraph* graph,
                          const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals,
                          double threshold=0.8);
