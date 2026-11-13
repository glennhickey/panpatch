#include <unordered_set>
#include <cassert>
#include <zlib.h>
#include <map>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include "panpatch.hpp"
#include "handlegraph/path_position_handle_graph.hpp"
#include <limits>
#include <cstdint>

//#define debug
//#define ultra_debug

// the patch report accumulated during a run (declared in panpatch.hpp); main stamps chrom/hap,
// finalizes the accept/reject decision, prints it as a table, and clears it after each haplotype.
vector<PatchRecord> g_patch_records;

// total sequence length (bp) of a path
int64_t path_bp(const PathHandleGraph* g, path_handle_t p) {
    int64_t n = 0;
    g->for_each_step_in_path(p, [&](step_handle_t s) { n += g->get_length(g->get_handle_of_step(s)); });
    return n;
}

void sort_paths_by_name(const PathHandleGraph* graph, vector<path_handle_t>& paths) {
    // path names are unique, so this is a total order
    vector<pair<string, path_handle_t>> named;
    named.reserve(paths.size());
    for (const path_handle_t& p : paths) named.emplace_back(graph->get_path_name(p), p);
    sort(named.begin(), named.end(), [](const pair<string, path_handle_t>& a, const pair<string, path_handle_t>& b) {
        return a.first < b.first;
    });
    for (size_t i = 0; i < named.size(); ++i) paths[i] = named[i].second;
}

vector<path_handle_t> sample_paths(const PathHandleGraph* graph, const string& sample) {
    vector<path_handle_t> paths;
    graph->for_each_path_of_sample(sample, [&](path_handle_t p) { paths.push_back(p); });
    sort_paths_by_name(graph, paths);
    return paths;
}

vector<path_handle_t> default_sample_paths(const PathHandleGraph* graph, const string& default_sample,
                                           int64_t target_haplotype) {
    vector<path_handle_t> paths;
    for (const path_handle_t& p : sample_paths(graph, default_sample)) {
        size_t hap = graph->get_haplotype(p);
        if (hap == 0 || hap == PathMetadata::NO_HAPLOTYPE || hap == (size_t)target_haplotype) {
            paths.push_back(p);
        }
    }
    return paths;
}

// offset of a step along its path.  only used to put the steps of one path that visits a node more
// than once in a stable order, so the slow walk (for a graph without positions) is rarely taken
static int64_t step_offset(const PathHandleGraph* graph, const step_handle_t& step) {
    if (const PathPositionHandleGraph* pos_graph = dynamic_cast<const PathPositionHandleGraph*>(graph)) {
        return pos_graph->get_position_of_step(step);
    }
    int64_t offset = 0;
    path_handle_t path = graph->get_path_handle_of_step(step);
    for (step_handle_t s = graph->path_begin(path); s != step; s = graph->get_next_step(s)) {
        offset += graph->get_length(graph->get_handle_of_step(s));
    }
    return offset;
}

ErrorRegions error_bed_to_paths(const PathHandleGraph* graph, const BedRegions& error_bed) {
    ErrorRegions regions;
    if (error_bed.empty()) return regions;
    graph->for_each_path_handle([&](path_handle_t p) {
        // the full path name, else that name without its subrange (a BED naming SAMPLE#HAP#CONTIG#0 for a
        // path SAMPLE#HAP#CONTIG#0[a-b]), else its PanSN sample#haplotype#locus (which leaves out the phase
        // block and subrange the graph may add: cactus names a donor contig SAMPLE#HAP#CONTIG#0, or CONTIG[a-b])
        auto it = error_bed.find(graph->get_path_name(p));
        if (it == error_bed.end() && graph->get_subrange(p) != PathMetadata::NO_SUBRANGE) {
            it = error_bed.find(PathMetadata::create_path_name(graph->get_sense(p), graph->get_sample_name(p),
                                                               graph->get_locus_name(p), graph->get_haplotype(p),
                                                               graph->get_phase_block(p), PathMetadata::NO_SUBRANGE));
        }
        if (it == error_bed.end() && graph->get_sample_name(p) != PathMetadata::NO_SAMPLE_NAME) {
            size_t hap = graph->get_haplotype(p);
            string key = graph->get_sample_name(p) + "#" +
                (hap == PathMetadata::NO_HAPLOTYPE ? string() : to_string(hap) + "#") + graph->get_locus_name(p);
            it = error_bed.find(key);
        }
        if (it == error_bed.end()) return;
        int64_t offset = 0;
        if (graph->get_subrange(p) != PathMetadata::NO_SUBRANGE) offset = graph->get_subrange(p).first;
        int64_t len = path_bp(graph, p);
        vector<pair<int64_t, int64_t>> ivs;
        for (const auto& iv : it->second) {
            int64_t s = max((int64_t)0, iv.first - offset), e = min(len, iv.second - offset);
            if (s < e) ivs.push_back(make_pair(s, e));
        }
        if (ivs.empty()) return;
        sort(ivs.begin(), ivs.end());
        vector<pair<int64_t, int64_t>>& merged = regions[p];
        for (const auto& iv : ivs) {
            if (!merged.empty() && iv.first <= merged.back().second) merged.back().second = max(merged.back().second, iv.second);
            else merged.push_back(iv);
        }
    });
    return regions;
}

bool step_in_error(const PathHandleGraph* graph, const step_handle_t& step, const ErrorRegions& error_regions) {
    if (error_regions.empty()) return false;
    auto it = error_regions.find(graph->get_path_handle_of_step(step));
    if (it == error_regions.end()) return false;
    int64_t s = step_offset(graph, step);
    int64_t e = s + graph->get_length(graph->get_handle_of_step(step));
    const vector<pair<int64_t, int64_t>>& ivs = it->second;
    // the first interval that ends after the node starts
    auto j = upper_bound(ivs.begin(), ivs.end(), s,
                         [](int64_t v, const pair<int64_t, int64_t>& iv) { return v < iv.second; });
    return j != ivs.end() && j->first < e;
}

// the distinct paths of a list of intervals, in the order they first appear (ie output order).  used
// to walk the per-contig maps below, which are unordered, in a stable order
static vector<path_handle_t> paths_in_order(const PathHandleGraph* graph,
                                            const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals) {
    vector<path_handle_t> order;
    unordered_set<path_handle_t> seen;
    for (const auto& iv : intervals) {
        path_handle_t p = graph->get_path_handle_of_step(get<0>(iv));
        if (seen.insert(p).second) order.push_back(p);
    }
    return order;
}

unordered_map<path_handle_t, double> compute_overlap_identity(const PathHandleGraph* graph,
                                                              const vector<path_handle_t>& tgt_paths,
                                                              const vector<path_handle_t>& other_paths,
                                                              int64_t w) {
    unordered_map<path_handle_t, int64_t> total_window_overlaps;
    unordered_map<path_handle_t, int64_t> total_window_counts;
    unordered_map<path_handle_t, int64_t> window_overlaps;

    unordered_set<path_handle_t> other_path_set(other_paths.begin(), other_paths.end());

    for (path_handle_t tgt_path : tgt_paths) {
        int64_t start_offset = 0;
        int64_t cur_win_length = 0;
        for (step_handle_t step = graph->path_begin(tgt_path);
             step != graph->path_end(tgt_path);) {
            handle_t handle = graph->get_handle_of_step(step);
            int64_t len = graph->get_length(handle);
#ifdef ultra_debug
            cerr << " node " << graph->get_id(handle) << " has len " << len << endl;        
            cerr << "starting step at " << graph->get_id(handle) << " offset " << start_offset << endl;
#endif
            // scan forward "w" bases to find last position of our window
            step_handle_t next_step;
            int64_t next_offset;
            for (next_step = step;
                 next_step != graph->path_end(tgt_path);
                 next_step = graph->get_next_step(next_step)) {
                handle_t next_handle = graph->get_handle_of_step(next_step);
                int64_t next_len = graph->get_length(next_handle);

                int64_t room = step == next_step ? next_len - start_offset : next_len;
                if (cur_win_length + room >= w) {
                    // we can cover our window by ending in this node
                    next_offset = w - cur_win_length;
#ifdef ultra_debug
                    cerr << "room " << room << endl;
                    cerr << "next offset = " << w << " - " << cur_win_length << " = " << next_offset << endl;
#endif
                    if (step == next_step) {
                        next_offset += start_offset;
                    }
                    cur_win_length = w;
                    assert(next_offset <= next_len);
#ifdef ultra_debug
                    cerr << "cutting window scan at " << graph->get_id(graph->get_handle_of_step(next_step))
                         << " offset " << next_offset << endl;
#endif
                    break;
                } else {
                    // we keep looking
                    cur_win_length += room;
                }
            }

            // re-walk the window, this time counting the number of bases in
            // each other path we cover.
            int64_t tot_cur_len = 0;
            for (step_handle_t cur_step = step;
                 cur_step != graph->path_end(tgt_path);
                 cur_step = graph->get_next_step(cur_step)) {
                int64_t cur_len = graph->get_length(graph->get_handle_of_step(cur_step));
                if (cur_step == step) {
                    cur_len -= start_offset;
                }
                if (cur_step == next_step) {
                    cur_len -= (graph->get_length(graph->get_handle_of_step(cur_step)) - next_offset);
                }
                tot_cur_len += cur_len;

                handle_t handle = graph->get_handle_of_step(cur_step);
                // todo: cycles handled fairly naively, but don't think at the resolution
                // we're looking at really matters.
                unordered_set<path_handle_t> other_paths;
                graph->for_each_step_on_handle(handle, [&](step_handle_t other_step) {
                    path_handle_t other_path = graph->get_path_handle_of_step(other_step);
                    if (other_path_set.count(other_path)) {
                        other_paths.insert(other_path);
                    }
                });

                for (path_handle_t other_path : other_paths) {
                    window_overlaps[other_path] += cur_len;
                }
                if (cur_step == next_step) {
                    break;
                }
            }
            assert(tot_cur_len == cur_win_length);
            // add it to the total
            for (const auto& path_count : window_overlaps) {
                assert(path_count.second <= cur_win_length);
                total_window_overlaps[path_count.first] += path_count.second;
                total_window_counts[path_count.first] += cur_win_length;
            
            }

            // move to next window
            step = next_step;
            start_offset = next_offset;
            window_overlaps.clear();
            cur_win_length = 0;
        }
    }
    
    unordered_map<path_handle_t, double> coverage_map;
    for (const auto& path_overlaps : total_window_overlaps) {
        coverage_map[path_overlaps.first] = (double)path_overlaps.second / (double)total_window_counts[path_overlaps.first];
#ifdef debug
        cerr << " COV MAP " << graph->get_path_name(path_overlaps.first) << " -> " << coverage_map[path_overlaps.first] << endl;
#endif
    }
    
    return coverage_map;
}

unordered_map<string, vector<path_handle_t>> select_sample_covers(const PathHandleGraph* graph,
                                                                  const unordered_map<path_handle_t, double>& coverage_map) {

    // visit the paths by name, not in hash order: the floating-point sums below then always add up
    // the same way, and each sample's selected paths come out sorted
    vector<path_handle_t> paths;
    for (const auto& path_cov : coverage_map) {
        paths.push_back(path_cov.first);
    }
    sort_paths_by_name(graph, paths);

    // haplotypes are ordered, so a tie in mean coverage goes to the lowest haplotype
    unordered_map<string, map<int64_t, double>> sample_hap_coverage;
    unordered_map<string, map<int64_t, int64_t>> sample_hap_count;
    for (const path_handle_t& path : paths) {
        string sample = graph->get_sample_name(path);
        int64_t haplotype = graph->get_haplotype(path);
        int64_t path_length = 0;
        graph->for_each_step_in_path(path, [&](const step_handle_t& step) {
            path_length += graph->get_length(graph->get_handle_of_step(step));
        });
        sample_hap_coverage[sample][haplotype] += (double)path_length * coverage_map.at(path);
        sample_hap_count[sample][haplotype] += path_length;
    }

    // todo: this is prety unsophistacated and could be tricked by some weird edge cases
    // but should be fine for very high-quality assemblies (like we deal with
    unordered_map<string, int64_t> sample_to_hap;
    for (const auto& sample_cov : sample_hap_coverage) {
        int64_t best_haplotype = -1;
        double best_mean_cov = 0;
        for (const auto& hap_cov : sample_cov.second) {
            double mean_cov = (double)hap_cov.second / (double)sample_hap_count[sample_cov.first][hap_cov.first];
            if (mean_cov > best_mean_cov) {
                best_haplotype = hap_cov.first;
                best_mean_cov = mean_cov;
            }
        }
        sample_to_hap[sample_cov.first] = best_haplotype;
    }

    unordered_map<string, vector<path_handle_t>> result;
    for (const path_handle_t& path : paths) {
        string sample = graph->get_sample_name(path);
        int64_t haplotype = graph->get_haplotype(path);
        if (haplotype == sample_to_hap[sample]) {
            result[sample].push_back(path);
        }
    }

    return result;
}

BedRegions parse_bed_file(const string& bed_filename) {
    BedRegions regions;
    ifstream bed_file(bed_filename);
    if (!bed_file) {
        cerr << "[panpatch] error: Unable to open BED file " << bed_filename << endl;
        exit(1);
    }
    string line;
    int64_t line_num = 0;
    while (getline(bed_file, line)) {
        ++line_num;
        if (line.empty() || line[0] == '#' || line.substr(0, 5) == "track" || line.substr(0, 7) == "browser") {
            continue;
        }
        istringstream ss(line);
        string contig;
        int64_t start, end;
        if (ss >> contig >> start >> end) {
            if (start < 0 || end < 0 || start >= end) {
                cerr << "[panpatch] warning: skipping malformed BED line " << line_num << ": " << line << endl;
                continue;
            }
            regions[contig].push_back(make_pair(start, end));
        } else {
            cerr << "[panpatch] warning: skipping malformed BED line " << line_num << ": " << line << endl;
        }
    }
    return regions;
}

static bool interval_overlaps_excluded(int64_t start, int64_t end, const ExcludedRefRegions& excluded) {
    if (excluded.empty()) return false;
    // find first excluded region where region.end > start
    auto it = upper_bound(excluded.begin(), excluded.end(), make_pair(start, start),
                          [](const pair<int64_t, int64_t>& a, const pair<int64_t, int64_t>& b) {
                              return a.second < b.second;
                          });
    // check if this region overlaps [start, end)
    if (it != excluded.end() && it->first < end) {
        return true;
    }
    return false;
}

ExcludedRefRegions bed_to_ref_regions(const PathHandleGraph* graph,
                                      const vector<path_handle_t>& tgt_paths,
                                      const unordered_map<int64_t, int64_t>& ref_anchors,
                                      const BedRegions& bed_regions) {
    ExcludedRefRegions excluded;

    for (const path_handle_t& tgt_path : tgt_paths) {
        string path_name = graph->get_path_name(tgt_path);

        // check if this path has any BED regions
        if (!bed_regions.count(path_name)) {
            continue;
        }
        const vector<pair<int64_t, int64_t>>& bed_intervals = bed_regions.at(path_name);

        // walk the target path, building ordered map of target_pos -> ref_pos for anchor nodes
        map<int64_t, int64_t> tgt_pos_to_ref_pos;
        int64_t tgt_pos = 0;
        graph->for_each_step_in_path(tgt_path, [&](step_handle_t step) {
            handle_t handle = graph->get_handle_of_step(step);
            int64_t node_id = graph->get_id(handle);
            if (ref_anchors.count(node_id)) {
                tgt_pos_to_ref_pos[tgt_pos] = ref_anchors.at(node_id);
            }
            tgt_pos += graph->get_length(handle);
        });

        if (tgt_pos_to_ref_pos.empty()) {
            cerr << "[panpatch] warning: no anchor mapping found for BED contig " << path_name << endl;
            continue;
        }

        // for each BED interval, find the corresponding ref-position range
        for (const auto& bed_interval : bed_intervals) {
            int64_t bed_start = bed_interval.first;
            int64_t bed_end = bed_interval.second;

            // find anchors that overlap the BED region [bed_start, bed_end)
            // lower_bound: first anchor at or after bed_start
            auto it_start = tgt_pos_to_ref_pos.lower_bound(bed_start);
            // we also want the anchor just before bed_start if it exists
            if (it_start != tgt_pos_to_ref_pos.begin()) {
                auto prev = std::prev(it_start);
                // include previous anchor if the node it represents could extend into the BED region
                it_start = prev;
            }
            // upper_bound: first anchor strictly after bed_end
            auto it_end = tgt_pos_to_ref_pos.lower_bound(bed_end);

            if (it_start == tgt_pos_to_ref_pos.end()) {
                continue;
            }

            // collect all ref positions in this range and take min/max
            // (handles reversed paths where ref positions may not be monotonic)
            int64_t ref_min = numeric_limits<int64_t>::max();
            int64_t ref_max = numeric_limits<int64_t>::min();
            for (auto it = it_start; it != it_end; ++it) {
                ref_min = min(ref_min, it->second);
                ref_max = max(ref_max, it->second);
            }
            if (ref_min <= ref_max) {
                excluded.push_back(make_pair(ref_min, ref_max));
            }
        }
    }

    // sort and merge overlapping regions
    sort(excluded.begin(), excluded.end());
    ExcludedRefRegions merged;
    for (const auto& region : excluded) {
        if (!merged.empty() && region.first <= merged.back().second) {
            merged.back().second = max(merged.back().second, region.second);
        } else {
            merged.push_back(region);
        }
    }
    return merged;
}

unordered_map<int64_t, int64_t> find_anchors(const PathHandleGraph* graph,
                                             const path_handle_t& ref_path,
                                             const vector<path_handle_t>& tgt_paths,
                                             const unordered_set<path_handle_t>& path_set) {
    unordered_map<int64_t, int64_t> anchors;
    int64_t pos = 0;
    
    pair<int64_t, int64_t> prev_anchor;
    pair<int64_t, int64_t> prev_prev_anchor;
    unordered_set<path_handle_t> prev_anchor_paths;
    unordered_set<path_handle_t> prev_prev_anchor_paths;
    bool ref_is_target = false;
    for (const path_handle_t& tgt_path : tgt_paths) {
        if (tgt_path == ref_path) {
            ref_is_target = true;
            break;
        }
    }

    // todo: we can optimise this a lot by targeting edges of tgt paths, assembly gaps etc...
    // but for now we just return everything
    graph->for_each_step_in_path(ref_path, [&](step_handle_t step) {
        handle_t handle = graph->get_handle_of_step(step);
        unordered_set<path_handle_t> covered_paths;
        graph->for_each_step_on_handle(handle, [&](step_handle_t other_step) {
            path_handle_t other_path = graph->get_path_handle_of_step(other_step);
            if (path_set.count(other_path) && other_path != ref_path) {
                covered_paths.insert(graph->get_path_handle_of_step(other_step));
            }
        });
        // we currently have no use for ref-only anchors unless they are endpoints
        if (!covered_paths.empty() || step == graph->path_begin(ref_path) || step == graph->path_back(ref_path)) {
            // we currently restrict anchors to target paths.
            // todo: this prevents nested patches (but supporting nesting requires more logic than just removing this check)
            for (const path_handle_t& tgt_path : tgt_paths) {
                if (covered_paths.count(tgt_path) || ref_is_target) {
                    anchors.insert(make_pair(graph->get_id(handle), pos));
                    break;
                }
            }
        }
        pos += graph->get_length(handle);
    });
    return anchors;
}

pair<step_handle_t, bool> find_next_anchor_on_path(const PathHandleGraph* graph,
                                                   const unordered_map<int64_t, int64_t>& anchors,
                                                   step_handle_t step,
                                                   int64_t pos,
                                                   int64_t direction,
                                                   bool cross_Ns,
                                                   int64_t bound,
                                                   const ErrorRegions* error_regions,
                                                   const unordered_map<int64_t, bool>* anchor_rev,
                                                   bool strict) {

    // NOTE:
    // As currently implemented, anchors are all on the reference genome.  This means
    // we can take simplifying steps by making sure we're going in the right direction
    // by looking at the reference positions.

    // Since rearrangements may break the order of our target genome along the reference,
    // we let it go "bound" anchors in the wrong direction before giving up search
    //
    // Relaxing this assumption (which should be done) would require moving to
    // a slightly less trivial graph search. 
    path_handle_t path = graph->get_path_handle_of_step(step);

    function<bool(handle_t)> has_n = [&](handle_t handle) {
        string s = graph->get_sequence(handle);
        for (char c : s) {
            if (c == 'n' || c == 'N') {
#ifdef debug
                cerr << "aborting N (pos=" << pos << ")" << endl;
#endif
                return true;
            }
        }
        return false;
    };
    // a step the interval may not use: an N in its node (an assembly gap), or an --error-bed interval on
    // its own path.  the second is per path, not per node: where the target's flagged sequence is
    // aligned to a donor's clean copy of it, the shared node is blocked for the target and open to the donor
    auto blocked = [&](const step_handle_t& s) {
        return has_n(graph->get_handle_of_step(s)) || (error_regions && step_in_error(graph, s, *error_regions));
    };
    // strict (a donor): an anchor that is not ahead of pos is sequence the patch already has (or would skip
    // back over), and one the path traverses the opposite way to the reference is an inversion between
    // them; a donor run through either would splice a rearrangement into the patch, so it stops there.
    // (k-mer recovery cannot see an inversion: it counts canonical k-mers)
    auto strict_stop = [&](const handle_t& h, bool backward) {
        nid_t id = graph->get_id(h);
        if (anchors.at(id) <= pos) return true;
        if (anchor_rev) {
            auto it = anchor_rev->find(id);
            bool traversed_rev = backward ? !graph->get_is_reverse(h) : graph->get_is_reverse(h);
            if (it != anchor_rev->end() && traversed_rev != it->second) return true;
        }
        return false;
    };

    // the interval starts on this step's node, so it must be N-free too: where the target's and a
    // donor's N runs (gaps) were aligned to each other, the shared N node can be an anchor, and a donor
    // interval starting there would carry its N's into the patch.  an --error-bed interval does not
    // block the start: the anchor is a node this path shares with another, so its sequence is the same
    // whichever of them is threaded through it, and only the steps after it decide what goes in.  this
    // is what lets the target take over again on the last shared node of its flagged region, so a swap
    // stops there instead of also taking the donor's allele of the next, unflagged, bubble
    if (!cross_Ns && has_n(graph->get_handle_of_step(step))) {
        return make_pair(graph->path_end(path), false);
    }

    // search forward
    int64_t fcount = 0;
    int64_t bcount = 0;
    if (direction >= 0) {
        step_handle_t next_step = graph->get_next_step(step);
        for (; next_step != graph->path_end(path); next_step = graph->get_next_step(next_step)) {
            handle_t next_handle = graph->get_handle_of_step(next_step);
            ++fcount;
            if (!cross_Ns && blocked(next_step)) {
                break;
            }
            if (anchors.count(graph->get_id(next_handle))) {
                if (strict && strict_stop(next_handle, false)) {
                    break;
                }
                if (anchors.at(graph->get_id(next_handle)) > pos) {
#ifdef ultra_debug
                    cerr << "Find anchor INPUT=" << graph->get_id(graph->get_handle_of_step(step)) << ":"
                         << graph->get_is_reverse(graph->get_handle_of_step(step)) << " pos=" << pos
                         << " path=" << graph->get_path_name(graph->get_path_handle_of_step(step))
                         << " OUTPUT=" << graph->get_id(graph->get_handle_of_step(next_step)) << ":"
                         << graph->get_is_reverse(graph->get_handle_of_step(next_step)) << " pos="
                         << anchors.at(graph->get_id(next_handle)) << endl;
#endif
                    return make_pair(next_step, false);
                } else if (fcount > bound) {
                    break;
                }
            }
        }
    }
    // search backward
    if (direction <= 0) {
        step_handle_t next_step = graph->get_previous_step(step);
        for (; next_step != graph->path_front_end(path); next_step = graph->get_previous_step(next_step)) {
            handle_t next_handle = graph->get_handle_of_step(next_step);
            ++bcount;
            if (!cross_Ns && blocked(next_step)) {
                break;
            }
            if (anchors.count(graph->get_id(next_handle))) {
                if (strict && strict_stop(next_handle, true)) {
                    break;
                }
                if (anchors.at(graph->get_id(next_handle)) > pos) {
#ifdef ultra_debug
                    cerr << "Find anchor INPUT=" << graph->get_id(graph->get_handle_of_step(step)) << ":"
                         << graph->get_is_reverse(graph->get_handle_of_step(step)) << " pos=" << pos
                         << " path=" << graph->get_path_name(graph->get_path_handle_of_step(step))
                         << " OUTPUT=" << graph->get_id(graph->get_handle_of_step(next_step)) << ":"
                         << graph->get_is_reverse(graph->get_handle_of_step(next_step)) << " pos="
                         << anchors.at(graph->get_id(next_handle)) << endl;
#endif
                    return make_pair(next_step, true);
                } else if (bcount > bound) {
                    break;
                }
            }
        }
    }
#ifdef ultra_debug
    cerr << "NA fail fc="  << fcount << " bc=" << bcount << " " << graph->get_path_name(path) << endl;
#endif
    return make_pair(graph->path_end(path), false);
}


vector<tuple<step_handle_t, step_handle_t, bool>> thread_intervals(const PathHandleGraph* graph,
                                                                   const path_handle_t& ref_path,
                                                                   const unordered_map<int64_t, int64_t>& ref_anchors,
                                                                   const vector<path_handle_t>& tgt_paths,
                                                                   const vector<path_handle_t>& other_paths,
                                                                   const ExcludedRefRegions& excluded_regions,
                                                                   const ErrorRegions* error_regions) {

    // fewer than two reference anchors: the target shares <=1 node with the reference path (common
    // for acrocentric short arms aligned against an external reference), so there is nothing to
    // thread between -- return empty and let greedy_patch revert this contig to its input sequence
    if (ref_anchors.size() <= 1 || tgt_paths.empty()) {
        cerr << "[panpatch] warning: only " << ref_anchors.size() << " reference anchor(s) on "
             << graph->get_path_name(ref_path) << "; emitting input unchanged" << endl;
        return {};
    }

    // build target path set for O(1) lookup when checking excluded regions
    unordered_set<path_handle_t> tgt_path_set(tgt_paths.begin(), tgt_paths.end());

    unordered_map<path_handle_t, int64_t> path_rank;
    for (int64_t i = 0; i < tgt_paths.size(); ++i) {
        path_rank[tgt_paths[i]] = i;
    }
    // important note: these must be already sorted in order of sample priortiy    
    for (int64_t i = 0; i < other_paths.size(); ++i) {
        path_rank[other_paths[i]] = i + tgt_paths.size();
    }

    vector<tuple<step_handle_t, step_handle_t, bool>> interval_cover;

    // find the first anchor (note, it requires a full scan since we don't have
    // it indexed anywhere)
    int64_t cur_pos = numeric_limits<int64_t>::max();
    handle_t cur_handle;
    graph->for_each_step_in_path(ref_path, [&](step_handle_t step) {
        handle_t handle = graph->get_handle_of_step(step);        
        if (ref_anchors.count(graph->get_id(handle))) {
            int64_t handle_pos = ref_anchors.at(graph->get_id(handle));
            if (handle_pos < cur_pos) {
                assert(cur_pos == numeric_limits<int64_t>::max());
                cur_pos = handle_pos;
                cur_handle = handle;
            }
        }
    });
    // the orientation the reference traverses each anchor's node in (on its first visit, the one
    // find_anchors positions it by): a donor is followed only in the direction that traverses the anchors
    // the same way (the strand rule below)
    unordered_map<int64_t, bool> anchor_rev;
    graph->for_each_step_in_path(ref_path, [&](step_handle_t step) {
        handle_t handle = graph->get_handle_of_step(step);
        if (ref_anchors.count(graph->get_id(handle)) && !anchor_rev.count(graph->get_id(handle))) {
            anchor_rev[graph->get_id(handle)] = graph->get_is_reverse(handle);
        }
    });
    bool cur_backward = graph->get_is_reverse(cur_handle);   // (the reference's own first anchor may be traversed in reverse)
    unordered_map<path_handle_t, int64_t> tgt_ranks;
    // the step we reached the current anchor on (none for the first anchor)
    bool have_cur_step = false;
    step_handle_t cur_step;
        
    while (true) {
        // find the steps on the handle and sort them using the path priority
        // [paths that don't have a priority are ignored which is important]
        // todo: could be in weird regions we need to assess contiguiuty with
        // previous anchor, but not trying for now. 
        vector<pair<int64_t, step_handle_t>> sorted_steps;
        graph->for_each_step_on_handle(cur_handle, [&](step_handle_t step) {
            path_handle_t path = graph->get_path_handle_of_step(step);
            if (path_rank.count(path)) {
                int64_t rank = tgt_ranks.count(path) ? tgt_ranks.at(path) : path_rank.at(path);
                sorted_steps.push_back(make_pair(rank, step));
            }
        });
        // the graph lists the steps on a node in no particular order, so break ties in rank here:
        // by path name, then (for a path that visits the node more than once) the step we arrived
        // on first, so we carry on along the path from where we are, then the rest by path offset
        sort(sorted_steps.begin(), sorted_steps.end(),
             [&](const pair<int64_t, step_handle_t>& a, const pair<int64_t, step_handle_t>& b) {
            if (a.first != b.first) {
                return a.first < b.first;
            }
            path_handle_t a_path = graph->get_path_handle_of_step(a.second);
            path_handle_t b_path = graph->get_path_handle_of_step(b.second);
            if (a_path != b_path) {
                return graph->get_path_name(a_path) < graph->get_path_name(b_path);
            }
            if (a.second == b.second) {
                return false;
            }
            if (have_cur_step && (a.second == cur_step || b.second == cur_step)) {
                return a.second == cur_step;
            }
            return step_offset(graph, a.second) < step_offset(graph, b.second);
        });

        // how many times each path visits this node
        unordered_map<path_handle_t, int64_t> visits;
        for (const auto& rank_step : sorted_steps) {
            ++visits[graph->get_path_handle_of_step(rank_step.second)];
        }

        bool found_next = false;
        for (int64_t iteration = 0; iteration < 2 && !found_next; ++iteration) {
            for (const auto& rank_step : sorted_steps) {
                const step_handle_t& step = rank_step.second;
#ifdef debug
                cerr << "i=" << iteration << " cur_pos " << cur_pos << " cur step "
                     << graph->get_path_name(graph->get_path_handle_of_step(step))
                     << " cur_back " << cur_backward << " " << graph->get_id(graph->get_handle_of_step(step)) << ":"
                     << graph->get_is_reverse(graph->get_handle_of_step(step)) << endl;
#endif
                // the second pass, when no path reaches the next anchor without crossing an N or an
                // --error-bed interval, is for the target alone: it keeps its own gap (or its own flagged
                // sequence) rather than stopping here.  a donor crossing its N's or its own flagged
                // sequence would put them in the patch
                bool cross_gaps = iteration > 0;
                if (cross_gaps && !tgt_path_set.count(graph->get_path_handle_of_step(step))) {
                    continue;
                }
                // the path we are following goes on only from the step we reached this anchor on.  a path
                // can visit the anchor's node more than once (a donor through a satellite array): taking up
                // another of its visits here would put two pieces of it side by side that are not contiguous
                // on it, which smooth_intervals then merges into one interval running the wrong way (and
                // even kept apart, they would splice the donor's sequence out of order)
                path_handle_t step_path = graph->get_path_handle_of_step(step);
                bool following = have_cur_step && step_path == graph->get_path_handle_of_step(cur_step);
                if (following && step != cur_step) {
                    continue;
                }
                // nor is a path taken up here, from another, if it visits this node more than once: which of
                // its visits the patch should go on from is ambiguous, and the wrong one drops or repeats
                // sequence (a donor with a stray copy of the node goes back over what is already patched)
                if (!following && visits[step_path] > 1) {
                    continue;
                }
                bool step_backward = graph->get_is_reverse(graph->get_handle_of_step(step));
                bool is_tgt = tgt_path_set.count(step_path);
                int64_t direction = 0;
                if (!is_tgt) {
                    // the strand rule: a donor is followed in the direction that traverses this node the way
                    // the reference does.  (searching both ways took whichever reached an anchor first, so a
                    // donor carrying an inversion or a rotation of the target spliced it in backwards)
                    auto ar = anchor_rev.find(graph->get_id(cur_handle));
                    direction = (ar == anchor_rev.end() || step_backward == ar->second) ? 1 : -1;
                } else if (!interval_cover.empty() && graph->get_path_handle_of_step(get<0>(interval_cover.back())) ==
                           step_path) {
                    // we constrain the direction if we're already following the path, otherwise it can go either way
                    direction = get<2>(interval_cover.back()) ? -1 : 1;
                }
                pair<step_handle_t, bool> next_anchor = find_next_anchor_on_path(graph, ref_anchors, step, cur_pos,
                                                                                 direction, cross_gaps, 10000,
                                                                                 error_regions, &anchor_rev, !is_tgt);
                if (next_anchor.first != graph->path_end(graph->get_path_handle_of_step(step))) {
                    // check if this interval overlaps an excluded region for non-target paths
                    path_handle_t step_path = graph->get_path_handle_of_step(step);
                    if (!excluded_regions.empty() && !tgt_path_set.count(step_path)) {
                        int64_t next_pos = ref_anchors.at(graph->get_id(graph->get_handle_of_step(next_anchor.first)));
                        if (interval_overlaps_excluded(min(cur_pos, next_pos), max(cur_pos, next_pos), excluded_regions)) {
#ifdef debug
                            cerr << "Skipping interval on " << graph->get_path_name(step_path)
                                 << " at ref pos " << cur_pos << "-" << next_pos
                                 << " due to excluded region" << endl;
#endif
                            continue;
                        }
                    }
                    interval_cover.push_back(make_tuple(step, next_anchor.first, next_anchor.second));
#ifdef debug
                    const auto& interval = interval_cover.back();
                    cerr << "Adding interval " << graph->get_path_name(graph->get_path_handle_of_step(get<0>(interval))) << " " 
                         << graph->get_id(graph->get_handle_of_step(get<0>(interval))) << ":"
                         << graph->get_is_reverse(graph->get_handle_of_step(get<0>(interval))) << " - "
                         << graph->get_id(graph->get_handle_of_step(get<1>(interval))) << ":"
                         << graph->get_is_reverse(graph->get_handle_of_step(get<1>(interval)))
                         << " rev=" <<get<2>(interval) << " range=" << cur_pos << "-" << ref_anchors.at(graph->get_id(cur_handle))
                         << endl;
                    check_intervals(graph, {interval});
#endif
                    found_next = true;
                    // slide position to the next anchor
                    have_cur_step = true;
                    cur_step = next_anchor.first;
                    cur_handle = graph->get_handle_of_step(next_anchor.first);
                    cur_pos = ref_anchors.at(graph->get_id(cur_handle));
                    cur_backward = graph->get_is_reverse(cur_handle);
                    path_handle_t next_path = graph->get_path_handle_of_step(next_anchor.first);
                    if (path_rank[next_path] < (int64_t)tgt_paths.size() && !tgt_ranks.count(next_path)) {
                        // preference for first-visited target paths when choosing fork: ranked ahead of every
                        // path_rank (which starts at 0), in the order they were first visited.  (this used to
                        // compare against tgt_ranks.size(), 0 at first, so no path was ever ranked)
                        tgt_ranks[next_path] = (int64_t)tgt_ranks.size() - (int64_t)tgt_paths.size() - 1;
                    }
                    break;
                }
            }
        }
        if (!found_next) {
            break;
        }
    }

    return interval_cover;
}

vector<tuple<step_handle_t, step_handle_t, bool>> smooth_intervals(const PathHandleGraph* graph,
                                                                   const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals) {
    vector<tuple<step_handle_t, step_handle_t, bool>> smoothed_intervals;
    int64_t i = 0;
    for (int64_t j = 1; j < intervals.size(); ++j) {
        path_handle_t i_path = graph->get_path_handle_of_step(get<0>(intervals[i]));
        path_handle_t j_path = graph->get_path_handle_of_step(get<0>(intervals[j]));
        if (i_path == j_path && get<2>(intervals[i]) == get<2>(intervals[j])) {
            // merge interval i with j, just by doing nothing
            continue;
        } else {
#ifdef debug
            cerr << "i: " << graph->get_id(graph->get_handle_of_step(get<0>(intervals[i]))) <<"-"
                 << graph->get_id(graph->get_handle_of_step(get<1>(intervals[i]))) << ":"
                 << get<2>(intervals[i]) << endl;
            cerr << "j: " << graph->get_id(graph->get_handle_of_step(get<0>(intervals[j]))) <<"-"
                 << graph->get_id(graph->get_handle_of_step(get<1>(intervals[j]))) << ":"
                 << get<2>(intervals[j]) << endl;
#endif
            // j cannot be merged with i, so we write interval from i to j-1 inclusive
            smoothed_intervals.push_back(make_tuple(get<0>(intervals[i]), get<1>(intervals[j-1]), get<2>(intervals[i])));
            i = j;
        }
    }
    // add last interval
    int64_t j = intervals.size();
    smoothed_intervals.push_back(make_tuple(get<0>(intervals[i]), get<1>(intervals[j-1]), get<2>(intervals[i])));

#ifdef debug
    for (const auto& interval : smoothed_intervals) {
        cerr << "Smoothed interval " << graph->get_path_name(graph->get_path_handle_of_step(get<0>(interval))) << " "
             << graph->get_id(graph->get_handle_of_step(get<0>(interval))) << ":"
             << graph->get_is_reverse(graph->get_handle_of_step(get<0>(interval))) << " - "
             << graph->get_id(graph->get_handle_of_step(get<1>(interval))) << ":"
             << graph->get_is_reverse(graph->get_handle_of_step(get<1>(interval)))
             << " rev=" <<get<2>(interval) << endl;
    }
#endif

    return smoothed_intervals;
}

vector<tuple<step_handle_t, step_handle_t, bool>> extend_intervals(const PathHandleGraph* graph,
                                                                   const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals,
                                                                   const vector<path_handle_t>& tgt_paths,
                                                                   const ErrorRegions* error_regions) {
    vector<tuple<step_handle_t, step_handle_t, bool>> extended_intervals = intervals;
    unordered_set<path_handle_t> tgt_path_set(tgt_paths.begin(), tgt_paths.end());

    // a terminal interval on the target is extended to the end of its path.  one on a donor (where
    // threading ended on a donor) is extended only up to its first N or --error-bed interval: past that
    // is an assembly gap or sequence flagged as an error, which would otherwise come into the patch with
    // the donor's tail.  a patch that comes up short of a target telomere this way is reverted by
    // discards_target_telomere
    auto node_has_n = [&](const step_handle_t& s) {
        string seq = graph->get_sequence(graph->get_handle_of_step(s));
        return seq.find_first_of("Nn") != string::npos || (error_regions && step_in_error(graph, s, *error_regions));
    };
    // walk from s toward the start (forward=false) or end of its path, as far as the path (target) or
    // its first N (donor) allows
    auto extend_to = [&](step_handle_t s, bool forward) -> step_handle_t {
        path_handle_t p = graph->get_path_handle_of_step(s);
        if (tgt_path_set.count(p)) {
            return forward ? graph->get_previous_step(graph->path_end(p)) : graph->path_begin(p);
        }
        while (true) {
            step_handle_t nx = forward ? graph->get_next_step(s) : graph->get_previous_step(s);
            if (nx == (forward ? graph->path_end(p) : graph->path_front_end(p)) || node_has_n(nx)) break;
            s = nx;
        }
        return s;
    };

    // extend the front (toward the path's start, or its end for a reverse interval)
    step_handle_t first_step = get<0>(extended_intervals[0]);
    bool first_backward = get<2>(extended_intervals[0]);
    first_step = extend_to(first_step, first_backward);
    extended_intervals[0] = make_tuple(first_step, get<1>(extended_intervals[0]), first_backward);

    // extend the back (toward the path's end, or its start for a reverse interval)
    step_handle_t last_step = get<1>(extended_intervals.back());
    bool last_backward = get<2>(extended_intervals.back());
    last_step = extend_to(last_step, !last_backward);
    extended_intervals.back() = make_tuple(get<0>(extended_intervals.back()), last_step, last_backward);

#ifdef debug
    for (const auto& interval : extended_intervals) {
        cerr << "Extended interval " << graph->get_path_name(graph->get_path_handle_of_step(get<0>(interval))) << " " 
             << graph->get_id(graph->get_handle_of_step(get<0>(interval))) << ":"
             << graph->get_is_reverse(graph->get_handle_of_step(get<0>(interval))) << " - "
             << graph->get_id(graph->get_handle_of_step(get<1>(interval))) << ":"
             << graph->get_is_reverse(graph->get_handle_of_step(get<1>(interval)))
             << " rev=" <<get<2>(interval) << endl;
    }
#endif

    return extended_intervals;
}

// Find where a telomere run anchored at a tip ends/begins within sequence[search_start, search_end).
//   scan_forward : telomere anchored at search_start; returns the index where the run ends (-1 if none).
//  !scan_forward : telomere anchored at search_end;   returns the index where the run begins (-1 if none).
// Slides a 500bp window inward from the tip, extending the run while window density stays >= 0.7.
static int64_t find_telomere_boundary(const string& sequence, int64_t search_start, int64_t search_end,
                                      bool scan_forward) {
    const int64_t window_size = 500;
    const double min_density = 0.7;
    if (scan_forward) {
        int64_t telomere_end = search_start;
        for (int64_t win_start = search_start; win_start + window_size < search_end; win_start += 100) {
            int64_t repeats = 0;
            for (int64_t pos = win_start; pos < min(win_start + window_size, search_end - 6); ++pos) {
                if (sequence.substr(pos, 6) == "TTAGGG" || sequence.substr(pos, 6) == "CCCTAA") { ++repeats; pos += 5; }
            }
            double density = 6.0 * (double)repeats / (double)window_size;
            if (density >= min_density) telomere_end = win_start + window_size;
            else if (telomere_end > search_start) break;
        }
        return telomere_end > search_start ? telomere_end : -1;
    } else {
        int64_t telomere_start = search_end;
        for (int64_t win_end = search_end; win_end - window_size > search_start; win_end -= 100) {
            int64_t win_start = max(search_start, win_end - window_size);
            int64_t repeats = 0;
            for (int64_t pos = win_start; pos < win_end - 6; ++pos) {
                if (sequence.substr(pos, 6) == "TTAGGG" || sequence.substr(pos, 6) == "CCCTAA") { ++repeats; pos += 5; }
            }
            double density = 6.0 * (double)repeats / (double)window_size;
            if (density >= min_density) telomere_start = win_start;
            else if (telomere_start < search_end) break;
        }
        return telomere_start < search_end ? telomere_start : -1;
    }
}

// Canonical "does this region carry a telomere?" test, shared by telomere validation and patching.
// When find_boundary is set, the region is first narrowed to the telomeric run anchored at the tip
// (is_right_end => tip at max_end, else tip at start); the run must then be >=500bp and reach the
// density threshold.  Averaging over the actual run (rather than accepting any single dense window)
// is what makes a telomere buried under terminal junk, or a short/degraded telomere, read as absent.
static bool seq_has_telomere(const string& sequence, int64_t start, int64_t max_end,
                             bool find_boundary, bool is_right_end, double threshold,
                             double* out_density = nullptr) {
    int64_t actual_start = start;
    int64_t actual_end = max_end;
    if (find_boundary) {
        if (is_right_end) {
            int64_t telomere_start = find_telomere_boundary(sequence, start, max_end, false);
            if (telomere_start >= start && telomere_start < max_end) actual_start = telomere_start;
        } else {
            int64_t telomere_end = find_telomere_boundary(sequence, start, max_end, true);
            if (telomere_end > start) actual_end = telomere_end;
        }
    }
    int64_t fw_count = 0, r_count = 0;
    for (int64_t pos = actual_start; pos < actual_end - 6; ++pos) {
        if (sequence.substr(pos, 6) == "TTAGGG") { ++fw_count; pos += 5; }
        else if (sequence.substr(pos, 6) == "CCCTAA") { ++r_count; pos += 5; }
    }
    int64_t region_len = actual_end - actual_start;
    if (region_len < 500) { if (out_density) *out_density = 0.0; return false; }
    double fw_density = 6. * ((double)fw_count / (double)region_len);
    double r_density = 6. * ((double)r_count / (double)region_len);
    if (out_density) *out_density = max(fw_density, r_density);
    return (fw_density >= threshold || r_density >= threshold);
}

// lenient telomere presence: true if ANY 500bp window in s clears the density threshold. Unlike
// seq_has_telomere (which requires a clean terminal telomere), this just detects that telomeric
// repeats exist somewhere in the region - used only to explain *why* a capless tip wasn't patched
// (e.g. a telomere buried under terminal junk, or a degraded/fragmented one).
static bool has_telomeric_window(const string& s, double threshold) {
    const int64_t W = 500;
    int64_t n = (int64_t)s.size();
    for (int64_t i = 0; i + W <= n; i += 100) {
        int64_t fw = 0, rv = 0;
        for (int64_t p = i; p < i + W - 6; ) {
            if (s.compare(p, 6, "TTAGGG") == 0) { ++fw; p += 6; }
            else if (s.compare(p, 6, "CCCTAA") == 0) { ++rv; p += 6; }
            else ++p;
        }
        if (6.0 * (double)max(fw, rv) / (double)W >= threshold) return true;
    }
    return false;
}

// true if s holds a telomeric run of at least min_run bp: consecutive 500bp windows (100bp apart) that each
// clear the density threshold, as has_telomeric_window measures one.  (a single window would also catch the
// short interstitial telomere-like repeats a genome has; an array of a kb or more is a telomere)
static bool has_telomeric_run(const string& s, double threshold, int64_t min_run) {
    const int64_t W = 500, STEP = 100;
    int64_t n = (int64_t)s.size(), run_start = -1;
    for (int64_t i = 0; i + W <= n; i += STEP) {
        int64_t fw = 0, rv = 0;
        for (int64_t p = i; p < i + W - 6; ) {
            if (s.compare(p, 6, "TTAGGG") == 0) { ++fw; p += 6; }
            else if (s.compare(p, 6, "CCCTAA") == 0) { ++rv; p += 6; }
            else ++p;
        }
        if (6.0 * (double)max(fw, rv) / (double)W >= threshold) {
            if (run_start < 0) run_start = i;
            if (i + W - run_start >= min_run) return true;
        } else {
            run_start = -1;
        }
    }
    return false;
}

static inline uint64_t hash64(uint64_t x) {
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33; return x;
}

// insert a ~1/16 subsample of canonical 31-mer hashes of s into out
static void sample_kmers(const string& s, unordered_set<uint64_t>& out) {
    const int K = 31;
    const uint64_t SAMPLE_MASK = 0xF;  // keep hashes with low 4 bits zero -> ~1/16
    int n = (int)s.size();
    if (n < K) return;
    uint64_t fwd = 0, rev = 0, mask = (1ULL << (2 * K)) - 1;
    int valid = 0;
    for (int i = 0; i < n; ++i) {
        int c;
        switch (s[i]) {
            case 'A': case 'a': c = 0; break;
            case 'C': case 'c': c = 1; break;
            case 'G': case 'g': c = 2; break;
            case 'T': case 't': c = 3; break;
            default: c = -1;
        }
        if (c < 0) { valid = 0; fwd = rev = 0; continue; }
        fwd = ((fwd << 2) | (uint64_t)c) & mask;
        rev = (rev >> 2) | ((uint64_t)(3 - c) << (2 * (K - 1)));
        if (++valid >= K) {
            uint64_t h = hash64(min(fwd, rev));
            if ((h & SAMPLE_MASK) == 0) out.insert(h);
        }
    }
}

// the (start position, hash) of each k-mer sample_kmers would keep from s, in order along s
static vector<pair<int64_t, uint64_t>> sampled_kmer_positions(const string& s) {
    const int K = 31;
    const uint64_t SAMPLE_MASK = 0xF;   // the same ~1/16 subsample as sample_kmers
    vector<pair<int64_t, uint64_t>> sampled;
    uint64_t fwd = 0, rev = 0, mask = (1ULL << (2 * K)) - 1;
    int valid = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        int c;
        switch (s[i]) {
            case 'A': case 'a': c = 0; break;
            case 'C': case 'c': c = 1; break;
            case 'G': case 'g': c = 2; break;
            case 'T': case 't': c = 3; break;
            default: c = -1;
        }
        if (c < 0) { valid = 0; fwd = rev = 0; continue; }
        fwd = ((fwd << 2) | (uint64_t)c) & mask;
        rev = (rev >> 2) | ((uint64_t)(3 - c) << (2 * (K - 1)));
        if (++valid >= K) {
            uint64_t h = hash64(min(fwd, rev));
            if ((h & SAMPLE_MASK) == 0) sampled.push_back(make_pair((int64_t)i + 1 - K, h));
        }
    }
    return sampled;
}

int64_t load_no_patch_fasta(const string& path, NoPatch& no_patch) {
    ifstream in(path);
    if (!in) return -1;
    string line, seq;
    int64_t n = 0;
    auto flush = [&]() { if (!seq.empty()) { sample_kmers(seq, no_patch.kmers); seq.clear(); } };
    while (getline(in, line)) {
        if (!line.empty() && line[0] == '>') { flush(); ++n; continue; }
        // bases only: a '\r' (a CRLF file) or a space would otherwise break the k-mers at every line
        for (char c : line) if (isalpha((unsigned char)c)) seq += (char)toupper(c);
    }
    flush();
    return n;
}

void add_no_patch_sequence(const string& seq, NoPatch& no_patch) {
    string up(seq);
    for (char& c : up) c = (char)toupper(c);
    sample_kmers(up, no_patch.kmers);
}

double no_patch_fraction(const string& s, const NoPatch* no_patch) {
    const int K = 31;
    // a sequence (or window) with fewer distinct sampled k-mers than this (~250 bp) is too short to judge
    const size_t MIN_SAMPLED = 16;
    if (!no_patch || no_patch->kmers.empty()) return -1.0;
    vector<pair<int64_t, uint64_t>> sampled = sampled_kmer_positions(s);
    // the share of the distinct k-mers among sampled[a, b) that are no-patch k-mers
    auto fraction = [&](size_t a, size_t b) -> double {
        unordered_set<uint64_t> seen;
        size_t hit = 0;
        for (size_t j = a; j < b; ++j) {
            if (seen.insert(sampled[j].second).second && no_patch->kmers.count(sampled[j].second)) ++hit;
        }
        return seen.size() < MIN_SAMPLED ? -1.0 : (double)hit / (double)seen.size();
    };
    double highest = fraction(0, sampled.size());
    // and in each window (overlapping by half; the last aligned to the end of s), where a unit inside a
    // longer graft is not diluted by the sequence around it
    int64_t W = no_patch->window;
    if (W > 0 && (int64_t)s.size() > W) {
        size_t first = 0, last = 0;
        for (int64_t a = 0; a < (int64_t)s.size(); a += max((int64_t)1, W / 2)) {
            int64_t start = min(a, (int64_t)s.size() - W);
            while (first < sampled.size() && sampled[first].first < start) ++first;
            last = max(last, first);
            while (last < sampled.size() && sampled[last].first + K <= start + W) ++last;
            highest = max(highest, fraction(first, last));
            if (start + W >= (int64_t)s.size()) break;
        }
    }
    return highest;
}

// fraction of the (sampled) canonical 31-mer positions in s that are distinct: near 1 for unique sequence,
// near 0 for a tandem array (an alpha or HSat satellite repeats a few k-mers over and over), so a low value
// means k-mer recovery cannot tell a faithful copy from any other array of the same family.  -1 if s has no
// k-mers to sample
static double kmer_distinct_fraction(const string& s) {
    const int K = 31;
    const uint64_t SAMPLE_MASK = 0xF;   // the same ~1/16 subsample as sample_kmers
    unordered_set<uint64_t> seen;
    int64_t total = 0;
    uint64_t fwd = 0, rev = 0, mask = (1ULL << (2 * K)) - 1;
    int valid = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        int c;
        switch (s[i]) {
            case 'A': case 'a': c = 0; break;
            case 'C': case 'c': c = 1; break;
            case 'G': case 'g': c = 2; break;
            case 'T': case 't': c = 3; break;
            default: c = -1;
        }
        if (c < 0) { valid = 0; fwd = rev = 0; continue; }
        fwd = ((fwd << 2) | (uint64_t)c) & mask;
        rev = (rev >> 2) | ((uint64_t)(3 - c) << (2 * (K - 1)));
        if (++valid >= K) {
            uint64_t h = hash64(min(fwd, rev));
            if ((h & SAMPLE_MASK) == 0) { ++total; seen.insert(h); }
        }
    }
    return total ? (double)seen.size() / (double)total : -1.0;
}

// the lowest kmer_distinct_fraction over windows of `window` bp (overlapping by half) along s.  measured over
// a whole graft, a satellite array is diluted by the unique sequence around it (a 60 kb HSat1B array inside a
// 100 kb swap leaves the graft's fraction near 0.6), so the array shows only in the windows that fall inside
// it.  it does not replace the whole-sequence measure: a repeat whose copies lie further apart than half a
// window (a tandem duplication with a long unit) shows over the whole graft but in no window.  -1 if window
// is 0, s is no longer than one window, or no window has enough k-mers to judge (MIN_SAMPLED, ~1.6 kb of
// sequence at the 1/16 subsample)
static double min_window_distinct_fraction(const string& s, int64_t window) {
    if (window <= 0 || (int64_t)s.size() <= window) return -1.0;
    const int K = 31;
    const int64_t MIN_SAMPLED = 100;
    vector<pair<int64_t, uint64_t>> sampled = sampled_kmer_positions(s);
    double lowest = -1.0;
    int64_t step = max((int64_t)1, window / 2);
    size_t first = 0;
    for (int64_t a = 0; a < (int64_t)s.size(); a += step) {
        // the last window is aligned to the end of s, so its tail is judged on a full window too
        int64_t start = min(a, (int64_t)s.size() - window);
        while (first < sampled.size() && sampled[first].first < start) ++first;
        unordered_set<uint64_t> seen;
        int64_t total = 0;
        for (size_t j = first; j < sampled.size() && sampled[j].first + K <= start + window; ++j) {
            ++total; seen.insert(sampled[j].second);
        }
        if (total >= MIN_SAMPLED) {
            double f = (double)seen.size() / (double)total;
            if (lowest < 0 || f < lowest) lowest = f;
        }
        if (start + window >= (int64_t)s.size()) break;
    }
    return lowest;
}

// percent of the distinct (sampled) k-mers in `removed` that also occur in `added`. a cheap,
// strand-independent proxy for how much of the replaced target sequence the graft recapitulates.
static double kmer_recovery(const string& removed, const string& added) {
    unordered_set<uint64_t> aset, rset;
    sample_kmers(added, aset);
    sample_kmers(removed, rset);
    if (rset.empty()) return 0.0;
    int64_t found = 0;
    for (uint64_t h : rset) if (aset.count(h)) ++found;
    return 100.0 * (double)found / (double)rset.size();
}

// Telomere / contig-end patching.
//
// extend_intervals() only extends the two terminal intervals along their *own* paths to the
// contigs' own ends (and a donor's no further than its first N).  If the target assembly is simply missing a telomere
// there, nothing can add it.  This routine grafts a missing telomere in from a foreign
// cover (selected donor assembly) when one is available:
//
//   For each terminal end of the assembly that lacks a telomere, walk inward from the tip
//   along the terminal target interval looking for the nearest node also visited by a
//   foreign cover; if that cover continues outward (in the assembly's frame) and ends in a
//   real telomere, and the end it would make passes the -T telomere check, hand off to it
//   there.  The target's (divergent, capless) sequence beyond the handoff node is replaced
//   by the foreign cover's run to its telomere.
//
// This is graph-coherent: the handoff is a node genuinely shared by both paths, so the join
// is supported by the graph rather than a blind concatenation.
//
// It is also kept small.  Telomeres are missing from an assembly because reads represent them
// poorly, so read-based QC can rarely vouch for what is grafted at a contig end.  A handoff may
// replace at most max_handoff bp of target (-M), and graft at most max_subtelomere bp of donor
// sequence ahead of the donor's telomere (--max-subtelomere-graft): a telomere patch repairs or
// supplies the telomere array, and does not bring in a subtelomere the target lacks.
vector<tuple<step_handle_t, step_handle_t, bool>> extend_to_telomeres(
        const PathHandleGraph* graph,
        const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals,
        const unordered_map<string, vector<path_handle_t>>& sample_covers,
        const vector<string>& sample_names,
        double telo_threshold,
        int64_t max_handoff,
        int64_t max_subtelomere,
        bool verbose,
        bool report_failures,
        const ErrorRegions* error_regions) {

    static const int64_t OUTER = 20000;            // inspect this many bp at a contig tip for a telomere
    static const int64_t TIP_CHECK = 50000;        // validate_telomeres' (-T) window at a patched contig's tip
    static const int64_t REPORT_MARGIN = 3000000;  // keep searching this far past the cap to report a skipped handoff
    static const int64_t MAX_FOREIGN_EXT = 8000000;// follow a foreign cover at most this far to reach its telomere

    if (intervals.empty()) return intervals;

    // candidate foreign cover paths, in sample-priority order (exclude the target sample)
    vector<path_handle_t> foreign;
    for (size_t i = 1; i < sample_names.size(); ++i) {
        auto it = sample_covers.find(sample_names[i]);
        if (it != sample_covers.end()) {
            for (const path_handle_t& p : it->second) foreign.push_back(p);
        }
    }
    if (foreign.empty()) return intervals;
    unordered_set<path_handle_t> foreign_set(foreign.begin(), foreign.end());

    // sequence of the outermost OUTER bp at one end of a path (from_end -> the path_back tip)
    auto terminal_seq = [&](const path_handle_t& F, bool from_end) -> string {
        string s;
        if (from_end) {
            step_handle_t st = graph->path_back(F);
            while (true) {
                s = graph->get_sequence(graph->get_handle_of_step(st)) + s;
                if ((int64_t)s.size() >= OUTER || st == graph->path_begin(F)) break;
                st = graph->get_previous_step(st);
            }
        } else {
            step_handle_t st = graph->path_begin(F);
            while (true) {
                s += graph->get_sequence(graph->get_handle_of_step(st));
                if ((int64_t)s.size() >= OUTER) break;
                step_handle_t nx = graph->get_next_step(st);
                if (nx == graph->path_end(F)) break;
                st = nx;
            }
        }
        return s;
    };

    // precompute, once per cover, whether each of its two ends carries a telomere
    unordered_map<path_handle_t, bool> telo_begin, telo_end;
    for (const path_handle_t& F : foreign) {
        // terminal_seq(F,false) is in path order with the path_begin tip first; terminal_seq(F,true)
        // has the path_back tip last (is_right_end).
        string tseq_b = terminal_seq(F, false);
        string tseq_e = terminal_seq(F, true);
        telo_begin[F] = seq_has_telomere(tseq_b, 0, (int64_t)tseq_b.size(), true, false, telo_threshold);
        telo_end[F]   = seq_has_telomere(tseq_e, 0, (int64_t)tseq_e.size(), true, true,  telo_threshold);
        if (verbose) cerr << "[panpatch] telomere-patch cover " << graph->get_path_name(F)
                          << " telomere begin=" << telo_begin[F] << " end=" << telo_end[F] << endl;
    }

    // Try to find a telomere patch for one terminal interval.
    // Returns true and sets out_s_i (handoff step on the target path) and out_fi (new foreign
    // interval to splice on) if a patch was found.
    auto find_end_patch = [&](const tuple<step_handle_t, step_handle_t, bool>& term,
                              bool is_front,
                              step_handle_t& out_s_i,
                              tuple<step_handle_t, step_handle_t, bool>& out_fi) -> bool {
        path_handle_t P = graph->get_path_handle_of_step(get<0>(term));
        // an end that is already a donor's (an accepted tip graft) gets no telomere graft on top of it: the
        // -M and subtelomere limits would be measured on the donor, and the report would name it the target
        if (graph->get_sample_name(P) != sample_names[0]) return false;
        bool rev = get<2>(term);
        step_handle_t tip_step   = is_front ? get<0>(term) : get<1>(term);
        step_handle_t inner_bound = is_front ? get<1>(term) : get<0>(term);

        // oriented handle of a step in the assembly's 5'->3' frame
        auto asm_handle = [&](const step_handle_t& s) -> handle_t {
            handle_t h = graph->get_handle_of_step(s);
            return rev ? graph->flip(h) : h;
        };
        // step toward the interval interior
        auto inward = [&](const step_handle_t& s) -> step_handle_t {
            if (!is_front) return rev ? graph->get_next_step(s) : graph->get_previous_step(s);
            else           return rev ? graph->get_previous_step(s) : graph->get_next_step(s);
        };

        // would handing off at target step s, to a graft read outward from the handoff (as added_seq
        // is built below), leave a telomere at this end?  Measured as validate_telomeres (-T) measures
        // the patched contig: the outermost TIP_CHECK bp in assembly 5'->3' order.  Read outward, a
        // front graft is the reverse strand of the assembly, so it is reverse-complemented there.
        auto patched_end_capped = [&](const step_handle_t& s, const string& graft) -> bool {
            // the graft's outermost TIP_CHECK bp, then the target from the (shared) handoff node inward
            int64_t have = min((int64_t)graft.size(), TIP_CHECK);
            string tip_part = graft.substr(graft.size() - have);
            if (is_front) {
                reverse(tip_part.begin(), tip_part.end());
                for (char& c : tip_part) {
                    switch (c) {
                        case 'A': c = 'T'; break; case 'C': c = 'G'; break; case 'G': c = 'C'; break; case 'T': c = 'A'; break;
                        case 'a': c = 't'; break; case 'c': c = 'g'; break; case 'g': c = 'c'; break; case 't': c = 'a'; break;
                    }
                }
            }
            vector<string> inner_part;  // target nodes, handoff node first
            for (step_handle_t t = s; have < TIP_CHECK; t = inward(t)) {
                inner_part.push_back(graph->get_sequence(asm_handle(t)));
                have += inner_part.back().size();
                if (t == inner_bound) break;
            }
            string end_seq;
            if (is_front) {
                end_seq = tip_part;
                for (const string& node_seq : inner_part) end_seq += node_seq;
            } else {
                for (auto it = inner_part.rbegin(); it != inner_part.rend(); ++it) end_seq += *it;
                end_seq += tip_part;
            }
            int64_t from = is_front ? 0 : max((int64_t)0, (int64_t)end_seq.size() - TIP_CHECK);
            int64_t to = is_front ? min((int64_t)end_seq.size(), TIP_CHECK) : (int64_t)end_seq.size();
            return seq_has_telomere(end_seq, from, to, true, !is_front, telo_threshold);
        };

        // 1) if the tip already carries a telomere, nothing to do.
        // Build the outermost OUTER bp of this end in true assembly 5'->3' order so the boundary
        // detection in seq_has_telomere agrees with validate_telomeres (front tip -> sequence start,
        // is_right_end=false; back tip -> sequence end, is_right_end=true).  (Walking from the tip
        // inward visits nodes in reverse order for the back end, so prepend there.)
        bool buried = false;  // capless tip, but telomeric repeats are present nearby
        {
            string tip_seq;
            step_handle_t s = tip_step;
            while ((int64_t)tip_seq.size() < OUTER) {
                string node_seq = graph->get_sequence(asm_handle(s));
                if (is_front) tip_seq += node_seq;          // front tip at index 0
                else          tip_seq = node_seq + tip_seq;  // back tip at the end
                if (s == inner_bound) break;
                s = inward(s);
            }
            bool has = seq_has_telomere(tip_seq, 0, (int64_t)tip_seq.size(), true, !is_front, telo_threshold);
            if (verbose) cerr << "[panpatch] telomere-patch " << (is_front ? "front" : "back")
                              << " tip of " << graph->get_path_name(P) << ": telomere=" << has << endl;
            if (has) return false;
            buried = has_telomeric_window(tip_seq, telo_threshold);
        }

        // 2) walk inward looking for a shared-node handoff to a foreign cover with a telomere
        int64_t walked = 0;
        step_handle_t s_i = tip_step;

        // covers (by the direction they are read in) whose run out to their telomere holds an N (an
        // assembly gap) or crosses an
        // --error-bed interval on the cover.  a graft from them would put those N's, or that flagged
        // sequence, in the patched assembly.  further in, the same cover's run out is only longer, so it
        // is skipped at every later handoff too.  if no other cover completes the end, the first such
        // cover is reported, so a missing telomere is never silently down to N's or a flag
        unordered_set<path_handle_t> n_tainted_fwd, n_tainted_rev;
        PatchRecord n_reject;
        bool have_n_reject = false;
        // covers (by direction) whose graft would carry more than max_subtelomere bp of sequence ahead of
        // their telomere.  further in, the same cover's graft only grows on that (subtelomeric) side, so it
        // too is skipped at every later handoff, and the first one is reported if nothing completes the end
        unordered_set<path_handle_t> sub_tainted_fwd, sub_tainted_rev;
        PatchRecord sub_reject;
        bool have_sub_reject = false;
        auto report_skipped = [&]() {
            if (have_n_reject) g_patch_records.push_back(n_reject);
            if (have_sub_reject) g_patch_records.push_back(sub_reject);
        };
        // a rejected handoff to F that would replace the target's outermost `walked` bp
        auto rejected_handoff = [&](const path_handle_t& F, int64_t walked, const string& reason) {
            PatchRecord pr;
            pr.type = "telomere"; pr.target = graph->get_path_name(P); pr.target_bp = path_bp(graph, P);
            pr.donor = graph->get_path_name(F); pr.donor_bp = path_bp(graph, F);
            pr.replaced_bp = walked; pr.accepted = false;
            if (is_front != rev) { pr.target_start = 0; pr.target_end = walked; }
            else { pr.target_start = pr.target_bp - walked; pr.target_end = pr.target_bp; }
            pr.reason = reason;
            return pr;
        };
        while (true) {
            handle_t a_i = asm_handle(s_i);
            handle_t under = graph->get_handle_of_step(s_i);
            handle_t out_handle = is_front ? graph->flip(a_i) : a_i;

            // collect the foreign covers present on this node in one pass.  a cover can visit the node
            // more than once, and the graph lists the steps in no particular order, so try each visit
            // in path order (the first usable one wins)
            unordered_map<path_handle_t, vector<step_handle_t>> steps_here;
            graph->for_each_step_on_handle(under, [&](step_handle_t fs) {
                path_handle_t fp = graph->get_path_handle_of_step(fs);
                if (foreign_set.count(fp)) steps_here[fp].push_back(fs);
            });

            // the handoffs to try here: the covers in priority order, each one's visits in path order
            vector<pair<path_handle_t, step_handle_t>> handoffs;
            for (const path_handle_t& F : foreign) {
                auto sit = steps_here.find(F);
                if (sit == steps_here.end()) continue;  // cover F does not visit this node
                vector<step_handle_t>& F_steps = sit->second;
                if (F_steps.size() > 1) {
                    sort(F_steps.begin(), F_steps.end(), [&](const step_handle_t& a, const step_handle_t& b) {
                        return step_offset(graph, a) < step_offset(graph, b);
                    });
                }
                for (const step_handle_t& s_F : F_steps) {
                    handoffs.push_back(make_pair(F, s_F));
                }
            }

            for (const pair<path_handle_t, step_handle_t>& handoff : handoffs) {
                const path_handle_t& F = handoff.first;
                const step_handle_t& s_F = handoff.second;

                // pick the direction along F that continues outward in the assembly frame
                handle_t gF = graph->get_handle_of_step(s_F);
                bool f_forward;
                if (gF == out_handle) f_forward = true;
                else if (gF == graph->flip(out_handle)) f_forward = false;
                else continue;

                // the outward end of F in this direction must be telomeric (precomputed)
                if (!(f_forward ? telo_end[F] : telo_begin[F])) continue;

                // its run out to that end must be free of N's (see n_tainted_* above)
                unordered_set<path_handle_t>& n_tainted = f_forward ? n_tainted_fwd : n_tainted_rev;
                if (n_tainted.count(F)) continue;
                unordered_set<path_handle_t>& sub_tainted = f_forward ? sub_tainted_fwd : sub_tainted_rev;
                if (sub_tainted.count(F)) continue;

                // `walked` is the length of the target's capless tail this handoff would replace.
                bool within_cap = (walked <= max_handoff);

                // follow F outward to its (telomeric) end, bounded; this confirms reachability and
                // gives us the terminal step. accumulate the grafted sequence only when we will
                // actually patch (within cap), for the recovery stat.
                int64_t ext = 0;
                string added_seq;
                step_handle_t fs = f_forward ? graph->get_next_step(s_F) : graph->get_previous_step(s_F);
                step_handle_t f_last = s_F;
                // the handoff node comes out on the donor's side of the splice, so it counts too: where
                // the target's and the donor's N runs were aligned, it is a shared N node
                bool reached = false, any = false;
                bool has_n = graph->get_sequence(graph->get_handle_of_step(s_F)).find_first_of("Nn") != string::npos;
                bool flagged = !has_n && error_regions && step_in_error(graph, s_F, *error_regions);
                while (!has_n && !flagged) {
                    bool at_end = f_forward ? (fs == graph->path_end(F)) : (fs == graph->path_front_end(F));
                    if (at_end) { reached = true; break; }
                    handle_t fh = graph->get_handle_of_step(fs);
                    string fseq = graph->get_sequence(f_forward ? fh : graph->flip(fh));
                    if (fseq.find_first_of("Nn") != string::npos) { has_n = true; break; }
                    if (error_regions && step_in_error(graph, fs, *error_regions)) { flagged = true; break; }
                    if (within_cap) added_seq += fseq;
                    ext += fseq.size();
                    f_last = fs; any = true;
                    if (ext > MAX_FOREIGN_EXT) break;
                    fs = f_forward ? graph->get_next_step(fs) : graph->get_previous_step(fs);
                }
                if (has_n || flagged) {
                    n_tainted.insert(F);
                    if (!have_n_reject) {
                        n_reject = rejected_handoff(F, walked, has_n ?
                                                    "donor's run to its telomere contains N (an assembly gap)" :
                                                    "donor's run to its telomere crosses one of its --error-bed intervals");
                        have_n_reject = true;
                    }
                    if (verbose) cerr << "[panpatch] telomere-patch " << (is_front ? "front" : "back")
                                      << " handoff to " << graph->get_path_name(F) << " at " << walked
                                      << "bp in would graft " << (has_n ? "N's" : "flagged sequence") << "; skipping that cover" << endl;
                    continue;
                }
                if (!reached || !any) continue;  // not a usable handoff; keep looking

                // the patched end must carry a telomere.  where the target's own (short or degraded) array
                // is aligned to the donor's, the nearest shared node is inside the arrays and the graft is
                // a few bp that leave the end as capless as it was (-T would then revert the contig), so
                // keep looking: the next cover here, then further in, where the graft takes the donor's array
                if (within_cap && !patched_end_capped(s_i, added_seq)) {
                    if (verbose) cerr << "[panpatch] telomere-patch " << (is_front ? "front" : "back")
                                      << " handoff to " << graph->get_path_name(F) << " at " << walked
                                      << "bp in leaves the end capless; looking further in" << endl;
                    continue;
                }

                if (!within_cap) {
                    // nearest usable handoff is too far in: replacing this much target is risky, so skip
                    g_patch_records.push_back(rejected_handoff(F, walked, "handoff would replace " + std::to_string(walked) +
                                                               "bp, over --max-telomere-patch " + std::to_string(max_handoff)));
                    report_skipped();
                    return false;
                }

                // the graft must be the donor's telomere array, not a subtelomere the target lacks: added_seq
                // is read outward, so the donor's telomere is at its end, and what precedes the telomeric run
                // found there is subtelomere (all of it, for a graft too short to hold a run)
                if (max_subtelomere > 0) {
                    int64_t run_start = find_telomere_boundary(added_seq, 0, (int64_t)added_seq.size(), false);
                    int64_t sub_bp = run_start < 0 ? (int64_t)added_seq.size() : run_start;
                    if (sub_bp > max_subtelomere) {
                        sub_tainted.insert(F);
                        if (!have_sub_reject) {
                            sub_reject = rejected_handoff(F, walked, "donor's graft carries " + std::to_string(sub_bp) +
                                                          "bp ahead of its telomere, over --max-subtelomere-graft " +
                                                          std::to_string(max_subtelomere));
                            have_sub_reject = true;
                        }
                        if (verbose) cerr << "[panpatch] telomere-patch " << (is_front ? "front" : "back")
                                          << " handoff to " << graph->get_path_name(F) << " at " << walked
                                          << "bp in would graft " << sub_bp << "bp of subtelomere; skipping that cover" << endl;
                        continue;
                    }
                }

                // replaced target tail (tip .. just inside the handoff node) for the recovery stat
                string removed_seq;
                for (step_handle_t s = tip_step; s != s_i; s = inward(s)) {
                    removed_seq += graph->get_sequence(asm_handle(s));
                }
                double recovery = kmer_recovery(removed_seq, added_seq);
                {
                    PatchRecord pr;
                    pr.type = "telomere"; pr.target = graph->get_path_name(P); pr.target_bp = path_bp(graph, P);
                    pr.donor = graph->get_path_name(F); pr.donor_bp = path_bp(graph, F);
                    pr.replaced_bp = walked; pr.kmer = recovery;
                    // localize the replaced tail in target-contig fwd coords: the tip sits at the contig
                    // front when (is_front != rev) -- a reverse interval swaps the front/back convention --
                    // so it replaced [0, walked); otherwise it is the back end, [L - walked, L)
                    if (is_front != rev) { pr.target_start = 0; pr.target_end = walked; }
                    else { pr.target_start = pr.target_bp - walked; pr.target_end = pr.target_bp; }
                    g_patch_records.push_back(pr);
                }

                out_s_i = s_i;
                out_fi  = is_front ? make_tuple(f_last, s_F, f_forward)
                                   : make_tuple(s_F, f_last, !f_forward);
                return true;
            }

            if (s_i == inner_bound) break;
            walked += graph->get_length(under);
            if (walked > max_handoff + REPORT_MARGIN) break;
            s_i = inward(s_i);
        }

        // capless end with no usable donor handoff in range.  report_failures (set under -T, where a
        // missing telomere is why a contig reverts) records why -- which ends were left as-is, and
        // whether it's a simple gap (no telomere anywhere, no donor) or a telomere present but
        // buried/degraded at the tip.  when telomere patching is just the default (not required), this
        // is a non-event, so stay quiet.  (a cover skipped for N's or for its subtelomere is always
        // reported: it is a donor that would have completed the end.)
        report_skipped();
        if (report_failures) {
            PatchRecord pr;
            pr.type = "telomere"; pr.target = graph->get_path_name(P); pr.target_bp = path_bp(graph, P);
            pr.accepted = false;
            pr.reason = buried ? "telomere present but buried/degraded; no clean donor telomere here"
                               : "no telomere at this end and no donor reaches one";
            g_patch_records.push_back(pr);
        }
        return false;
    };

    step_handle_t s_i_front, s_i_back;
    tuple<step_handle_t, step_handle_t, bool> fi_front, fi_back;
    size_t front_rows = g_patch_records.size();
    bool fp = find_end_patch(intervals.front(), true,  s_i_front, fi_front);
    bool bp = find_end_patch(intervals.back(),  false, s_i_back,  fi_back);
    if (!fp && !bp) return intervals;
    // on a contig shorter than the two handoffs' reach, both ends' handoffs can be on the one interval and
    // cross each other, which would leave an interval running backwards: take neither
    if (fp && bp && intervals.size() == 1) {
        int64_t of = step_offset(graph, s_i_front), ob = step_offset(graph, s_i_back);
        if (get<2>(intervals[0]) ? of < ob : of > ob) {
            for (size_t r = front_rows; r < g_patch_records.size(); ++r) {
                if (g_patch_records[r].type == "telomere" && g_patch_records[r].accepted) {
                    g_patch_records[r].accepted = false;
                    g_patch_records[r].reason = "the contig is too short for both ends' telomere handoffs (they cross)";
                }
            }
            return intervals;
        }
    }

    // splice in the foreign extension(s); the foreign interval owns its (telomere) tip and the
    // shared handoff node is kept exactly once at the boundary.
    vector<tuple<step_handle_t, step_handle_t, bool>> out;
    size_t n = intervals.size();
    if (n == 1) {
        auto I = intervals[0];
        step_handle_t a = get<0>(I), b = get<1>(I);
        bool r = get<2>(I);
        if (fp) a = s_i_front;
        if (bp) b = s_i_back;
        if (fp) out.push_back(fi_front);
        out.push_back(make_tuple(a, b, r));
        if (bp) out.push_back(fi_back);
    } else {
        if (fp) out.push_back(fi_front);
        {
            auto I = intervals[0];
            if (fp) I = make_tuple(s_i_front, get<1>(I), get<2>(I));
            out.push_back(I);
        }
        for (size_t i = 1; i + 1 < n; ++i) out.push_back(intervals[i]);
        {
            auto I = intervals[n - 1];
            if (bp) I = make_tuple(get<0>(I), s_i_back, get<2>(I));
            out.push_back(I);
        }
        if (bp) out.push_back(fi_back);
    }

    // Truncating a target interval to the handoff node can leave it zero-length when the handoff
    // lands on its inner boundary. Such a non-last interval emits nothing (the boundary node is
    // re-emitted by its neighbour, so the sequence is unaffected), but it would still print a
    // spurious zero-length BED line - drop it.
    vector<tuple<step_handle_t, step_handle_t, bool>> trimmed;
    for (size_t i = 0; i < out.size(); ++i) {
        bool empty = (get<0>(out[i]) == get<1>(out[i]));
        bool is_last = (i + 1 == out.size());
        if (empty && !is_last) continue;
        trimmed.push_back(out[i]);
    }
    return trimmed;
}

void print_intervals(const PathHandleGraph* graph,
                     const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals,
                     ostream& out) {

    // manually index all paths in the interval cover
    // todo: use some kind of position overlay!
    unordered_map<path_handle_t, unordered_map<step_handle_t, int64_t>> path_index;
    for (const auto& interval : intervals) {
        path_index[graph->get_path_handle_of_step(get<0>(interval))] = {};
    }
    for (auto& path_map : path_index) {
        int64_t pos = 0;
        graph->for_each_step_in_path(path_map.first, [&](step_handle_t step) {
            path_map.second[step] = pos;
            pos += graph->get_length(graph->get_handle_of_step(step));
        });
    }

    // hack:  all intervals are open-ended except the last
    // todo:  fix upstream!
    for (int64_t i = 0; i < intervals.size(); ++i) {
        const auto& interval = intervals[i];
        path_handle_t path = graph->get_path_handle_of_step(get<0>(interval));
        int64_t pos_1, pos_2;
        if (get<2>(interval) == false) {
            pos_1 = path_index[path][get<0>(interval)];
            pos_2 = path_index[path][get<1>(interval)];
            if (i == intervals.size() - 1) {
                pos_2 += graph->get_length(graph->get_handle_of_step(get<1>(interval)));
            }
        } else {
            pos_1 = path_index[path][get<1>(interval)] + graph->get_length(graph->get_handle_of_step(get<1>(interval)));
            pos_2 = path_index[path][get<0>(interval)] + graph->get_length(graph->get_handle_of_step(get<0>(interval)));
            if (i == intervals.size() -1) {
                pos_1 -= graph->get_length(graph->get_handle_of_step(get<1>(interval)));
            }
        }
        out << graph->get_path_name(path) << "\t" << pos_1 << "\t" << pos_2
             << "\t" << (get<2>(interval) ? '-' : '+') << endl;

    }
}

string intervals_to_sequence(const PathHandleGraph* graph,
                             const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals) {
    
    string seq;
    // hack:  all intervals are open-ended except the last
    // todo:  fix upstream!
    for (int64_t i = 0; i < intervals.size(); ++i) {
        const auto& interval = intervals[i];
#ifdef debug
        cerr << "Interval " << graph->get_id(graph->get_handle_of_step(get<0>(interval))) << ":"
             << graph->get_is_reverse(graph->get_handle_of_step(get<0>(interval))) << " - "
             << graph->get_id(graph->get_handle_of_step(get<1>(interval))) << ":"
             << graph->get_is_reverse(graph->get_handle_of_step(get<1>(interval)))
             << " rev=" <<get<2>(interval) << endl;
#endif
        if (get<2>(interval) == false) {            
            step_handle_t last_step = get<1>(interval);
            if (i == intervals.size() - 1) {
                last_step = graph->get_next_step(last_step);
            }
            for (step_handle_t step = get<0>(interval); step != last_step; step = graph->get_next_step(step)) {
                seq += graph->get_sequence(graph->get_handle_of_step(step));
            }
        } else {
            step_handle_t last_step = get<1>(interval);
            if (i == intervals.size() - 1) {
                last_step = graph->get_previous_step(last_step);
            }
            for (step_handle_t step = get<0>(interval); step != last_step; step = graph->get_previous_step(step)) {
                seq += graph->get_sequence(graph->flip(graph->get_handle_of_step(step)));
            }            
        }
    }
    return seq;
}

vector<tuple<step_handle_t, step_handle_t, bool>> greedy_patch(const PathHandleGraph* graph,
                                                               const path_handle_t& ref_path,
                                                               const vector<path_handle_t>& tgt_paths,
                                                               const vector<string>& sample_names,
                                                               const unordered_map<string, vector<path_handle_t>>& sample_covers,
                                                               const BedRegions& bed_regions,
                                                               bool patch_ends,
                                                               double telo_threshold,
                                                               int64_t max_telomere_patch,
                                                               int64_t max_subtelomere_graft,
                                                               bool verbose,
                                                               bool report_telomere_failures,
                                                               unordered_set<path_handle_t>* threaded_tgt_paths,
                                                               const ErrorRegions* error_regions,
                                                               bool patch_gap,
                                                               int64_t max_graft_pieces,
                                                               const NoPatch* no_patch) {

    size_t records_at_entry = g_patch_records.size();   // this call's report rows start here

#ifdef debug
    cerr << "greedy patch\n"
         << " ref_path = " << graph->get_path_name(ref_path) << endl
         << " tgt_paths =";
    for (const path_handle_t& tgt_path : tgt_paths) {
        cerr << " " << graph->get_path_name(tgt_path);
    }
    cerr << endl << " other_paths =";
    for (int64_t i = 1; i < sample_names.size(); ++i) {
        const string& sample_name = sample_names[i];
        if (sample_covers.count(sample_name)) {
            for (const path_handle_t& other_path : sample_covers.at(sample_name)) {
                cerr << " " << graph->get_path_name(other_path);
            }
        }
    }
    cerr << endl;
#endif

    // find the anchors along the reference path using the relevant paths.  the target's own paths are
    // always included: the cover chosen for the target sample can be a different haplotype (on a tie
    // in coverage, or when the target's contigs overlap each other), and without them there would be
    // no anchors at all
    unordered_set<path_handle_t> relevant_paths = {ref_path};
    relevant_paths.insert(tgt_paths.begin(), tgt_paths.end());
    for (const auto& sample_cover : sample_covers) {
        for (const path_handle_t& path : sample_cover.second) {
            relevant_paths.insert(path);
        }
    }
    unordered_map<int64_t, int64_t> ref_anchors = find_anchors(graph, ref_path, tgt_paths, relevant_paths);

    // convert BED exclusion regions to reference coordinates
    ExcludedRefRegions excluded_regions;
    if (!bed_regions.empty()) {
        excluded_regions = bed_to_ref_regions(graph, tgt_paths, ref_anchors, bed_regions);
#ifdef debug
        cerr << "excluded " << excluded_regions.size() << " ref regions from BED file" << endl;
        for (const auto& r : excluded_regions) {
            cerr << "  [" << r.first << ", " << r.second << ")" << endl;
        }
#endif
    }

#ifdef debug

    cerr << "number of anchors found " << ref_anchors.size() << endl;
#endif


    // todo: stop copying paths lists into so many different structures!    
    vector<path_handle_t> other_paths;
    for (int64_t i = 1; i < sample_names.size(); ++i) {
        const string& sample_name = sample_names[i];
        if (sample_covers.count(sample_name)) {
            for (const path_handle_t& other_path : sample_covers.at(sample_name)) {
                other_paths.push_back(other_path);
            }
        }
    }
                
    vector<tuple<step_handle_t, step_handle_t, bool>> patched_intervals = thread_intervals(graph,
                                                                                           ref_path,
                                                                                           ref_anchors,
                                                                                           tgt_paths,
                                                                                           other_paths,
                                                                                           excluded_regions,
                                                                                           error_regions);

    // the target contigs the threading went through: a telomere graft or trim below can take one out of
    // the final intervals again, but it was still patched over, unlike a contig never threaded at all
    if (threaded_tgt_paths) {
        unordered_set<path_handle_t> tgt_set(tgt_paths.begin(), tgt_paths.end());
        for (const auto& iv : patched_intervals) {
            path_handle_t p = graph->get_path_handle_of_step(get<0>(iv));
            if (tgt_set.count(p)) threaded_tgt_paths->insert(p);
        }
    }

    if (patched_intervals.empty()) {
        cerr << "[panpatch] warning: unable to patch assembly on " << graph->get_locus_name(ref_path)
             << " for " << graph->get_sample_name(tgt_paths.front()) << "#"
         << graph->get_haplotype(tgt_paths.front()) << ":" << endl;
        // say so in the report, whatever the revert then writes (the input, or with -e the default sample)
        PatchRecord pr;
        pr.type = "gap-fill"; pr.target = graph->get_path_name(tgt_paths.front()); pr.target_bp = path_bp(graph, tgt_paths.front());
        pr.accepted = false;
        pr.reason = "the target could not be threaded against the reference (it shares too few anchors with it)";
        g_patch_records.push_back(pr);
        return patched_intervals;
    }

    check_intervals(graph, patched_intervals);

#ifdef debug
    cerr << "number of patched intervals found " << patched_intervals.size() << endl;
#endif
    vector<tuple<step_handle_t, step_handle_t, bool>> smoothed_intervals = smooth_intervals(graph, patched_intervals);

    check_intervals(graph, smoothed_intervals);

#ifdef debug
    cerr << "number of smoothed intervals found " << smoothed_intervals.size() << endl;
#endif

    vector<tuple<step_handle_t, step_handle_t, bool>> extended_intervals = extend_intervals(graph, smoothed_intervals, tgt_paths,
                                                                                            error_regions);

    check_intervals(graph, extended_intervals);

    // tip grafts: a donor at either end of the threading (the target could not get past an N gap or a
    // flagged region near its tip) carries the patch out past the target's outermost piece, to the donor's
    // own end or its first N or flagged step.  nothing anchors the far end of such a graft, so it is held to
    // the limits of a telomere graft: it may replace at most max_telomere_patch (-M) bp of the target's end,
    // may add at most max_subtelomere_graft bp more than it replaces, and may not take off a telomere the
    // target already has.  (unchecked, a flag on a contig's last few kb appended the donor's whole tail.)  one
    // that fails is dropped -- the target keeps its own end -- and gets a rejected row; one that passes gets
    // an accepted gap-fill row
    {
        unordered_set<path_handle_t> tset(tgt_paths.begin(), tgt_paths.end());
        auto judge_tip = [&](bool at_front) {
            int n = (int)extended_intervals.size();
            int tk = -1;
            for (int k = 0; k < n; ++k) {
                int kk = at_front ? k : n - 1 - k;
                if (tset.count(graph->get_path_handle_of_step(get<0>(extended_intervals[kk])))) { tk = kk; break; }
            }
            if (tk < 0 || tk == (at_front ? 0 : n - 1)) return;   // no target piece, or the target is at this end
            path_handle_t P = graph->get_path_handle_of_step(get<0>(extended_intervals[tk]));
            path_handle_t D = graph->get_path_handle_of_step(get<0>(extended_intervals[at_front ? 0 : n - 1]));
            bool rev = get<2>(extended_intervals[tk]);
            // the target piece's forward range
            step_handle_t lo_step = rev ? get<1>(extended_intervals[tk]) : get<0>(extended_intervals[tk]);
            step_handle_t hi_step = rev ? get<0>(extended_intervals[tk]) : get<1>(extended_intervals[tk]);
            int64_t lo = step_offset(graph, lo_step);
            int64_t hi = step_offset(graph, hi_step) + graph->get_length(graph->get_handle_of_step(hi_step));
            int64_t L = path_bp(graph, P);
            // the output's front is the target's start for a forward piece, its end for a reverse one
            bool target_start_side = (at_front != rev);
            int64_t t_start = target_start_side ? 0 : hi, t_end = target_start_side ? lo : L;
            int64_t replaced = t_end - t_start;
            // the donor sequence the graft adds, less the node it shares with the target piece
            vector<tuple<step_handle_t, step_handle_t, bool>> tip(
                at_front ? extended_intervals.begin() : extended_intervals.begin() + tk + 1,
                at_front ? extended_intervals.begin() + tk : extended_intervals.end());
            step_handle_t junction = at_front ? get<1>(tip.back()) : get<0>(tip.front());
            string tip_seq = intervals_to_sequence(graph, tip);
            int64_t added = (int64_t)tip_seq.size() - (int64_t)graph->get_length(graph->get_handle_of_step(junction));
            // the target's own end the graft would replace
            string tail;
            if (replaced > 0) {
                int64_t off = 0;
                graph->for_each_step_in_path(P, [&](step_handle_t st) {
                    int64_t l = graph->get_length(graph->get_handle_of_step(st));
                    int64_t a = max(off, t_start), b = min(off + l, t_end);
                    if (a < b) tail += graph->get_sequence(graph->get_handle_of_step(st)).substr(a - off, b - a);
                    off += l;
                    return off < t_end;
                });
            }
            double np_added = no_patch_fraction(tip_seq, no_patch), np_replaced = no_patch_fraction(tail, no_patch);
            double np = max(np_added, np_replaced);
            string why;
            if (np >= 0 && np >= no_patch->fraction) {
                ostringstream ws; ws << fixed << setprecision(1) << (np_replaced > np_added ? "replaces" : "adds")
                                     << " sequence that shares " << 100.0 * np
                                     << "% of its k-mers with the no-patch sequences (rDNA, --no-patch-fasta)";
                why = ws.str();
            } else if (!patch_gap) {
                // a tip graft is gap filling: --patch-types without gap excludes it, as excise does interior ones
                why = "excluded by --patch-types (gap filling is off)";
            } else if (max_graft_pieces > 0 && (int64_t)tip.size() > max_graft_pieces) {
                why = "is spliced from " + to_string(tip.size()) + " pieces of donor sequence, over --max-graft-pieces " +
                    to_string(max_graft_pieces) + " (no donor has it in one piece)";
            } else if (replaced > max_telomere_patch) {   // (as for a telomere graft, 0 lets it replace nothing)
                why = "replaces " + to_string(replaced) + "bp of the target's end, over --max-telomere-patch " +
                    to_string(max_telomere_patch);
            } else if (max_subtelomere_graft > 0 && added - replaced > max_subtelomere_graft) {
                why = "adds " + to_string(added) + "bp of donor sequence for " + to_string(replaced) +
                    "bp of the target's end, more than --max-subtelomere-graft " + to_string(max_subtelomere_graft) +
                    " beyond it";
            } else if (replaced > 0) {
                // the target's own end it would replace (<= -M bp): does it carry a telomere?
                if (seq_has_telomere(tail, 0, (int64_t)tail.size(), true, !target_start_side, telo_threshold)) {
                    why = "would take off the telomere at the target's end";
                }
            }
            PatchRecord pr;
            pr.type = "gap-fill"; pr.target = graph->get_path_name(P); pr.target_bp = L;
            pr.donor = graph->get_path_name(D); pr.donor_bp = path_bp(graph, D);
            pr.target_start = t_start; pr.target_end = t_end; pr.replaced_bp = replaced;
            pr.accepted = why.empty();
            if (!why.empty()) pr.reason = "tip graft " + why + "; the target keeps its own end";
            g_patch_records.push_back(pr);
            if (!why.empty()) {
                // drop the donor's tip and run the target piece out to its own end, as extend_intervals
                // does for a target that reaches it
                if (at_front) {
                    get<0>(extended_intervals[tk]) = rev ? graph->path_back(P) : graph->path_begin(P);
                    extended_intervals.erase(extended_intervals.begin(), extended_intervals.begin() + tk);
                } else {
                    get<1>(extended_intervals[tk]) = rev ? graph->path_begin(P) : graph->path_back(P);
                    extended_intervals.erase(extended_intervals.begin() + tk + 1, extended_intervals.end());
                }
            }
        };
        judge_tip(true);
        judge_tip(false);
    }

    if (patch_ends) {
        extended_intervals = extend_to_telomeres(graph, extended_intervals, sample_covers,
                                                 sample_names, telo_threshold, max_telomere_patch,
                                                 max_subtelomere_graft, verbose,
                                                 report_telomere_failures, error_regions);
        check_intervals(graph, extended_intervals);
    }

    // backstop: every interval must run forward along its path in its own direction, or no sequence can be
    // read from it (intervals_to_sequence would walk off the path's end).  should one not, give up on this
    // patch -- the contig is kept as it is -- rather than crash or write a patch spliced out of order
    for (const auto& iv : extended_intervals) {
        int64_t a = step_offset(graph, get<0>(iv)), b = step_offset(graph, get<1>(iv));
        if (get<2>(iv) ? a < b : a > b) {
            path_handle_t p = graph->get_path_handle_of_step(get<0>(iv));
            cerr << "[panpatch] warning: patch of " << graph->get_path_name(tgt_paths.front()) << " has an interval of "
                 << graph->get_path_name(p) << " running " << a << " -> " << b << (get<2>(iv) ? " (reverse)" : "")
                 << "; keeping the contig as it is" << endl;
            g_patch_records.resize(records_at_entry);
            PatchRecord pr;
            pr.type = "gap-fill"; pr.target = graph->get_path_name(tgt_paths.front()); pr.target_bp = path_bp(graph, tgt_paths.front());
            pr.donor = graph->get_path_name(p); pr.donor_bp = path_bp(graph, p);
            pr.accepted = false;
            pr.reason = "patch splices this path out of order (internal); patch abandoned";
            g_patch_records.push_back(pr);
            return {};
        }
    }

    return extended_intervals;
}

// Sanity guard for repeat-region misjoins.
//
// Detects when a contig is used in two or more disjoint pieces and another contig OF THE SAME
// SAMPLE is spliced into the interior between them.  This is the signature of the acrocentric /
// pericentromeric misjoins: the threading bounces through ambiguous satellite/segdup anchors and
// splices the target assembly's own spare fragments into the middle of a contig that already spans
// the region, producing scrambled / collapsed output.
//
// Legitimate operations are unaffected:
//  - a genuine gap-fill bridges with a FOREIGN donor (different sample), so the interior material
//    is not same-sample and is allowed;
//  - end-to-end scaffolds and telomere patches use each contig contiguously (a single block), so
//    there is no interior to splice into.
// Restricted to the target sample's own contigs: the misjoin we care about is the target assembly's
// spare fragments spliced into the target's main contig. (This also avoids a donor-vs-donor case --
// e.g. one donor grafted at both telomeres with a same-donor gap-fill between -- reverting a valid
// telomere-completed assembly.)
static bool splices_same_sample_interior(const PathHandleGraph* graph,
                                         const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals,
                                         const string& target_sample,
                                         string& detail) {
    // first and last position at which each contig appears in the output (list) order; an interval
    // between those two positions that is on a different contig is "interior" material
    unordered_map<path_handle_t, pair<int, int>> span;
    vector<path_handle_t> idx_path(intervals.size());
    for (int i = 0; i < (int)intervals.size(); ++i) {
        path_handle_t p = graph->get_path_handle_of_step(get<0>(intervals[i]));
        idx_path[i] = p;
        auto it = span.find(p);
        if (it == span.end()) span[p] = make_pair(i, i);
        else it->second.second = i;
    }
    // contigs in output order, so the one reported is the first, not whichever the hash table gives
    for (const path_handle_t& C : paths_in_order(graph, intervals)) {
        int lo = span.at(C).first, hi = span.at(C).second;
        if (hi <= lo) continue;  // contig used as a single contiguous block: no interior to splice
        string c_sample = graph->get_sample_name(C);
        if (c_sample != target_sample) continue;  // only judge the target's own contigs
        for (int j = lo + 1; j < hi; ++j) {
            if (idx_path[j] == C) continue;            // another piece of the same contig
            if (graph->get_sample_name(idx_path[j]) == c_sample) {
                detail = "contig " + graph->get_path_name(C)
                       + " was used non-contiguously with same-sample fragment "
                       + graph->get_path_name(idx_path[j]) + " spliced into its interior";
                return true;
            }
        }
    }
    return false;
}

// Second guard (complements splices_same_sample_interior, which only catches same-sample
// interior splices).  When a *foreign* donor is spliced into the interior of a target contig (a
// gap-fill / replacement), measure how much of the replaced target sequence the graft actually
// recapitulates with the same canonical-31-mer recovery used for telomere grafts.  A graft that
// shares almost none of the replaced sequence's k-mers came from a different locus/paralog (a
// repeat-region misjoin), not a faithful fill -- revert it.  Only interior replacements are judged;
// terminal telomere grafts legitimately have low recovery (divergent subtelomeres) and are not
// considered here.  Restricted to the target sample's own contigs (avoids donor-vs-donor cases).
static bool interior_graft_low_recovery(const PathHandleGraph* graph,
                                        const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals,
                                        const string& target_sample,
                                        double min_recovery, int64_t min_replaced,
                                        const unordered_map<path_handle_t, int64_t>& excised_nonN,
                                        string& detail) {
    // step -> forward position index per contig used
    unordered_map<path_handle_t, unordered_map<step_handle_t, int64_t>> path_index;
    for (const auto& iv : intervals) path_index[graph->get_path_handle_of_step(get<0>(iv))] = {};
    for (auto& pm : path_index) {
        int64_t pos = 0;
        graph->for_each_step_in_path(pm.first, [&](step_handle_t s) {
            pm.second[s] = pos; pos += graph->get_length(graph->get_handle_of_step(s));
        });
    }
    // used forward ranges + list-index span per contig
    unordered_map<path_handle_t, vector<pair<int64_t, int64_t>>> used;
    unordered_map<path_handle_t, pair<int, int>> span;
    for (int i = 0; i < (int)intervals.size(); ++i) {
        const auto& iv = intervals[i];
        path_handle_t p = graph->get_path_handle_of_step(get<0>(iv));
        int64_t a, b;
        if (!get<2>(iv)) { a = path_index[p][get<0>(iv)]; b = path_index[p][get<1>(iv)] + graph->get_length(graph->get_handle_of_step(get<1>(iv))); }
        else             { a = path_index[p][get<1>(iv)]; b = path_index[p][get<0>(iv)] + graph->get_length(graph->get_handle_of_step(get<0>(iv))); }
        used[p].push_back(make_pair(a, b));
        auto it = span.find(p); if (it == span.end()) span[p] = make_pair(i, i); else it->second.second = i;
    }
    // contigs in output order, so the one reported is the first, not whichever the hash table gives
    for (const path_handle_t& C : paths_in_order(graph, intervals)) {
        if (graph->get_sample_name(C) != target_sample) continue;     // only the target's own contigs
        auto& r = used.at(C);
        if (r.size() < 2) continue;
        sort(r.begin(), r.end());
        vector<pair<int64_t, int64_t>> gaps;
        for (size_t i = 1; i < r.size(); ++i) { int64_t ga = r[i-1].second, gb = r[i].first; if (gb > ga) gaps.push_back(make_pair(ga, gb)); }
        if (gaps.empty()) continue;
        // foreign intervals spliced between C's first and last appearance
        int lo = span[C].first, hi = span[C].second;
        vector<tuple<step_handle_t, step_handle_t, bool>> foreign;
        for (int j = lo + 1; j < hi; ++j) {
            path_handle_t pj = graph->get_path_handle_of_step(get<0>(intervals[j]));
            if (pj == C) continue;
            if (graph->get_sample_name(pj) != target_sample) foreign.push_back(intervals[j]);
        }
        if (foreign.empty()) continue;   // same-sample interior is handled by the other guard
        // removed = the target sequence skipped over in C's interior
        string removed;
        int64_t pos = 0;
        graph->for_each_step_in_path(C, [&](step_handle_t s) {
            handle_t h = graph->get_handle_of_step(s); int64_t len = graph->get_length(h);
            int64_t na = pos, nb = pos + len; pos = nb;
            for (auto& g : gaps) { int64_t oa = max(na, g.first), ob = min(nb, g.second);
                if (oa < ob) removed += graph->get_sequence(h).substr(oa - na, ob - oa); }
        });
        int64_t nonN = 0; for (char c : removed) if (c != 'N' && c != 'n') ++nonN;
        // nothing left to judge by content: the grafts still in place replaced only N's (assembly gaps),
        // which excise_bad_interior_grafts judged by flank anchoring.  the
        // excised count below would otherwise lift this over min_replaced, and with no k-mers to sample
        // kmer_recovery is 0.0, so one excised graft would revert the whole contig, good patches and all
        {
            unordered_set<uint64_t> removed_kmers;
            sample_kmers(removed, removed_kmers);
            if (removed_kmers.empty()) continue;
        }
        // count non-N already excised from this contig by the partial-patch pass toward the gate, so
        // excising one graft doesn't drop the remainder below the threshold and disable this backstop
        // for a second still-bad graft (a regression vs. the pre-excise lumped full-revert)
        auto eit = excised_nonN.find(C); if (eit != excised_nonN.end()) nonN += eit->second;
        if (nonN < min_replaced) continue;   // too little real sequence replaced to judge reliably
        string added = intervals_to_sequence(graph, foreign);
        double rec = kmer_recovery(removed, added);
        if (rec < min_recovery) {
            ostringstream rs; rs << fixed << setprecision(1) << rec;
            detail = "contig " + graph->get_path_name(C) + " interior (" + to_string(removed.size())
                   + "bp) was replaced by a foreign graft sharing only " + rs.str() + "% of its k-mers";
            return true;
        }
    }
    return false;
}

// Telomere-preservation guard.  A scaffold or graft must never discard a real telomere: if a
// target contig is capped at a natural end but the patch uses that contig starting (or ending)
// well past the cap, the contig was already complete there and the join is redundant/erroneous.
// (Surfaced by the no-CHM13 self-reference comparison: CHM13's divergent subtelomere made panpatch
// trim a telomere-bearing tip to bolt on an overlapping same-haplotype fragment -- HG01074 chr14,
// where the main contig alone was already T2T.)  Legitimate operations are unaffected: telomere
// patches trim a *capless* tip, gap-fills trim only the interior, end-to-end scaffolds use each
// contig in full.  Restricted to the target sample's own contigs.
static bool discards_target_telomere(const PathHandleGraph* graph,
                                     const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals,
                                     const string& target_sample, double telo_threshold,
                                     string& detail) {
    const int64_t OUTER = 20000;   // telomere window at a contig's natural end
    unordered_map<path_handle_t, unordered_map<step_handle_t, int64_t>> idx;
    for (const auto& iv : intervals) {
        path_handle_t p = graph->get_path_handle_of_step(get<0>(iv));
        if (graph->get_sample_name(p) == target_sample) idx[p] = {};
    }
    for (auto& pm : idx) {
        int64_t pos = 0;
        graph->for_each_step_in_path(pm.first, [&](step_handle_t s) { pm.second[s] = pos; pos += graph->get_length(graph->get_handle_of_step(s)); });
    }
    unordered_map<path_handle_t, pair<int64_t, int64_t>> used;   // min start, max end (forward contig coords)
    unordered_map<path_handle_t, int64_t> length;
    for (const auto& iv : intervals) {
        path_handle_t p = graph->get_path_handle_of_step(get<0>(iv));
        if (!idx.count(p)) continue;
        int64_t a, b;
        if (!get<2>(iv)) { a = idx[p][get<0>(iv)]; b = idx[p][get<1>(iv)] + graph->get_length(graph->get_handle_of_step(get<1>(iv))); }
        else             { a = idx[p][get<1>(iv)]; b = idx[p][get<0>(iv)] + graph->get_length(graph->get_handle_of_step(get<0>(iv))); }
        auto it = used.find(p);
        if (it == used.end()) used[p] = make_pair(a, b);
        else { it->second.first = min(it->second.first, a); it->second.second = max(it->second.second, b); }
    }
    for (auto& kv : idx) {
        int64_t len = 0;
        graph->for_each_step_in_path(kv.first, [&](step_handle_t s) { len += graph->get_length(graph->get_handle_of_step(s)); });
        length[kv.first] = len;
    }
    // contigs in output order, so the one reported is the first, not whichever the hash table gives
    for (const path_handle_t& C : paths_in_order(graph, intervals)) {
        if (!used.count(C)) continue;   // not one of the target's own contigs
        int64_t lo = used.at(C).first, hi = used.at(C).second, len = length[C];
        // any trim of an end counts: a telomere array is typically 5-15 kb, so a trim under OUTER can take
        // all of it (this used to look only at trims of OUTER or more)
        if (lo > 0) {                            // 5' end trimmed -> does the natural 5' carry a telomere?
            string tip;
            for (step_handle_t s = graph->path_begin(C); ; s = graph->get_next_step(s)) {
                tip += graph->get_sequence(graph->get_handle_of_step(s));
                if ((int64_t)tip.size() >= OUTER || s == graph->path_back(C)) break;
            }
            if (seq_has_telomere(tip, 0, (int64_t)tip.size(), true, false, telo_threshold)) {
                detail = "patch trimmed the telomere-bearing 5' end of " + graph->get_path_name(C)
                       + " (used from " + to_string(lo) + "bp)";
                return true;
            }
        }
        if (len - hi > 0) {                      // 3' end trimmed -> does the natural 3' carry a telomere?
            string tip;
            for (step_handle_t s = graph->path_back(C); ; s = graph->get_previous_step(s)) {
                tip = graph->get_sequence(graph->get_handle_of_step(s)) + tip;
                if ((int64_t)tip.size() >= OUTER || s == graph->path_begin(C)) break;
            }
            if (seq_has_telomere(tip, 0, (int64_t)tip.size(), true, true, telo_threshold)) {
                detail = "patch trimmed the telomere-bearing 3' end of " + graph->get_path_name(C)
                       + " (used up to " + to_string(hi) + "bp of " + to_string(len) + ")";
                return true;
            }
        }
    }
    return false;
}

bool patch_reuses_sequence(const PathHandleGraph* graph,
                           const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals,
                           const string& target_sample,
                           const vector<path_handle_t>& tgt_paths,
                           const unordered_set<path_handle_t>& threaded_tgt_paths,
                           string& detail) {
    // the forward [lo, hi) range of interval k: intervals are half-open, all but the last, which is closed
    auto range = [&](size_t k) {
        const auto& iv = intervals[k];
        bool last = k + 1 == intervals.size();
        handle_t h0 = graph->get_handle_of_step(get<0>(iv)), h1 = graph->get_handle_of_step(get<1>(iv));
        int64_t o0 = step_offset(graph, get<0>(iv)), o1 = step_offset(graph, get<1>(iv));
        if (!get<2>(iv)) return make_pair(o0, o1 + (last ? (int64_t)graph->get_length(h1) : 0));
        return make_pair(o1 + (last ? 0 : (int64_t)graph->get_length(h1)), o0 + (int64_t)graph->get_length(h0));
    };
    // each path's non-empty ranges, in patch order: (lo, hi, interval index)
    unordered_map<path_handle_t, vector<tuple<int64_t, int64_t, size_t>>> by_path;
    for (size_t k = 0; k < intervals.size(); ++k) {
        pair<int64_t, int64_t> r = range(k);
        if (r.first < r.second) by_path[graph->get_path_handle_of_step(get<0>(intervals[k]))].push_back(make_tuple(r.first, r.second, k));
    }
    for (const path_handle_t& p : paths_in_order(graph, intervals)) {
        auto it = by_path.find(p);
        if (it == by_path.end()) continue;
        const auto& v = it->second;
        if (graph->get_sample_name(p) == target_sample) {
            bool strand = get<2>(intervals[get<2>(v[0])]);
            for (size_t k = 1; k < v.size(); ++k) {
                if (get<2>(intervals[get<2>(v[k])]) != strand) {
                    detail = "the patch uses target contig " + graph->get_path_name(p) + " on both strands";
                    return true;
                }
                // in patch order, the pieces must go along the contig: forward, each from where the last ended
                if (strand ? get<1>(v[k]) > get<0>(v[k - 1]) : get<0>(v[k]) < get<1>(v[k - 1])) {
                    detail = "the patch uses target contig " + graph->get_path_name(p) + " out of its own order (at " +
                        to_string(get<0>(v[k])) + "bp)";
                    return true;
                }
            }
        }
        vector<tuple<int64_t, int64_t, size_t>> sorted_v = v;
        sort(sorted_v.begin(), sorted_v.end());
        for (size_t k = 1; k < sorted_v.size(); ++k) {
            if (get<0>(sorted_v[k]) < get<1>(sorted_v[k - 1])) {
                detail = "the patch uses " + to_string(min(get<1>(sorted_v[k - 1]), get<1>(sorted_v[k])) - get<0>(sorted_v[k])) +
                    "bp of " + graph->get_path_name(p) + " twice (from " + to_string(get<0>(sorted_v[k])) + "bp)";
                return true;
            }
        }
    }
    // target contigs the patch does not use are written through as they are: the patch's donor sequence
    // must not stand in for them too.  a node of such a contig counts if a donor interval of the patch
    // visits it and none of the patch's own target contigs does (so a fragment that merely duplicates part
    // of the patched contig, already a duplicate in the input, does not count)
    unordered_set<path_handle_t> in_patch;
    for (const auto& iv : intervals) in_patch.insert(graph->get_path_handle_of_step(get<0>(iv)));
    vector<path_handle_t> unused;
    for (const path_handle_t& p : tgt_paths) {
        if (graph->get_sample_name(p) == target_sample && !in_patch.count(p) && !threaded_tgt_paths.count(p)) unused.push_back(p);
    }
    if (unused.empty()) return false;
    unordered_set<nid_t> donor_nodes;
    for (size_t k = 0; k < intervals.size(); ++k) {
        const auto& iv = intervals[k];
        if (graph->get_sample_name(graph->get_path_handle_of_step(get<0>(iv))) == target_sample) continue;
        for (step_handle_t s = get<0>(iv); ; s = get<2>(iv) ? graph->get_previous_step(s) : graph->get_next_step(s)) {
            donor_nodes.insert(graph->get_id(graph->get_handle_of_step(s)));
            if (s == get<1>(iv)) break;
        }
    }
    const int64_t MIN_SHARED = 1000;
    for (const path_handle_t& X : unused) {
        int64_t shared = 0;
        graph->for_each_step_in_path(X, [&](step_handle_t s) {
            handle_t h = graph->get_handle_of_step(s);
            if (!donor_nodes.count(graph->get_id(h))) return;
            bool on_tgt = false;
            graph->for_each_step_on_handle(h, [&](step_handle_t o) {
                path_handle_t op = graph->get_path_handle_of_step(o);
                if (op != X && in_patch.count(op) && graph->get_sample_name(op) == target_sample) { on_tgt = true; return false; }
                return true;
            });
            if (!on_tgt) shared += graph->get_length(h);
        });
        if (shared >= MIN_SHARED) {
            detail = "the patch's donor sequence covers " + to_string(shared) + "bp of target contig " + graph->get_path_name(X) +
                ", which it does not use (and which is written through as it is)";
            return true;
        }
    }
    return false;
}

void report_unfilled_gaps(const PathHandleGraph* graph,
                          const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals,
                          const string& target_sample) {
    // each target contig's forward ranges in the output: half-open intervals, all but the last, which is closed
    unordered_map<path_handle_t, vector<pair<int64_t, int64_t>>> ranges;
    for (size_t k = 0; k < intervals.size(); ++k) {
        const auto& iv = intervals[k];
        path_handle_t p = graph->get_path_handle_of_step(get<0>(iv));
        if (graph->get_sample_name(p) != target_sample) continue;
        bool last = k + 1 == intervals.size();
        int64_t o0 = step_offset(graph, get<0>(iv)), o1 = step_offset(graph, get<1>(iv));
        int64_t l0 = graph->get_length(graph->get_handle_of_step(get<0>(iv)));
        int64_t l1 = graph->get_length(graph->get_handle_of_step(get<1>(iv)));
        pair<int64_t, int64_t> r = !get<2>(iv) ? make_pair(o0, o1 + (last ? l1 : 0))
                                               : make_pair(o1 + (last ? 0 : l1), o0 + l0);
        if (r.first < r.second) ranges[p].push_back(r);
    }
    for (const path_handle_t& P : paths_in_order(graph, intervals)) {
        auto rit = ranges.find(P);
        if (rit == ranges.end()) continue;
        string name = graph->get_path_name(P);
        int64_t L = path_bp(graph, P);
        // the contig's N runs
        vector<pair<int64_t, int64_t>> gaps;
        int64_t off = 0, run_start = -1;
        graph->for_each_step_in_path(P, [&](step_handle_t s) {
            string seq = graph->get_sequence(graph->get_handle_of_step(s));
            for (size_t c = 0; c < seq.size(); ++c) {
                bool n = seq[c] == 'N' || seq[c] == 'n';
                if (n && run_start < 0) run_start = off + (int64_t)c;
                if (!n && run_start >= 0) { gaps.push_back(make_pair(run_start, off + (int64_t)c)); run_start = -1; }
            }
            off += seq.size();
        });
        if (run_start >= 0) gaps.push_back(make_pair(run_start, off));
        for (const auto& g : gaps) {
            bool in_output = false;
            for (const auto& r : rit->second) if (r.first < g.second && g.first < r.second) { in_output = true; break; }
            if (!in_output) continue;
            bool explained = false;
            for (const PatchRecord& pr : g_patch_records) {
                if (pr.target == name && pr.target_start >= 0 && pr.target_start < g.second && g.first < pr.target_end) {
                    explained = true; break;
                }
            }
            if (explained) continue;
            PatchRecord pr;
            pr.type = "gap-fill"; pr.target = name; pr.target_bp = L;
            pr.accepted = false;
            pr.reason = "no donor could be threaded across this gap; kept as it is";
            pr.target_start = g.first; pr.target_end = g.second;
            g_patch_records.push_back(pr);
        }
    }
}

// every node id a path traverses (used to test donor/target homology by shared nodes)
static unordered_set<nid_t> path_node_set(const PathHandleGraph* g, path_handle_t p) {
    unordered_set<nid_t> s;
    g->for_each_step_in_path(p, [&](step_handle_t st) { s.insert(g->get_id(g->get_handle_of_step(st))); });
    return s;
}

// Fraction (0-100) of a `window`-bp walk along a path -- starting at `start`, stepping `nxt` (next) or
// !nxt (previous) -- that lands on nodes in `donor_nodes`.  Measuring a fraction over a *fixed window*
// (not a contiguous run) is tolerant of SNP bubbles: a SNP splits one homologous node into a small
// bubble, costing only a few bp.  A faithful fill stays homologous to the target's own flank (high
// fraction); a repeat-region misjoin diverges into different sequence (low fraction).  Crucially this
// works even for N-gap fills, where k-mer recovery cannot judge (the replaced region has no sequence).
static double flank_fraction(const PathHandleGraph* g, step_handle_t start, bool nxt,
                             const unordered_set<nid_t>& donor_nodes, int64_t window) {
    // step off the shared junction node first: it is trivially shared by both paths and would inflate
    // the fraction for short / contig-end flanks
    if (nxt ? !g->has_next_step(start) : !g->has_previous_step(start)) return 0.0;
    step_handle_t s = nxt ? g->get_next_step(start) : g->get_previous_step(start);
    int64_t walked = 0, shared = 0;
    while (walked < window) {
        handle_t h = g->get_handle_of_step(s);
        int64_t use = min((int64_t)g->get_length(h), window - walked);
        walked += use;
        if (donor_nodes.count(g->get_id(h))) shared += use;
        if (nxt ? !g->has_next_step(s) : !g->has_previous_step(s)) break;
        s = nxt ? g->get_next_step(s) : g->get_previous_step(s);
    }
    return walked ? 100.0 * shared / walked : 0.0;
}

// Record scaffold joins (where the patch spans more than one of the target's own contigs) as report
// rows, and guard foreign-bridged joins.  A join [contig A | foreign run | contig B] is a foreign
// bridge: the bridging donor must anchor to A's and B's own sequence over flank_window bp (>= min_flank
// %) on both sides, else it is a wrong-locus misjoin and the whole scaffold must be reverted.  A join
// where A and B are simply adjacent (no foreign bridge) is a plain target-to-target scaffold and is
// recorded but not flank-guarded (telomere validation, if requested, gates it instead).  Returns true
// (+ bridge_detail) if any foreign bridge fails the flank check, so the caller can revert.
bool record_scaffolds(const PathHandleGraph* graph,
                      const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals,
                      const string& target_sample, double min_flank, int64_t flank_window,
                      string& bridge_detail, int64_t max_graft_pieces, const NoPatch* no_patch) {
    bool bad = false;
    unordered_set<path_handle_t> seen;   // target contigs already recorded (so [A|B|A] splices give one row)
    int prev_tgt = -1;       // index of the most recent target-sample interval
    int foreign_start = -1;  // index where a foreign run began since prev_tgt (-1 = none)
    for (int i = 0; i < (int)intervals.size(); ++i) {
        path_handle_t p = graph->get_path_handle_of_step(get<0>(intervals[i]));
        if (graph->get_sample_name(p) != target_sample) { if (foreign_start < 0) foreign_start = i; continue; }
        if (prev_tgt >= 0) {
            path_handle_t pc = graph->get_path_handle_of_step(get<0>(intervals[prev_tgt]));
            if (p != pc && !seen.count(p)) {   // a join into a not-yet-seen target contig
                PatchRecord pr; pr.type = "scaffold";
                pr.target = graph->get_path_name(pc); pr.target_bp = path_bp(graph, pc);
                if (foreign_start >= 0 && foreign_start < i) {
                    // foreign-bridged: flank-check the bridge against both contigs' own sequence
                    // each flank against the donor next to it (a bridge can take more than one donor)
                    path_handle_t donor = graph->get_path_handle_of_step(get<0>(intervals[foreign_start]));
                    path_handle_t donor_r = graph->get_path_handle_of_step(get<0>(intervals[i - 1]));
                    unordered_set<nid_t> dn = path_node_set(graph, donor);
                    unordered_set<nid_t> dn_r = donor_r == donor ? dn : path_node_set(graph, donor_r);
                    double fl = flank_fraction(graph, get<1>(intervals[prev_tgt]), get<2>(intervals[prev_tgt]), dn, flank_window);
                    double fr = flank_fraction(graph, get<0>(intervals[i]), !get<2>(intervals[i]), dn_r, flank_window);
                    pr.donor = graph->get_path_name(donor); pr.donor_bp = path_bp(graph, donor);
                    pr.flankL = fl; pr.flankR = fr;
                    pr.accepted = (fl >= min_flank && fr >= min_flank);
                    if (!pr.accepted) {
                        ostringstream rs; rs << fixed << setprecision(1)
                           << "foreign bridge anchoring " << fl << "%/" << fr << "% < " << min_flank << "% (wrong-locus join)";
                        pr.reason = rs.str();
                        bridge_detail = "contig " + pr.target + " scaffolded to " + graph->get_path_name(p)
                                      + " by foreign " + pr.donor + ": " + rs.str();
                        bad = true;
                    } else if (max_graft_pieces > 0 && i - foreign_start > max_graft_pieces) {
                        // a bridge spliced from several pieces of donor sequence, as for an interior graft
                        pr.accepted = false;
                        pr.reason = "foreign bridge spliced from " + to_string(i - foreign_start) +
                            " pieces of donor sequence, over --max-graft-pieces " + to_string(max_graft_pieces);
                        bridge_detail = "contig " + pr.target + " scaffolded to " + graph->get_path_name(p) + ": " + pr.reason;
                        bad = true;
                    } else if (no_patch && !no_patch->kmers.empty()) {
                        // a bridge of sequence --no-patch-fasta names, as for an interior graft
                        vector<tuple<step_handle_t, step_handle_t, bool>> bridge(intervals.begin() + foreign_start, intervals.begin() + i);
                        double np = no_patch_fraction(intervals_to_sequence(graph, bridge), no_patch);
                        if (np >= 0 && np >= no_patch->fraction) {
                            ostringstream rs; rs << fixed << setprecision(1) << "foreign bridge shares " << 100.0 * np
                                                 << "% of its k-mers with the no-patch sequences (rDNA, --no-patch-fasta)";
                            pr.accepted = false;
                            pr.reason = rs.str();
                            bridge_detail = "contig " + pr.target + " scaffolded to " + graph->get_path_name(p) + ": " + pr.reason;
                            bad = true;
                        }
                    }
                } else {   // adjacent target-to-target scaffold (no foreign bridge)
                    pr.donor = graph->get_path_name(p); pr.donor_bp = path_bp(graph, p);
                }
                g_patch_records.push_back(pr);
            }
        }
        seen.insert(p);
        prev_tgt = i;
        foreign_start = -1;
    }
    return bad;
}

// Partial-patch cleanup (run before revert_bad_patch).  For each foreign interior graft of the shape
// [C-piece | foreign run | same-C-piece], excise it -- drop the foreign run and merge the two flanking
// C-pieces, splicing C's own sequence back in -- when it is a repeat-region misjoin by either test:
//   * k-mer content: it replaced >= min_replaced non-N bp of C but shares < min_recovery of C's k-mers
//     (the donor came from a different locus/array); or
//   * flank anchoring: the donor is not homologous to C's own sequence over flank_window bp on one of
//     the two flanks (< min_flank %), i.e. the graft is not anchored at the right locus.  This catches
//     N-gap fills (which k-mer cannot judge) and any join that rides a long repeat before diverging.
// The rest of the patch (telomere completions, faithful fills) is kept instead of reverting the whole
// contig.  Anything that doesn't fit this clean shape is left to revert_bad_patch's full-revert guard.
void excise_bad_interior_grafts(const PathHandleGraph* graph,
                                vector<tuple<step_handle_t, step_handle_t, bool>>& intervals,
                                const string& target_sample,
                                double min_recovery, int64_t min_replaced,
                                double min_flank, int64_t flank_window,
                                bool patch_gap,
                                unordered_map<path_handle_t, int64_t>& excised_nonN,
                                int64_t max_size_change, double min_complexity, int64_t complexity_window,
                                int64_t complexity_min_bp, int64_t max_graft_pieces, int64_t max_swap_size_change,
                                const NoPatch* no_patch) {
    unordered_map<path_handle_t, unordered_map<step_handle_t, int64_t>> posidx;   // step -> forward pos, per contig (path is stable)
    auto pos_of = [&](path_handle_t C) -> unordered_map<step_handle_t, int64_t>& {
        auto it = posidx.find(C);
        if (it != posidx.end()) return it->second;
        auto& m = posidx[C]; int64_t p = 0;
        graph->for_each_step_in_path(C, [&](step_handle_t s) { m[s] = p; p += graph->get_length(graph->get_handle_of_step(s)); });
        return m;
    };
    vector<PatchRecord> kept_this_call;   // kept grafts from the current pass (final pass survives)
    bool changed = true;
    while (changed) {
        changed = false;
        kept_this_call.clear();
        for (int i = 0; i + 1 < (int)intervals.size() && !changed; ++i) {
            path_handle_t C = graph->get_path_handle_of_step(get<0>(intervals[i]));
            if (graph->get_sample_name(C) != target_sample) continue;
            int j = -1; bool only_foreign = true;
            for (int k = i + 1; k < (int)intervals.size(); ++k) {
                path_handle_t pk = graph->get_path_handle_of_step(get<0>(intervals[k]));
                if (pk == C) { j = k; break; }
                if (graph->get_sample_name(pk) == target_sample) { only_foreign = false; break; }   // a different target contig: not a simple foreign graft
            }
            if (j < 0 || !only_foreign || j == i + 1) continue;
            if (get<2>(intervals[i]) != get<2>(intervals[j])) continue;   // inconsistent orientation: leave to the full-revert guard
            auto& pos = pos_of(C);
            // forward [start,end) range of an interval, orientation-aware: for a reverse interval the
            // convention is swapped (intervals_to_sequence walks get<0> backward to get<1>), so get<1>
            // is the lower-position step and get<0> the higher.
            auto fwd_range = [&](const tuple<step_handle_t, step_handle_t, bool>& iv) {
                int64_t p0 = pos[get<0>(iv)], l0 = graph->get_length(graph->get_handle_of_step(get<0>(iv)));
                int64_t p1 = pos[get<1>(iv)], l1 = graph->get_length(graph->get_handle_of_step(get<1>(iv)));
                return get<2>(iv) ? make_pair(p1, p0 + l0) : make_pair(p0, p1 + l1);
            };
            pair<int64_t, int64_t> rgi = fwd_range(intervals[i]), rgj = fwd_range(intervals[j]);
            int64_t lo, hi;
            // The merge below is (get<0> of i, get<1> of j, orient of i): valid ONLY when the patch
            // visits the two pieces in their on-contig order for that orientation -- forward => piece i
            // is lower, reverse => piece i is higher (piece j lower).  A transposed revisit (order vs.
            // orientation disagree) would make intervals_to_sequence walk off the path end (crash), so
            // bail and leave it to the full-revert backstop.
            if (!get<2>(intervals[i])) {                     // forward: piece i must be the lower piece
                if (rgi.second > rgj.first) continue;
                lo = rgi.second; hi = rgj.first;
            } else {                                         // reverse: piece j must be the lower piece
                if (rgj.second > rgi.first) continue;
                lo = rgj.second; hi = rgi.first;
            }
            if (hi < lo) continue;                           // pieces overlap: not a clean graft, leave it to the backstop
            // hi == lo: the pieces are adjacent, so nothing of C was replaced.  the donor run between them
            // either only passes through the anchor nodes it shares with C (handled below: folded back
            // into C), or inserts sequence C lacks there, which is judged and reported like any graft
            string removed; int64_t p = 0;
            graph->for_each_step_in_path(C, [&](step_handle_t s) {
                handle_t h = graph->get_handle_of_step(s); int64_t l = graph->get_length(h);
                int64_t oa = max(p, lo), ob = min(p + l, hi); if (oa < ob) removed += graph->get_sequence(h).substr(oa - p, ob - oa); p += l; });
            int64_t nonN = 0; for (char c : removed) if (c != 'N' && c != 'n') ++nonN;

            // Metrics for the report: k-mer recovery (CONTENT, when real sequence was replaced) AND flank
            // anchoring (LOCUS, both sides) -- compute both so the table shows them.  The decision uses
            // whichever is authoritative: k-mer for a non-N replacement, flank for a mostly-N gap fill.
            double rec = -1.0, fl = -1.0, fr = -1.0;
            vector<tuple<step_handle_t, step_handle_t, bool>> foreign(intervals.begin() + i + 1, intervals.begin() + j);
            string added = intervals_to_sequence(graph, foreign);
            // the foreign run starts and ends on the anchor nodes it shares with C's pieces either side,
            // which `removed` (C between its pieces) leaves out, so leave them out of what was inserted too
            int64_t trim_front = 0, trim_back = 0;
            if (graph->get_id(graph->get_handle_of_step(get<0>(foreign.front()))) == graph->get_id(graph->get_handle_of_step(get<1>(intervals[i]))))
                trim_front = graph->get_length(graph->get_handle_of_step(get<0>(foreign.front())));
            if (graph->get_id(graph->get_handle_of_step(get<1>(foreign.back()))) == graph->get_id(graph->get_handle_of_step(get<0>(intervals[j]))))
                trim_back = graph->get_length(graph->get_handle_of_step(get<1>(foreign.back())));
            string inserted = (trim_front + trim_back < (int64_t)added.size())
                ? added.substr(trim_front, added.size() - trim_front - trim_back) : string();
            // (measured on the inserted sequence alone: the shared anchor nodes would count as recovered)
            if (nonN >= min_replaced) rec = kmer_recovery(removed, inserted);
            {
                unordered_set<nid_t> dL = path_node_set(graph, graph->get_path_handle_of_step(get<0>(intervals[i + 1])));
                fl = flank_fraction(graph, get<1>(intervals[i]), get<2>(intervals[i]), dL, flank_window);
                unordered_set<nid_t> dR = path_node_set(graph, graph->get_path_handle_of_step(get<0>(intervals[j - 1])));
                fr = flank_fraction(graph, get<0>(intervals[j]), !get<2>(intervals[j]), dR, flank_window);
            }
            bool bad = (nonN >= min_replaced) ? (rec < min_recovery) : (fl < min_flank || fr < min_flank);
            if (!patch_gap) bad = true;   // gap patching disabled (--patch-types): excise every gap-fill graft, restoring the target's own (gapped) sequence

            // grafts that neither test can vouch for, excised to stay conservative: one that makes the
            // contig much longer or shorter (a satellite array or rDNA gap swapped for the donor's), and one
            // that replaces or inserts a long low-complexity run, whose k-mers any array of the same
            // satellite family shares (k-mer recovery passes whether or not the copy is faithful).  the run
            // is looked for in complexity_window-bp windows, so an array is caught inside a larger swap
            if (removed.empty() && inserted.empty()) {
                // a donor run through nothing but the anchor nodes it shares with C: no patch at all (it
                // arises where C's --error-bed interval starts or ends at a shared node), so it is folded
                // back into C, without a report row, as if C had been threaded straight through
                tuple<step_handle_t, step_handle_t, bool> merged =
                    make_tuple(get<0>(intervals[i]), get<1>(intervals[j]), get<2>(intervals[i]));
                intervals.erase(intervals.begin() + i, intervals.begin() + j + 1);
                intervals.insert(intervals.begin() + i, merged);
                changed = true;
                continue;
            }
            int64_t size_change = (int64_t)inserted.size() - (int64_t)removed.size();
            bool too_big = max_size_change > 0 && llabs(size_change) > max_size_change;
            // a swap -- a graft over sequence with no N gap in it, or a pure insertion -- that changes the length
            // is a copy-number change (an array, a tandem duplication), which k-mer recovery cannot check.  an
            // N gap's length is only a placeholder, so a gap fill is held to max_size_change alone
            bool swap_too_big = max_swap_size_change > 0 && nonN == (int64_t)removed.size() &&
                llabs(size_change) > max_swap_size_change;
            // complexity of the replaced and inserted sequence, over the whole of each (cx_) and in its
            // lowest window (cxw_); either one low is enough
            double cx_removed = -1.0, cx_added = -1.0, cxw_removed = -1.0, cxw_added = -1.0;
            int64_t added_nonN = 0; for (char c : inserted) if (c != 'N' && c != 'n') ++added_nonN;
            // judged from complexity_min_bp, well below the k-mer test's min_replaced: a satellite swap of a
            // few kb is no easier to verify than a long one
            if (min_complexity > 0) {
                if (nonN >= complexity_min_bp) {
                    cx_removed = kmer_distinct_fraction(removed);
                    cxw_removed = min_window_distinct_fraction(removed, complexity_window);
                }
                if (added_nonN >= complexity_min_bp) {
                    cx_added = kmer_distinct_fraction(inserted);
                    cxw_added = min_window_distinct_fraction(inserted, complexity_window);
                }
            }
            auto low = [&](double cx) { return cx >= 0 && cx < min_complexity; };
            bool low_cx = low(cx_removed) || low(cx_added) || low(cxw_removed) || low(cxw_added);
            // reported as unverifiable only when the content/locus tests above passed it, so a graft they
            // reject keeps its own reason
            // a graft spliced from more than one piece of donor sequence (pieces of two donor contigs, or of
            // one contig on both strands or out of its order) spells sequence that no donor has in one
            // piece, and that neither test above can vouch for: each piece is anchored where it switches
            int64_t n_pieces = j - i - 1;
            bool chimeric = max_graft_pieces > 0 && n_pieces > max_graft_pieces;
            // a donor run that spells exactly the target sequence it replaces changes nothing (judged before the
            // guards below, which would otherwise count a no-op as a rejected patch).  in a graph
            // built from unmasked assemblies, a donor that agrees with the target across an --error-bed
            // interval threads the target's own nodes there: drop it, so the report counts only patches
            // that change the assembly
            // removed is on C's forward strand and inserted on the patch's: for a C threaded in reverse,
            // compare it with inserted's reverse complement
            auto comp = [](char c) {
                switch (toupper(c)) { case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A'; default: return (char)toupper(c); }
            };
            bool same_seq = !bad && removed.size() == inserted.size() &&
                (get<2>(intervals[i])
                 ? equal(removed.begin(), removed.end(), inserted.rbegin(),
                         [&](char a, char b) { return toupper(a) == comp(b); })
                 : equal(removed.begin(), removed.end(), inserted.begin(),
                         [](char a, char b) { return toupper(a) == toupper(b); }));
            // no-patch sequence (the rDNA unit, --no-patch-fasta): high-copy sequence whose k-mers any copy shares,
            // so k-mer recovery passes whichever copy is swapped in
            double np_removed = no_patch_fraction(removed, no_patch), np_inserted = no_patch_fraction(inserted, no_patch);
            bool no_patch_hit = (np_removed >= 0 && np_removed >= no_patch->fraction) ||
                                (np_inserted >= 0 && np_inserted >= no_patch->fraction);
            bool unverifiable = !bad && !same_seq && (no_patch_hit || too_big || swap_too_big || low_cx || chimeric);
            if (unverifiable || same_seq) bad = true;

            // record this graft for the report.  A kept graft goes to kept_this_call (the loop rescans
            // from the top after every excision, so only the final settled pass's kept grafts should be
            // reported -- avoids duplicates).  An excised graft is recorded now, since it is removed below
            // and will not be rescanned.
            PatchRecord pr;
            pr.type = "gap-fill"; pr.target = graph->get_path_name(C); pr.target_bp = path_bp(graph, C);
            path_handle_t dp = graph->get_path_handle_of_step(get<0>(intervals[i + 1]));
            pr.donor = graph->get_path_name(dp); pr.donor_bp = path_bp(graph, dp);
            pr.replaced_bp = removed.size(); pr.kmer = rec; pr.flankL = fl; pr.flankR = fr;
            pr.target_start = lo; pr.target_end = hi;   // the target region this graft replaced (fwd coords)
            pr.accepted = !bad;
            if (!bad) { kept_this_call.push_back(pr); continue; }   // faithful (content) or anchored (locus) -- keep
            {
                ostringstream rs; rs << fixed << setprecision(1);
                if (!patch_gap)    rs << "gap-fill excluded by --patch-types";
                else if (same_seq) rs << "donor sequence is identical to the target's: nothing to patch";
                else if (unverifiable && no_patch_hit) rs << "the " << (np_removed >= np_inserted ? "replaced" : "inserted")
                                      << " sequence shares " << 100.0 * max(np_removed, np_inserted)
                                      << "% of its k-mers with the no-patch sequences (rDNA, --no-patch-fasta) (unverifiable)";
                else if (unverifiable && chimeric) rs << "spliced from " << n_pieces << " pieces of donor sequence, over --max-graft-pieces "
                                      << max_graft_pieces << " (no donor has it in one piece) (unverifiable)";
                else if (unverifiable && too_big)  rs << "changes the contig's length by " << (size_change > 0 ? "+" : "") << size_change
                                      << "bp, over --max-graft-size-change " << max_size_change << " (unverifiable)";
                else if (unverifiable && swap_too_big) rs << "swap changes the contig's length by " << (size_change > 0 ? "+" : "") << size_change
                                      << "bp, over --max-swap-size-change " << max_swap_size_change << " (unverifiable)";
                else if (unverifiable) {
                    // the whole-sequence fraction when it is the low one, else the lowest window's
                    bool whole = low(cx_removed) || low(cx_added);
                    double cx = whole ? (low(cx_removed) ? cx_removed : cx_added)
                                      : (low(cxw_removed) ? cxw_removed : cxw_added);
                    rs << setprecision(2) << "low-complexity sequence (distinct k-mer fraction " << cx;
                    if (!whole) rs << " in a " << complexity_window << "bp window";
                    rs << " replaced/inserted < --min-graft-complexity " << min_complexity << ") (unverifiable)";
                }
                else if (rec >= 0) rs << "k-mer recovery " << rec << "% < " << min_recovery << "% (repeat-region misjoin)";
                else               rs << "flank anchoring " << fl << "%/" << fr << "% < " << min_flank << "% (wrong-locus join)";
                pr.reason = rs.str();
            }
            g_patch_records.push_back(pr);
            if (!same_seq) excised_nonN[C] += nonN;   // so the revert backstop's threshold still sees this replaced bp
            // the patch visits piece i then piece j; merging them spans the whole region in the patch's
            // own direction -- (get<0> of i, get<1> of j) works for forward and reverse alike.
            tuple<step_handle_t, step_handle_t, bool> merged =
                make_tuple(get<0>(intervals[i]), get<1>(intervals[j]), get<2>(intervals[i]));
            intervals.erase(intervals.begin() + i, intervals.begin() + j + 1);
            intervals.insert(intervals.begin() + i, merged);
            changed = true;
        }
    }
    // the final settled pass's kept grafts are the ones actually in the output
    for (auto& kp : kept_this_call) g_patch_records.push_back(kp);
}

bool revert_bad_patch(const PathHandleGraph* graph,
                      const path_handle_t& ref_path,
                      const vector<path_handle_t>& tgt_paths,
                      const vector<string>& sample_names,
                      const vector<tuple<step_handle_t, step_handle_t, bool>>& in_intervals,
                      vector<tuple<step_handle_t, step_handle_t, bool>>& out_intervals,
                      string default_sample,
                      double threshold,
                      double graft_recovery,
                      int64_t graft_min_bp,
                      double telo_threshold,
                      const unordered_map<path_handle_t, int64_t>& excised_nonN,
                      string& revert_reason) {

    out_intervals.clear();
    
    vector<path_handle_t> first_tgt_paths;
    int64_t tgt_length = 0;
    for (const path_handle_t& tgt_path : tgt_paths) {
        if (graph->get_sample_name(tgt_path) == sample_names[0]) {
            graph->for_each_step_in_path(tgt_path, [&](step_handle_t step) {
                tgt_length += graph->get_length(graph->get_handle_of_step(step));
            });
            first_tgt_paths.push_back(tgt_path);
        }
    }

    int64_t patch_length = intervals_to_sequence(graph, in_intervals).length();

    // --min-cover compares the patch with the target contigs it is made from: one it never uses (an
    // unplaced fragment binned on this chromosome) is written through as it is, so counting it in
    // reverted good patches on chromosomes with a large unused fragment.  (all of them when it uses none)
    {
        unordered_set<path_handle_t> in_patch;
        for (const auto& iv : in_intervals) in_patch.insert(graph->get_path_handle_of_step(get<0>(iv)));
        int64_t used_length = 0;
        for (const path_handle_t& tgt_path : first_tgt_paths) {
            if (in_patch.count(tgt_path)) used_length += path_bp(graph, tgt_path);
        }
        if (used_length > 0) tgt_length = used_length;
    }

    // we replace the patch with the input because it was too short
    bool to_revert = (double)patch_length / (double)tgt_length < threshold;
    if (to_revert) {
        ostringstream rr; rr << "patch covers only " << ((double)patch_length / (double)tgt_length) << " of target (< --min-cover)";
        revert_reason = rr.str();
    }

    // sanity guard: reject patches that replaced real interior sequence of a contig without an
    // assembly gap to justify it (repeat-region misjoin -- e.g. fragments spliced into satellite/
    // segdup through ambiguous anchors).  this overrides an otherwise-accepted (even telomere-valid)
    // patch and reverts to the input contigs.
    string interior_detail;
    if (splices_same_sample_interior(graph, in_intervals, sample_names[0], interior_detail)) {
        revert_reason = interior_detail + " (repeat-region misjoin)";
        to_revert = true;
    }

    // reject a foreign interior graft that recapitulates almost none of the target sequence it
    // replaced -- a repeat-region misjoin where the donor came from a different locus
    // (thresholds controlled by --graft-recovery / --graft-min-bp)
    string graft_detail;
    if (interior_graft_low_recovery(graph, in_intervals, sample_names[0], graft_recovery, graft_min_bp, excised_nonN, graft_detail)) {
        revert_reason = graft_detail + " (repeat-region misjoin)";
        to_revert = true;
    }

    // telomere-preservation: a patch must not discard a telomere the target already had
    string telo_detail;
    if (discards_target_telomere(graph, in_intervals, sample_names[0], telo_threshold, telo_detail)) {
        revert_reason = telo_detail + " (target was already capped)";
        to_revert = true;
    }

    // no-patch detection: the threaded output is just the target's own single contig (no foreign
    // sequence spliced in, no second target contig scaffolded) -> nothing was actually patched.
    // this runs even when -e/--default-sample is set: "nothing to patch" means the target contig is
    // emitted untouched (native orientation) below, NOT replaced -- -e substitutes only on a real
    // patch *failure* (the to_revert checks above), never here.  (gating this on !default_sample used
    // to skip detection entirely with -e, leaving the reference-threaded path in place, which
    // reverse-complements any target contig stored antisense to the reference.)
    bool no_patch = false;
    if (!to_revert) {
        bool patch_happened = false;
        unordered_set<path_handle_t> tgt_contigs_in_patch;
        for (const auto& interval : in_intervals) {
            path_handle_t interval_path = graph->get_path_handle_of_step(get<0>(interval));
            if (graph->get_sample_name(interval_path) != sample_names[0]) {
                // a different (non-target) sample contributed sequence: this is a patch
                patch_happened = true;
                break;
            }
            tgt_contigs_in_patch.insert(interval_path);
        }
        // scaffolding two or more of the target sample's own contigs into a single
        // sequence is also a patch, even when no foreign sequence was used to bridge them.
        // (if -T is given, this join is telomere-validated upstream in main, and reverted
        //  there if it fails; without -T there is no check and the join is kept as-is)
        if (tgt_contigs_in_patch.size() > 1) {
            patch_happened = true;
        }
        // we replace the patch with the target's own contig because there was no patch
        if (!patch_happened) {
            to_revert = true;
            no_patch = true;
            // check if target paths have any gaps (N bases)
            bool has_gaps = false;
            for (const path_handle_t& tgt_path : first_tgt_paths) {
                graph->for_each_step_in_path(tgt_path, [&](step_handle_t step) {
                    if (!has_gaps) {
                        string seq = graph->get_sequence(graph->get_handle_of_step(step));
                        for (char c : seq) {
                            if (c == 'N' || c == 'n') {
                                has_gaps = true;
                                return;
                            }
                        }
                    }
                });
                if (has_gaps) break;
            }
            revert_reason = has_gaps ? NO_DONOR_PATCH_REASON
                                     : "no patching required (no gaps)";
        }
    }

    if (to_revert) {
        // -e/--default-sample substitutes its own contig only for a genuine patch *failure*; a
        // no-patch revert always restores the target's own contig in native orientation.
        if (!default_sample.empty() && !no_patch) {
            vector<path_handle_t> default_paths = default_sample_paths(graph, default_sample,
                                                                       graph->get_haplotype(first_tgt_paths.front()));
            // the default sample stands in only where it has exactly one contig here; otherwise (none, or
            // several to choose between) the target's own contigs are kept, rather than aborting the run
            if (default_paths.size() == 1) {
                out_intervals.push_back(make_tuple(graph->path_begin(default_paths[0]), graph->path_back(default_paths[0]), false));
                return true;
            }
            cerr << "[panpatch] warning: -e " << default_sample << " has " << default_paths.size() << " contigs on "
                 << graph->get_locus_name(ref_path) << ", not 1; keeping the target's own contigs" << endl;
        }
        {
            for (const path_handle_t& tgt_path : first_tgt_paths) {
                out_intervals.push_back(make_tuple(graph->path_begin(tgt_path), graph->path_back(tgt_path), false));
            }
        }
        return true;
    }
    return false;
}

void check_intervals(const PathHandleGraph* graph,
                     const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals) {
    for (size_t idx = 0; idx < intervals.size(); ++idx) {
        const auto& interval = intervals[idx];
        path_handle_t interval_path = graph->get_path_handle_of_step(get<0>(interval));
        path_handle_t end_path = graph->get_path_handle_of_step(get<1>(interval));

        // Verify start and end are on the same path
        if (interval_path != end_path) {
            cerr << "ERROR: Interval " << idx << " has mismatched paths!" << endl;
            cerr << "  Start path: " << graph->get_path_name(interval_path) << endl;
            cerr << "  End path: " << graph->get_path_name(end_path) << endl;
            assert(false);
        }

#ifdef debug
        cerr << "Checking interval " << idx << ": " << graph->get_path_name(interval_path)
             << " " << graph->get_id(graph->get_handle_of_step(get<0>(interval))) << "-"
             << graph->get_id(graph->get_handle_of_step(get<1>(interval)))
             << " rev=" << get<2>(interval) << endl;
#endif
        // Note: Validation disabled because smooth_intervals may create intervals that span
        // across what appears to be disconnected segments. The actual sequence extraction
        // in intervals_to_sequence handles the closed-interval semantics correctly.
        // TODO: Investigate why smoothed intervals can have end steps that appear unreachable
    }
}

bool validate_telomeres(const PathHandleGraph* graph,
                        const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals,
                        double threshold,
                        bool require_ends,
                        bool verbose,
                        string& fail_reason) {

    if (intervals.empty()) {
        if (verbose) {
            cerr << "[panpatch] Telomere validation: no intervals to validate" << endl;
        }
        fail_reason = "telomere validation failed";
        return false;
    }

    // Get the full sequence
    string sequence = intervals_to_sequence(graph, intervals);
    int64_t seq_len = sequence.length();

    // Telomere detection uses the shared seq_has_telomere() helper (also used by the patcher, so the
    // two cannot disagree). Tips use boundary detection; the internal check does not.  A sequence too
    // short to analyze leaves all three flags false -- harmless unless require_ends demands both ends.
    bool has_start_telomere = false, has_end_telomere = false, has_internal_telomere = false;
    if (seq_len >= 2000) {
        // Search window - check up to 50kb from each end
        int64_t tip_check_len = min((int64_t)50000, seq_len / 2);
        has_start_telomere = seq_has_telomere(sequence, 0, tip_check_len, true, false, threshold);
        has_end_telomere = seq_has_telomere(sequence, max((int64_t)0, seq_len - tip_check_len), seq_len, true, true, threshold);
        // internal telomeres should NOT exist; don't find boundary here -- detect any telomeric sequence
        if (seq_len > 2 * tip_check_len) {
            has_internal_telomere = seq_has_telomere(sequence, tip_check_len, seq_len - tip_check_len, false, false, threshold);
        }
    }

    // a scaffold join must not bury a telomere: where the patch goes from one target contig to another
    // (directly, or across a donor bridge), neither contig's end at the join, nor the bridge between them,
    // may hold one.  the whole-interior density above cannot see a telomere array there (it would have to
    // be most of the interior), so look at each join, from JOIN_CHECK bp before it to JOIN_CHECK bp after
    if (!has_internal_telomere && intervals.size() > 1) {
        const int64_t JOIN_CHECK = 20000;
        // each interval's [start, end) in the patched sequence: half-open, all but the last, which is closed
        vector<pair<int64_t, int64_t>> out_range;
        int64_t cum = 0;
        for (size_t k = 0; k < intervals.size(); ++k) {
            const auto& iv = intervals[k];
            bool last = k + 1 == intervals.size();
            int64_t o0 = step_offset(graph, get<0>(iv)), o1 = step_offset(graph, get<1>(iv));
            int64_t l0 = graph->get_length(graph->get_handle_of_step(get<0>(iv)));
            int64_t l1 = graph->get_length(graph->get_handle_of_step(get<1>(iv)));
            int64_t len = !get<2>(iv) ? (o1 - o0 + (last ? l1 : 0)) : (o0 + l0 - o1 - (last ? 0 : l1));
            out_range.push_back(make_pair(cum, cum + max((int64_t)0, len)));
            cum += max((int64_t)0, len);
        }
        string target_sample = graph->get_sample_name(graph->get_path_handle_of_step(get<0>(intervals[0])));
        int prev_t = -1;
        for (size_t k = 0; k < intervals.size() && !has_internal_telomere; ++k) {
            path_handle_t p = graph->get_path_handle_of_step(get<0>(intervals[k]));
            if (graph->get_sample_name(p) != target_sample) continue;
            if (prev_t >= 0 && graph->get_path_handle_of_step(get<0>(intervals[prev_t])) != p) {
                int64_t a_end = min(out_range[prev_t].second, seq_len), b_start = min(out_range[k].first, seq_len);
                int64_t from = max((int64_t)0, a_end - JOIN_CHECK), to = min(seq_len, max(a_end, b_start) + JOIN_CHECK);
                if (to > from && has_telomeric_run(sequence.substr(from, to - from), threshold, 1000)) {
                    has_internal_telomere = true;
                    if (verbose) cerr << "[panpatch]   telomere at the join of " << graph->get_path_name(graph->get_path_handle_of_step(get<0>(intervals[prev_t])))
                                      << " and " << graph->get_path_name(p) << endl;
                }
            }
            prev_t = (int)k;
        }
    }

    if (verbose) {
        cerr << "[panpatch] Telomere validation (threshold=" << threshold << ", require_ends=" << require_ends << "):" << endl;
        cerr << "[panpatch]   Sequence length: " << seq_len << "bp" << endl;
        cerr << "[panpatch]   Start telomere: " << (has_start_telomere ? "FOUND" : "NOT FOUND") << endl;
        cerr << "[panpatch]   End telomere: " << (has_end_telomere ? "FOUND" : "NOT FOUND") << endl;
        cerr << "[panpatch]   Internal telomeres: " << (has_internal_telomere ? "FOUND (BAD)" : "NOT FOUND (GOOD)") << endl;
    }

    // an internal telomere is a scaffold misjoin -- always disqualifying (with -T or without)
    if (has_internal_telomere) {
        fail_reason = "internal telomere in patched contig (scaffold misjoin)";
        return false;
    }
    // require_ends (-T) additionally demands a real telomere at BOTH ends; telomere patching alone does not
    if (require_ends && !(has_start_telomere && has_end_telomere)) {
        fail_reason = "telomere validation failed";
        return false;
    }
    return true;
}

void print_contig_telomere_line(const string& name, const string& sequence, double threshold) {
    int64_t seq_len = sequence.length();

    if (seq_len < 50) {
        return;  // Too short to analyze
    }

    // Check tip regions - use up to 50kb search window
    int64_t tip_check_len = min((int64_t)50000, seq_len / 2);

    // use the shared seq_has_telomere (the patcher/validator's detector) so this printed status
    // cannot disagree with the patch decisions; report the density it measured
    double left_max_density = 0.0, right_max_density = 0.0;
    bool has_left  = seq_has_telomere(sequence, 0, tip_check_len, true, false, threshold, &left_max_density);
    bool has_right = seq_has_telomere(sequence, max((int64_t)0, seq_len - tip_check_len), seq_len, true, true, threshold, &right_max_density);

    // Log all contigs with telomere information
    cout << "#Contig " << name
         << " len=" << seq_len << "bp"
         << " left=" << (has_left ? "YES" : "NO")
         << "(" << fixed << setprecision(3) << left_max_density << ")"
         << " right=" << (has_right ? "YES" : "NO")
         << "(" << fixed << setprecision(3) << right_max_density << ")"
         << endl;
}

int telomere_report_fasta(const string& fasta_path, double threshold) {
    // gzopen reads plaintext transparently as well as gzip/bgzip, so this accepts a finished .fa.gz
    // (what a user has on hand) and the uncompressed file (what cactus-panpatch measures pre-bgzip)
    gzFile in = gzopen(fasta_path.c_str(), "rb");
    if (in == Z_NULL) {
        cerr << "[panpatch] error: Unable to open FASTA " << fasta_path << endl;
        return 1;
    }
    string name, sequence;
    auto flush_record = [&]() {
        if (!name.empty()) {
            print_contig_telomere_line(name, sequence, threshold);
        }
        name.clear();
        sequence.clear();
    };
    // a FASTA line is normally <=80bp, but nothing forbids one long line per record, so read in fixed
    // chunks and only treat a '\n' as a line end (a chunk that fills the buffer is a partial line)
    const size_t BUF = 1 << 16;
    vector<char> buf(BUF);
    string line;
    bool ok = true;
    while (gzgets(in, buf.data(), (int)BUF) != Z_NULL) {
        line += buf.data();
        if (line.empty() || line.back() != '\n') {
            continue;                                               // partial line: keep accumulating
        }
        line.pop_back();
        if (!line.empty() && line.back() == '\r') line.pop_back();   // tolerate CRLF
        if (!line.empty()) {
            if (line[0] == '>') {
                flush_record();
                // the record name is the header's first token, matching how contigs are named elsewhere
                size_t ws = line.find_first_of(" \t");
                name = line.substr(1, ws == string::npos ? string::npos : ws - 1);
            } else {
                // normalise case on read: a soft-masked assembly would otherwise never match the
                // (uppercase) telomere hexamers, and telomere presence must not depend on
                // repeat-masking.  Graph sequence is already uppercase, so this only matters here.
                for (char& ch : line) ch = toupper((unsigned char)ch);
                sequence += line;
            }
        }
        line.clear();
    }
    int errnum = 0;
    const char* gzerr = gzerror(in, &errnum);
    if (errnum != Z_OK && errnum != Z_STREAM_END) {
        cerr << "[panpatch] error: Failed reading " << fasta_path << ": " << gzerr << endl;
        ok = false;
    }
    gzclose(in);
    if (!ok) {
        return 1;
    }
    if (!line.empty()) {          // final line with no trailing newline
        if (line[0] != '>') {
            for (char& ch : line) ch = toupper((unsigned char)ch);
            sequence += line;
        }
    }
    flush_record();
    return 0;
}

