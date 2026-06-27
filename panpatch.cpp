#include <unordered_set>
#include <cassert>
#include <map>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include "panpatch.hpp"
#include <limits>
#include <cstdint>

//#define debug
//#define ultra_debug

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

    unordered_map<string, unordered_map<int64_t, double>> sample_hap_coverage;
    unordered_map<string, unordered_map<int64_t, int64_t>> sample_hap_count;
    for (const auto& path_cov : coverage_map) {
        string sample = graph->get_sample_name(path_cov.first);
        int64_t haplotype = graph->get_haplotype(path_cov.first);
        int64_t path_length = 0;
        graph->for_each_step_in_path(path_cov.first, [&](const step_handle_t& step) {
            path_length += graph->get_length(graph->get_handle_of_step(step));
        });        
        sample_hap_coverage[sample][haplotype] += (double)path_length * path_cov.second;
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
    for (const auto& path_cov : coverage_map) {
        string sample = graph->get_sample_name(path_cov.first);
        int64_t haplotype = graph->get_haplotype(path_cov.first);
        if (haplotype == sample_to_hap[sample]) {
            result[sample].push_back(path_cov.first);
        } 
    }

    return result;
}

multimap<pair<int64_t, int64_t>, path_handle_t> sort_overlapping_paths(const PathHandleGraph* graph,
                                                                       const path_handle_t& tgt_path,
                                                                       const vector<path_handle_t>& other_paths) {

    // note: this logic only works properly on acyclic reference path
    unordered_map<int64_t, int64_t> id2pos;
    int64_t pos = 0;
    graph->for_each_step_in_path(tgt_path, [&](step_handle_t step) {
        id2pos[graph->get_id(graph->get_handle_of_step(step))] = pos;
        pos += graph->get_length(graph->get_handle_of_step(step));
    });


    multimap<pair<int64_t, int64_t>, path_handle_t> result;

    for (path_handle_t other_path : other_paths) {
        int64_t min_pos = numeric_limits<int64_t>::max();
        int64_t max_pos = 0;
        graph->for_each_step_in_path(other_path, [&](step_handle_t step) {
            int64_t node_id = graph->get_id(graph->get_handle_of_step(step));
            if (id2pos.count(node_id)) {
                int64_t pos = id2pos[node_id];
                min_pos = min(min_pos, pos);
                max_pos = max(max_pos, pos);
            }
        });

        if (min_pos < max_pos) {
            pair<int64_t, int64_t> key = make_pair(min_pos, max_pos);
            result.insert(make_pair(key, other_path));
        }
    }

    return result;
}

pair<int64_t, int64_t> find_telomeres(const PathHandleGraph* graph,
                                      const path_handle_t path,
                                      double threshold) {

    // quick and dirty telomere checker!!
    
    static const int64_t min_len = 50;
    string path_str;
    graph->for_each_step_in_path(path, [&](step_handle_t step){
        path_str += graph->get_sequence(graph->get_handle_of_step(step));
    });

    // forward
    int64_t fw_count = 0;
    int64_t r_count = 0;
    int64_t pos;
    for (pos = 0; pos < path_str.length()-7; ++pos) {
        if (path_str.substr(pos, 6) == "TTAGGG") {
            ++fw_count;
            pos+= 5;
        } else if (path_str.substr(pos, 6) == "CCCTAA") {
            ++r_count;
            pos+= 5;
        }
        if (pos > min_len) {
            double fw_density = 6. * ((double)fw_count / (double) pos);
            double r_density = 6. * ((double)r_count / (double) pos);
            if (fw_density < threshold && r_density < threshold) {
                break;
            }
        }
    }

    // reverse
    fw_count = 0;
    r_count = 0;
    int64_t r_pos;
    for (r_pos = 0; r_pos < path_str.length()-7; ++r_pos) {
        if (path_str.substr(path_str.length() - 7 - r_pos, 6) == "TTAGGG") {
            ++fw_count;
            r_pos+= 5;
        } else if (path_str.substr(path_str.length() - 7 - r_pos, 6) == "CCCTAA") {
            ++r_count;
            r_pos+= 5;
        }
        if (r_pos > min_len) {
            double fw_density = 6. * ((double)fw_count / (double) r_pos);
            double r_density = 6. * ((double)r_count / (double) r_pos);
            if (fw_density < threshold && r_density < threshold) {
                break;
            }
        }
    }

    return make_pair(pos, r_pos);
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
                                                   int64_t bound) {

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

    // search forward
    int64_t fcount = 0;
    int64_t bcount = 0;
    if (direction >= 0) {
        step_handle_t next_step = graph->get_next_step(step);
        for (; next_step != graph->path_end(path); next_step = graph->get_next_step(next_step)) {
            handle_t next_handle = graph->get_handle_of_step(next_step);
            ++fcount;
            if (!cross_Ns && has_n(next_handle)) {
                break;
            }
            if (anchors.count(graph->get_id(next_handle))) {
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
            if (!cross_Ns && has_n(next_handle)) {
                break;
            }
            if (anchors.count(graph->get_id(next_handle))) {
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
                                                                   const ExcludedRefRegions& excluded_regions) {

    assert(ref_anchors.size() > 1);
    assert(tgt_paths.size() > 0);

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
    bool cur_backward = false;
    assert(cur_backward == graph->get_is_reverse(cur_handle));
    unordered_map<path_handle_t, int64_t> tgt_ranks;
        
    while (true) {
        // find the steps on the handle and sort them using the path priority
        // [paths that don't have a priority are ignored which is important]
        // todo: could be in weird regions we need to assess contiguiuty with
        // previous anchor, but not trying for now. 
        multimap<int64_t, step_handle_t> sorted_steps;
        graph->for_each_step_on_handle(cur_handle, [&](step_handle_t step) {
            path_handle_t path = graph->get_path_handle_of_step(step);
            if (path_rank.count(path)) {
                int64_t rank = tgt_ranks.count(path) ? tgt_ranks.at(path) : path_rank.at(path);
                sorted_steps.insert(make_pair(rank, step));
            }
        });

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
                bool cross_gaps = iteration > 0;
                bool step_backward = graph->get_is_reverse(graph->get_handle_of_step(step));
                int64_t direction = 0;
                if (!interval_cover.empty() && graph->get_path_handle_of_step(get<0>(interval_cover.back())) ==
                    graph->get_path_handle_of_step(step)) {
                    // we constrain the direction if we're already following the path, otherwise it can go either way
                    direction = get<2>(interval_cover.back()) ? -1 : 1;
                }
                pair<step_handle_t, bool> next_anchor = find_next_anchor_on_path(graph, ref_anchors, step, cur_pos,
                                                                                 direction, cross_gaps);
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
                    cur_handle = graph->get_handle_of_step(next_anchor.first);
                    cur_pos = ref_anchors.at(graph->get_id(cur_handle));
                    cur_backward = graph->get_is_reverse(cur_handle);
                    if (path_rank[graph->get_path_handle_of_step(next_anchor.first)] < tgt_ranks.size() &&
                        !tgt_ranks.count(graph->get_path_handle_of_step(next_anchor.first))) {
                        // preference for first-visited target paths when choosing fork
                        tgt_ranks[graph->get_path_handle_of_step(next_anchor.first)] = tgt_ranks.size();
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
                                                                   const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals) {
    vector<tuple<step_handle_t, step_handle_t, bool>> extended_intervals = intervals;

    // extend the front
    step_handle_t first_step = get<0>(extended_intervals[0]);
    bool first_backward = get<2>(extended_intervals[0]);
    if (!first_backward) {
        first_step = graph->path_begin(graph->get_path_handle_of_step(first_step));
    } else {
        first_step = graph->path_end(graph->get_path_handle_of_step(first_step));
        first_step = graph->get_previous_step(first_step);
    }
    extended_intervals[0] = make_tuple(first_step, get<1>(extended_intervals[0]), first_backward);

    // extend the back
    step_handle_t last_step = get<1>(extended_intervals.back());
    bool last_backward = get<2>(extended_intervals.back());
    if (!last_backward) {
        last_step = graph->path_end(graph->get_path_handle_of_step(last_step));
        last_step = graph->get_previous_step(last_step);        
    } else {
        last_step = graph->path_begin(graph->get_path_handle_of_step(last_step));
    }
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
                             bool find_boundary, bool is_right_end, double threshold) {
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
    if (region_len < 500) return false;
    double fw_density = 6. * ((double)fw_count / (double)region_len);
    double r_density = 6. * ((double)r_count / (double)region_len);
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
// extend_intervals() only extends the two terminal intervals along their *own* (target)
// paths to the contigs' own ends.  If the target assembly is simply missing a telomere
// there, nothing can add it.  This routine grafts a missing telomere in from a foreign
// cover (selected donor assembly) when one is available:
//
//   For each terminal end of the assembly that lacks a telomere, walk inward from the tip
//   along the terminal target interval looking for the nearest node also visited by a
//   foreign cover; if that cover continues outward (in the assembly's frame) and ends in a
//   real telomere, hand off to it there.  The target's (divergent, capless) sequence beyond
//   the handoff node is replaced by the foreign cover's run to its telomere.
//
// This is graph-coherent: the handoff is a node genuinely shared by both paths, so the join
// is supported by the graph rather than a blind concatenation.
vector<tuple<step_handle_t, step_handle_t, bool>> extend_to_telomeres(
        const PathHandleGraph* graph,
        const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals,
        const unordered_map<string, vector<path_handle_t>>& sample_covers,
        const vector<string>& sample_names,
        double telo_threshold,
        int64_t max_handoff,
        bool verbose) {

    static const int64_t OUTER = 20000;            // inspect this many bp at a contig tip for a telomere
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
        while (true) {
            handle_t a_i = asm_handle(s_i);
            handle_t under = graph->get_handle_of_step(s_i);
            handle_t out_handle = is_front ? graph->flip(a_i) : a_i;

            // collect the foreign covers present on this node in one pass (first step of each)
            unordered_map<path_handle_t, step_handle_t> steps_here;
            graph->for_each_step_on_handle(under, [&](step_handle_t fs) {
                path_handle_t fp = graph->get_path_handle_of_step(fs);
                if (foreign_set.count(fp) && !steps_here.count(fp)) steps_here[fp] = fs;
            });

            for (const path_handle_t& F : foreign) {
                auto sit = steps_here.find(F);
                if (sit == steps_here.end()) continue;  // cover F does not visit this node
                step_handle_t s_F = sit->second;

                // pick the direction along F that continues outward in the assembly frame
                handle_t gF = graph->get_handle_of_step(s_F);
                bool f_forward;
                if (gF == out_handle) f_forward = true;
                else if (gF == graph->flip(out_handle)) f_forward = false;
                else continue;

                // the outward end of F in this direction must be telomeric (precomputed)
                if (!(f_forward ? telo_end[F] : telo_begin[F])) continue;

                // `walked` is the length of the target's capless tail this handoff would replace.
                bool within_cap = (walked <= max_handoff);

                // follow F outward to its (telomeric) end, bounded; this confirms reachability and
                // gives us the terminal step. accumulate the grafted sequence only when we will
                // actually patch (within cap), for the recovery stat.
                int64_t ext = 0;
                string added_seq;
                step_handle_t fs = f_forward ? graph->get_next_step(s_F) : graph->get_previous_step(s_F);
                step_handle_t f_last = s_F;
                bool reached = false, any = false;
                while (true) {
                    bool at_end = f_forward ? (fs == graph->path_end(F)) : (fs == graph->path_front_end(F));
                    if (at_end) { reached = true; break; }
                    handle_t fh = graph->get_handle_of_step(fs);
                    if (within_cap) added_seq += graph->get_sequence(f_forward ? fh : graph->flip(fh));
                    ext += graph->get_length(fh);
                    f_last = fs; any = true;
                    if (ext > MAX_FOREIGN_EXT) break;
                    fs = f_forward ? graph->get_next_step(fs) : graph->get_previous_step(fs);
                }
                if (!reached || !any) continue;  // not a usable handoff; keep looking

                if (!within_cap) {
                    // nearest usable handoff is too far in: replacing this much target is risky, so skip
                    cout << "#Telomere not patched (" << (is_front ? "front" : "back")
                         << "): nearest donor handoff (" << graph->get_path_name(F) << ") would replace "
                         << walked << "bp of target, over the --max-telomere-patch cap of " << max_handoff
                         << "bp (rerun with -M " << walked << " to allow)" << endl;
                    return false;
                }

                // replaced target tail (tip .. just inside the handoff node) for the recovery stat
                string removed_seq;
                for (step_handle_t s = tip_step; s != s_i; s = inward(s)) {
                    removed_seq += graph->get_sequence(asm_handle(s));
                }
                double recovery = kmer_recovery(removed_seq, added_seq);
                // format the percentage in a local stream so we don't leave cout stuck in
                // fixed/setprecision state (those manipulators are sticky and would corrupt
                // later default-formatted floats, e.g. the revert ratio in revert_bad_patch)
                ostringstream rec_ss;
                rec_ss << fixed << setprecision(1) << recovery;
                cout << "#Telomere patch (" << (is_front ? "front" : "back") << "): donor="
                     << graph->get_path_name(F) << " replaced=" << walked << "bp grafted=" << ext
                     << "bp kmer_recovery=" << rec_ss.str() << "%" << endl;

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

        // capless end with no usable donor handoff in range: report why, so the user can see which
        // ends were left as-is and whether it's a simple gap (no telomere anywhere, no donor) or an
        // assembly issue beyond panpatch's scope (a telomere present but buried/degraded at the tip).
        cout << "#Telomere not patched (" << (is_front ? "front" : "back") << ") of "
             << graph->get_path_name(P) << ": ";
        if (buried) {
            cout << "telomeric repeats are present near the tip but not as a clean terminal telomere "
                    "(sequence extending past the telomere, or a degraded/fragmented telomere) - "
                    "beyond simple patching; no donor provides a clean telomere here" << endl;
        } else {
            cout << "no telomere at this end and no donor assembly reaches one near it" << endl;
        }
        return false;
    };

    step_handle_t s_i_front, s_i_back;
    tuple<step_handle_t, step_handle_t, bool> fi_front, fi_back;
    bool fp = find_end_patch(intervals.front(), true,  s_i_front, fi_front);
    bool bp = find_end_patch(intervals.back(),  false, s_i_back,  fi_back);
    if (!fp && !bp) return intervals;

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
                     const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals) {

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
        cout << graph->get_path_name(path) << "\t" << pos_1 << "\t" << pos_2
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
                                                               bool verbose) {


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

    // find the anchors along the reference path using the relevant paths
    unordered_set<path_handle_t> relevant_paths = {ref_path};
    for (const auto& sample_paths : sample_covers) {
        for (const path_handle_t& path : sample_paths.second) {
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
                                                                                           excluded_regions);

    if (patched_intervals.empty()) {
        cerr << "[panpatch] warning: unable to patch assembly on " << graph->get_locus_name(ref_path)
             << " for " << graph->get_sample_name(tgt_paths.front()) << "#"
         << graph->get_haplotype(tgt_paths.front()) << ":" << endl;
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

    vector<tuple<step_handle_t, step_handle_t, bool>> extended_intervals = extend_intervals(graph, smoothed_intervals);

    check_intervals(graph, extended_intervals);

    if (patch_ends) {
        extended_intervals = extend_to_telomeres(graph, extended_intervals, sample_covers,
                                                 sample_names, telo_threshold, max_telomere_patch, verbose);
        check_intervals(graph, extended_intervals);
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
    for (const auto& kv : span) {
        int lo = kv.second.first, hi = kv.second.second;
        if (hi <= lo) continue;  // contig used as a single contiguous block: no interior to splice
        string c_sample = graph->get_sample_name(kv.first);
        if (c_sample != target_sample) continue;  // only judge the target's own contigs
        for (int j = lo + 1; j < hi; ++j) {
            if (idx_path[j] == kv.first) continue;            // another piece of the same contig
            if (graph->get_sample_name(idx_path[j]) == c_sample) {
                detail = "contig " + graph->get_path_name(kv.first)
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
    for (auto& kv : used) {
        path_handle_t C = kv.first;
        if (graph->get_sample_name(C) != target_sample) continue;     // only the target's own contigs
        auto& r = kv.second;
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
        if (nonN < min_replaced) continue;   // too little real sequence replaced to judge reliably
        string added = intervals_to_sequence(graph, foreign);
        double rec = kmer_recovery(removed, added);
        ostringstream ss; ss << fixed << setprecision(1) << rec;
        cout << "#Interior graft: contig " << graph->get_path_name(C) << " replaced=" << removed.size()
             << "bp grafted=" << added.size() << "bp kmer_recovery=" << ss.str() << "%" << endl;
        if (rec < min_recovery) {
            detail = "contig " + graph->get_path_name(C) + " interior (" + to_string(removed.size())
                   + "bp) was replaced by a foreign graft sharing only " + ss.str() + "% of its k-mers";
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
    const int64_t OUTER = 20000;   // telomere window; also the minimum end-trim to consider the cap "dropped"
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
    for (auto& kv : used) {
        path_handle_t C = kv.first;
        int64_t lo = kv.second.first, hi = kv.second.second, len = length[C];
        if (lo >= OUTER) {                       // 5' end trimmed -> does the natural 5' carry a telomere?
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
        if (len - hi >= OUTER) {                 // 3' end trimmed -> does the natural 3' carry a telomere?
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
    int64_t walked = 0, shared = 0; step_handle_t s = start;
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
                                double min_flank, int64_t flank_window) {
    unordered_map<path_handle_t, unordered_map<step_handle_t, int64_t>> posidx;   // step -> forward pos, per contig (path is stable)
    auto pos_of = [&](path_handle_t C) -> unordered_map<step_handle_t, int64_t>& {
        auto it = posidx.find(C);
        if (it != posidx.end()) return it->second;
        auto& m = posidx[C]; int64_t p = 0;
        graph->for_each_step_in_path(C, [&](step_handle_t s) { m[s] = p; p += graph->get_length(graph->get_handle_of_step(s)); });
        return m;
    };
    bool changed = true;
    while (changed) {
        changed = false;
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
            if (rgi.second <= rgj.first) { lo = rgi.second; hi = rgj.first; }   // piece i is the lower piece
            else if (rgj.second <= rgi.first) { lo = rgj.second; hi = rgi.first; }  // piece j is the lower piece
            else continue;                                  // overlapping ranges -- unexpected, leave it
            if (hi <= lo) continue;                          // pieces adjacent: nothing of C was replaced
            string removed; int64_t p = 0;
            graph->for_each_step_in_path(C, [&](step_handle_t s) {
                handle_t h = graph->get_handle_of_step(s); int64_t l = graph->get_length(h);
                int64_t oa = max(p, lo), ob = min(p + l, hi); if (oa < ob) removed += graph->get_sequence(h).substr(oa - p, ob - oa); p += l; });
            int64_t nonN = 0; for (char c : removed) if (c != 'N' && c != 'n') ++nonN;

            // Apply whichever test is authoritative for this graft:
            //   * real (non-N) sequence was replaced -> k-mer recovery judges the CONTENT;
            //   * mostly-N gap (nothing to recapitulate) -> flank anchoring judges the LOCUS, i.e. the
            //     donor must stay homologous to C's own sequence over flank_window bp on both flanks.
            double rec = -1.0, fl = -1.0, fr = -1.0; bool bad = false;
            if (nonN >= min_replaced) {
                vector<tuple<step_handle_t, step_handle_t, bool>> foreign(intervals.begin() + i + 1, intervals.begin() + j);
                rec = kmer_recovery(removed, intervals_to_sequence(graph, foreign));
                bad = (rec < min_recovery);
            } else {
                unordered_set<nid_t> dL = path_node_set(graph, graph->get_path_handle_of_step(get<0>(intervals[i + 1])));
                fl = flank_fraction(graph, get<1>(intervals[i]), get<2>(intervals[i]), dL, flank_window);
                unordered_set<nid_t> dR = path_node_set(graph, graph->get_path_handle_of_step(get<0>(intervals[j - 1])));
                fr = flank_fraction(graph, get<0>(intervals[j]), !get<2>(intervals[j]), dR, flank_window);
                bad = (fl < min_flank || fr < min_flank);
            }
            if (!bad) continue;                              // faithful (content) or anchored (locus) -- keep

            ostringstream ss; ss << fixed << setprecision(1);
            if (rec >= 0) ss << "k-mer recovery " << rec << "%";
            else          ss << "flank anchoring " << fl << "%/" << fr << "%";
            cout << "#Interior graft excised: contig " << graph->get_path_name(C) << " interior ("
                 << removed.size() << "bp) -- " << ss.str()
                 << " (likely repeat-region misjoin) -- restored target sequence, kept other patches" << endl;
            // the patch visits piece i then piece j; merging them spans the whole region in the patch's
            // own direction -- (get<0> of i, get<1> of j) works for forward and reverse alike.
            tuple<step_handle_t, step_handle_t, bool> merged =
                make_tuple(get<0>(intervals[i]), get<1>(intervals[j]), get<2>(intervals[i]));
            intervals.erase(intervals.begin() + i, intervals.begin() + j + 1);
            intervals.insert(intervals.begin() + i, merged);
            changed = true;
        }
    }
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
                      double telo_threshold) {

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

    // we replace the patch with the input because it was too short
    bool to_revert = (double)patch_length / (double)tgt_length < threshold;
    if (to_revert) {
        cout << "#Reverting failed patch as it covers only " << ((double)patch_length / (double)tgt_length)
             << " of target" << endl;
    }

    // sanity guard: reject patches that replaced real interior sequence of a contig without an
    // assembly gap to justify it (repeat-region misjoin -- e.g. fragments spliced into satellite/
    // segdup through ambiguous anchors).  this overrides an otherwise-accepted (even telomere-valid)
    // patch and reverts to the input contigs.
    string interior_detail;
    if (splices_same_sample_interior(graph, in_intervals, sample_names[0], interior_detail)) {
        cout << "#Reverting patch: " << interior_detail << " (likely repeat-region misjoin)" << endl;
        to_revert = true;
    }

    // reject a foreign interior graft that recapitulates almost none of the target sequence it
    // replaced -- a repeat-region misjoin where the donor came from a different locus
    // (thresholds controlled by --graft-recovery / --graft-min-bp)
    string graft_detail;
    if (interior_graft_low_recovery(graph, in_intervals, sample_names[0], graft_recovery, graft_min_bp, graft_detail)) {
        cout << "#Reverting patch: " << graft_detail << " (likely repeat-region misjoin)" << endl;
        to_revert = true;
    }

    // telomere-preservation: a patch must not discard a telomere the target already had
    string telo_detail;
    if (discards_target_telomere(graph, in_intervals, sample_names[0], telo_threshold, telo_detail)) {
        cout << "#Reverting patch: " << telo_detail << " (target was already capped there)" << endl;
        to_revert = true;
    }

    if (!default_sample.length() && !to_revert) {
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
        // we replace the patch with the reference because there was no patch
        to_revert = !patch_happened;
        if (to_revert) {
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
            if (!has_gaps) {
                cout << "#No patching is required (the sequence contains no gaps)" << endl;
            } else {
                cout << "#Reverting to input assembly because no patches from other assemblies were found" << endl;
            }
        }
    }

    if (to_revert) {
        if (!default_sample.empty()) {
            vector<path_handle_t> default_paths;
            graph->for_each_path_of_sample(default_sample, [&](path_handle_t path_handle) {
                size_t hap = graph->get_haplotype(path_handle);
                if (hap == 0 || hap == PathMetadata::NO_HAPLOTYPE ||
                    hap == graph->get_haplotype(first_tgt_paths.front())) {
                    default_paths.push_back(path_handle);
                }
            });
            assert(default_paths.size() == 1); // todo: can relax but want to be 1 in all current use cases
            const path_handle_t& ref_path = default_paths[0];
            out_intervals.push_back(make_tuple(graph->path_begin(ref_path), graph->path_back(ref_path), false));
        } else {
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
                        bool verbose) {

    if (intervals.empty()) {
        if (verbose) {
            cerr << "[panpatch] Telomere validation: no intervals to validate" << endl;
        }
        return false;
    }

    // Get the full sequence
    string sequence = intervals_to_sequence(graph, intervals);
    int64_t seq_len = sequence.length();

    if (seq_len < 2000) {
        if (verbose) {
            cerr << "[panpatch] Telomere validation: sequence too short (" << seq_len << "bp)" << endl;
        }
        return false;
    }

    // Search window - check up to 50kb from each end
    int64_t tip_check_len = min((int64_t)50000, seq_len / 2);

    // Telomere detection uses the shared seq_has_telomere() helper (also used by the patcher, so
    // the two cannot disagree). Tips use boundary detection; the internal check does not.

    // Check for telomeres at the start (find boundary scanning forward)
    bool has_start_telomere = seq_has_telomere(sequence, 0, tip_check_len, true, false, threshold);

    // Check for telomeres at the end (find boundary scanning backward)
    bool has_end_telomere = seq_has_telomere(sequence, max((int64_t)0, seq_len - tip_check_len), seq_len, true, true, threshold);

    // Check for telomeres in the middle (internal telomeres - should NOT exist)
    // Don't find boundary here - we want to detect any telomeric sequence
    bool has_internal_telomere = false;
    if (seq_len > 2 * tip_check_len) {
        has_internal_telomere = seq_has_telomere(sequence, tip_check_len, seq_len - tip_check_len, false, false, threshold);
    }

    if (verbose) {
        cerr << "[panpatch] Telomere validation (threshold=" << threshold << "):" << endl;
        cerr << "[panpatch]   Sequence length: " << seq_len << "bp" << endl;
        cerr << "[panpatch]   Start telomere: " << (has_start_telomere ? "FOUND" : "NOT FOUND") << endl;
        cerr << "[panpatch]   End telomere: " << (has_end_telomere ? "FOUND" : "NOT FOUND") << endl;
        cerr << "[panpatch]   Internal telomeres: " << (has_internal_telomere ? "FOUND (BAD)" : "NOT FOUND (GOOD)") << endl;
    }

    return has_start_telomere && has_end_telomere && !has_internal_telomere;
}

void log_contig_telomeres(const PathHandleGraph* graph,
                          const vector<tuple<step_handle_t, step_handle_t, bool>>& intervals,
                          double threshold) {

    if (intervals.empty()) {
        return;
    }

    // Dummy string for lambda capture - will be set per contig
    string sequence;

    // Helper function to find where telomere region ends (direction=1) or starts (direction=-1)
    // Uses sliding window to detect where repeat density drops
    auto find_telomere_end = [&](int64_t start_pos, int direction, int64_t max_search) -> int64_t {
        const int64_t window_size = 500;
        const double min_density = 0.7;
        int64_t search_limit = direction > 0 ?
            min(start_pos + max_search, (int64_t)sequence.length()) :
            max(start_pos - max_search, (int64_t)0);

        if (direction > 0) {
            // Scan forward to find where telomere ends
            int64_t telomere_end = start_pos;
            for (int64_t win_start = start_pos; win_start + window_size < search_limit; win_start += 100) {
                int64_t repeats = 0;
                for (int64_t pos = win_start; pos < min(win_start + window_size, search_limit - 6); ++pos) {
                    if (sequence.substr(pos, 6) == "TTAGGG" || sequence.substr(pos, 6) == "CCCTAA") {
                        ++repeats;
                        pos += 5;
                    }
                }
                double density = 6.0 * (double)repeats / (double)window_size;
                if (density >= min_density) {
                    telomere_end = win_start + window_size;
                } else if (telomere_end > start_pos) {
                    break;
                }
            }
            return telomere_end > start_pos ? telomere_end : -1;
        } else {
            // Scan backward to find where telomere starts
            int64_t telomere_start = start_pos;
            for (int64_t win_end = start_pos; win_end - window_size > search_limit; win_end -= 100) {
                int64_t win_start = max(search_limit, win_end - window_size);
                int64_t repeats = 0;
                for (int64_t pos = win_start; pos < win_end - 6; ++pos) {
                    if (sequence.substr(pos, 6) == "TTAGGG" || sequence.substr(pos, 6) == "CCCTAA") {
                        ++repeats;
                        pos += 5;
                    }
                }
                double density = 6.0 * (double)repeats / (double)window_size;
                if (density >= min_density) {
                    telomere_start = win_start;
                } else if (telomere_start < start_pos) {
                    break;
                }
            }
            return telomere_start < start_pos ? telomere_start : -1;
        }
    };

    // Helper function to check for telomere density in a region
    // This version finds the actual telomere boundary instead of using fixed windows
    auto check_telomere_density = [&](int64_t start, int64_t end, bool find_boundary, bool is_right_end) -> pair<double, double> {
        int64_t fw_count = 0;
        int64_t r_count = 0;
        int64_t actual_start = start;
        int64_t actual_end = end;

        // If requested, find where telomere actually ends/starts
        if (find_boundary) {
            if (is_right_end) {
                // Scan backward from end to find where telomere starts
                int64_t telomere_start = find_telomere_end(end, -1, end - start);
                if (telomere_start >= start && telomere_start < end) {
                    actual_start = telomere_start;
                }
            } else {
                // Scan forward from start to find where telomere ends
                int64_t telomere_end = find_telomere_end(start, 1, end - start);
                if (telomere_end > start) {
                    actual_end = telomere_end;
                }
            }
        }

        for (int64_t pos = actual_start; pos < actual_end - 6; ++pos) {
            if (sequence.substr(pos, 6) == "TTAGGG") {
                ++fw_count;
                pos += 5;
            } else if (sequence.substr(pos, 6) == "CCCTAA") {
                ++r_count;
                pos += 5;
            }
        }

        int64_t region_len = actual_end - actual_start;
        if (region_len < 500) {  // Require at least 500bp of telomere
            return make_pair(0.0, 0.0);
        }

        double fw_density = 6. * ((double)fw_count / (double)region_len);
        double r_density = 6. * ((double)r_count / (double)region_len);

        return make_pair(fw_density, r_density);
    };

    // Group intervals by path to analyze each contig separately
    unordered_map<path_handle_t, vector<tuple<step_handle_t, step_handle_t, bool>>> path_intervals;
    for (const auto& interval : intervals) {
        path_handle_t path = graph->get_path_handle_of_step(get<0>(interval));
        path_intervals[path].push_back(interval);
    }

    // Analyze each contig
    for (const auto& path_int : path_intervals) {
        path_handle_t path = path_int.first;
        const auto& path_ints = path_int.second;

        sequence = intervals_to_sequence(graph, path_ints);
        int64_t seq_len = sequence.length();

        if (seq_len < 50) {
            continue;  // Too short to analyze
        }

        // Check tip regions - use up to 50kb search window
        int64_t tip_check_len = min((int64_t)50000, seq_len / 2);

        auto left_densities = check_telomere_density(0, tip_check_len, true, false);
        auto right_densities = check_telomere_density(max((int64_t)0, seq_len - tip_check_len), seq_len, true, true);

        double left_max_density = max(left_densities.first, left_densities.second);
        double right_max_density = max(right_densities.first, right_densities.second);

        bool has_left = left_max_density >= threshold;
        bool has_right = right_max_density >= threshold;

        // Log all contigs with telomere information
        cout << "#Contig " << graph->get_path_name(path)
             << " len=" << seq_len << "bp"
             << " left=" << (has_left ? "YES" : "NO")
             << "(" << fixed << setprecision(3) << left_max_density << ")"
             << " right=" << (has_right ? "YES" : "NO")
             << "(" << fixed << setprecision(3) << right_max_density << ")"
             << endl;
    }
}

