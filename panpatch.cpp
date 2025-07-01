#include <unordered_set>
#include <cassert>
#include <map>
#include "panpatch.hpp"
#include <limits>

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
                                                                   const vector<path_handle_t>& other_paths) {

    assert(ref_anchors.size() > 1);
    assert(tgt_paths.size() > 0);
    
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
                    if (!tgt_ranks.count(graph->get_path_handle_of_step(next_anchor.first))) {
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
                pos_2 -= graph->get_length(graph->get_handle_of_step(get<0>(interval)));
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
                                                               const unordered_map<string, vector<path_handle_t>>& sample_covers) {

       
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
                                                                                           other_paths);

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

    return extended_intervals;
}

bool revert_bad_patch(const PathHandleGraph* graph,
                      const path_handle_t& ref_path,
                      const vector<path_handle_t>& tgt_paths,
                      const vector<string>& sample_names,
                      const vector<tuple<step_handle_t, step_handle_t, bool>>& in_intervals,
                      vector<tuple<step_handle_t, step_handle_t, bool>>& out_intervals,
                      string default_sample,
                      double threshold) {

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

    if (!default_sample.length() && !to_revert) {
        bool patch_happened = false;
        for (const auto& interval : in_intervals) {
            if (graph->get_sample_name(graph->get_path_handle_of_step(get<0>(interval))) != sample_names[0]) {
                patch_happened = true;
                break;
            }
        }
        // we replace the patch with the reference because there was no patch
        to_revert = !patch_happened;
        if (to_revert) {
            cout << "#Reverting to input assembly because no patches from other assemblies were found" << endl;            
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
    for (const auto& interval : intervals) {
        if (get<2>(interval) == false) {
            for (step_handle_t step = get<0>(interval); step != get<1>(interval); step = graph->get_next_step(step)) {
                assert(step != graph->path_end(graph->get_path_handle_of_step(step)) &&
                       step != graph->path_front_end(graph->get_path_handle_of_step(step)));
            }
        } else {
            for (step_handle_t step = get<1>(interval); step != get<0>(interval); step = graph->get_next_step(step)) {
                assert(step != graph->path_end(graph->get_path_handle_of_step(step)) &&
                       step != graph->path_front_end(graph->get_path_handle_of_step(step)));
            }            
        }
    }    
}

