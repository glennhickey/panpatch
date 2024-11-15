#include <unordered_set>
#include <cassert>
#include <map>
#include "panpatch.hpp"
#include <limits>

//#define debug

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
#ifdef debug
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
#ifdef debug
                    cerr << "room " << room << endl;
                    cerr << "next offset = " << w << " - " << cur_win_length << " = " << next_offset << endl;
#endif
                    if (step == next_step) {
                        next_offset += start_offset;
                    }
                    cur_win_length = w;
                    assert(next_offset <= next_len);
#ifdef debug
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
        sample_hap_coverage[sample][haplotype] += path_cov.second;
        sample_hap_count[sample][haplotype] += 1;
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
                                              const unordered_set<path_handle_t>& path_set) {
    unordered_map<int64_t, int64_t> anchors;
    int64_t pos = 0;
    
    pair<int64_t, int64_t> prev_anchor;
    pair<int64_t, int64_t> prev_prev_anchor;
    unordered_set<path_handle_t> prev_anchor_paths;
    unordered_set<path_handle_t> prev_prev_anchor_paths;

    // todo: we can optimise this a lot by targeting edges of tgt paths...
    graph->for_each_step_in_path(ref_path, [&](step_handle_t step) {
        handle_t handle = graph->get_handle_of_step(step);
        unordered_set<path_handle_t> covered_paths;
        graph->for_each_step_on_handle(handle, [&](step_handle_t other_step) {
            path_handle_t other_path = graph->get_path_handle_of_step(other_step);
            if (path_set.count(other_path) && other_path != ref_path) {
                covered_paths.insert(graph->get_path_handle_of_step(other_step));
            }
        });
        // we currently have no use for ref-only anchors
        if (!covered_paths.empty()) {            
            if (prev_prev_anchor_paths.empty()) {
                prev_prev_anchor_paths = covered_paths;
                prev_prev_anchor = make_pair(graph->get_id(handle), pos);
            } else if (prev_anchor_paths.empty()) {
                prev_anchor_paths = covered_paths;
                prev_anchor = make_pair(graph->get_id(handle), pos);
            } else {
                bool p_same = covered_paths == prev_anchor_paths;
                bool pp_same = prev_prev_anchor_paths == prev_anchor_paths;
                if (p_same && pp_same) {
                    // we have three anchors the same, don't do anything but slide prev_anchor to here
                    prev_anchor = make_pair(graph->get_id(handle), pos);
                } else {
                    // save the first anchor and update the other two
                    anchors.insert(prev_prev_anchor);
                    prev_prev_anchor = prev_anchor;
                    prev_prev_anchor_paths = prev_anchor_paths;
                    prev_anchor = make_pair(graph->get_id(handle), pos);
                    prev_anchor_paths = covered_paths;
                }
            }
        }
        pos += graph->get_length(handle);
    });

    if (!prev_prev_anchor_paths.empty()) {
        anchors.insert(prev_prev_anchor);
    }
    if (!prev_anchor_paths.empty()) {
        anchors.insert(prev_anchor);
    }

    return anchors;
}

pair<step_handle_t, bool> find_next_anchor_on_path(const PathHandleGraph* graph,
                                                   const unordered_map<int64_t, int64_t>& anchors,
                                                   step_handle_t step,
                                                   int64_t pos) {

    path_handle_t path = graph->get_path_handle_of_step(step);
    
    // search forward
    step_handle_t next_step = graph->get_next_step(step);
    for (; next_step != graph->path_end(path); next_step = graph->get_next_step(next_step)) {
        handle_t next_handle = graph->get_handle_of_step(next_step);
        if (anchors.count(graph->get_id(next_handle)) &&
            anchors.at(graph->get_id(next_handle)) > pos) {
            return make_pair(next_step, false);
        }
    }

    // search backward
    next_step = graph->get_previous_step(step);
    for (; next_step != graph->path_front_end(path); next_step = graph->get_previous_step(next_step)) {
        handle_t next_handle = graph->get_handle_of_step(next_step);
        if (anchors.count(graph->get_id(next_handle)) &&
            anchors.at(graph->get_id(next_handle)) > pos) {
            return make_pair(next_step, true);
        }
    }

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
                cur_pos = handle_pos;
                cur_handle = handle;
            }
        }
    });
        
    while (true) {
        // find the steps on the handle and sort them using the path priority
        // [paths that don't have a priority are ignored which is important]
        // todo: could be in weird regions we need to assess contiguiuty with
        // previous anchor, but not trying for now. 
        multimap<int64_t, step_handle_t> sorted_steps;
        graph->for_each_step_on_handle(cur_handle, [&](step_handle_t step) {
            path_handle_t path = graph->get_path_handle_of_step(step);
            if (path_rank.count(path)) {
                sorted_steps.insert(make_pair(path_rank.at(path), step));
            }
        });

        bool found_next = false;
        for (const auto& rank_step : sorted_steps) {
            const step_handle_t& step = rank_step.second;

            pair<step_handle_t, bool> next_anchor = find_next_anchor_on_path(graph, ref_anchors, step, cur_pos);
            if (next_anchor.first != graph->path_end(graph->get_path_handle_of_step(step))) {
                interval_cover.push_back(make_tuple(step, next_anchor.first, next_anchor.second));
                found_next = true;
                // slide position to the next anchor
                cur_handle = graph->get_handle_of_step(next_anchor.first);
                cur_pos = ref_anchors.at(graph->get_id(cur_handle));
                break;
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

    return smoothed_intervals;    
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

    for (const auto& interval : intervals) {
        path_handle_t path = graph->get_path_handle_of_step(get<0>(interval));
        int64_t pos_1 = path_index[path][get<0>(interval)];
        int64_t pos_2 = path_index[path][get<1>(interval)];
        pos_2 += graph->get_length(graph->get_handle_of_step(get<1>(interval)));
        cout << graph->get_path_name(path) << "\t" << pos_1 << "\t" << pos_2 << endl;
    }
}


void greedy_patch(const PathHandleGraph* graph,
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

/*
    // start by snapping the target to the reference
    multimap<pair<int64_t, int64_t>, path_handle_t> snapped_target = sort_overlapping_paths(graph, ref_path, tgt_paths);

#ifdef debug
    cerr << "snapped target" << endl;
    for (const auto& region_path : snapped_target) {
        cerr << region_path.first.first << "-" << region_path.first.second << ": "
             << graph->get_path_name(region_path.second) << endl;
    }
    cerr << endl;
#endif
*/
    // find the anchors along the reference path using the relevant paths
    unordered_set<path_handle_t> relevant_paths = {ref_path};
    for (const auto& sample_paths : sample_covers) {
        for (const path_handle_t& path : sample_paths.second) {
            relevant_paths.insert(path);
        }
    }
    unordered_map<int64_t, int64_t> ref_anchors = find_anchors(graph, ref_path, relevant_paths);

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
        return;
    }

#ifdef debug
    cerr << "number of patched intervals found " << patched_intervals.size() << endl;
#endif
    vector<tuple<step_handle_t, step_handle_t, bool>> smoothed_intervals = smooth_intervals(graph, patched_intervals);

#ifdef debug
    cerr << "number of smoothed intervals found " << smoothed_intervals.size() << endl;
#endif
    cout << "Patched assembly on " << graph->get_locus_name(ref_path) << " for "
         << graph->get_sample_name(tgt_paths.front()) << "#"
         << graph->get_haplotype(tgt_paths.front()) << ":" << endl;
    print_intervals(graph, smoothed_intervals);
    cout << endl;
    
}
                  
