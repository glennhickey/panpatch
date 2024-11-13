#include <unordered_set>
#include <cassert>
#include <map>
#include "panpatch.hpp"

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

void greedy_patch(const PathHandleGraph* graph,
                  const path_handle_t& ref_path,
                  const vector<path_handle_t>& tgt_paths,                  
                  const vector<string>& sample_names,
                  const unordered_map<string, vector<path_handle_t>>& sample_covers) {

//#ifdef debug
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
//#endif
    
}
                  
