#include <unordered_set>
#include <cassert>
#include "panpatch.hpp"

#define debug

unordered_map<path_handle_t, double> compute_overlap_identity(const PathHandleGraph* graph,
                                                              const path_handle_t& ref_path,
                                                              int64_t w) {
    unordered_map<path_handle_t, int64_t> total_window_overlaps;
    unordered_map<path_handle_t, int64_t> total_window_counts;
    unordered_map<path_handle_t, int64_t> window_overlaps;

    int64_t start_offset = 0;
    int64_t cur_win_length = 0;
    for (step_handle_t step = graph->path_begin(ref_path);
         step != graph->path_end(ref_path);) {
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
             next_step != graph->path_end(ref_path);
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
             cur_step != graph->path_end(ref_path);
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
                if (other_path != ref_path) {
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

    unordered_map<path_handle_t, double> coverage_map;
    for (const auto& path_overlaps : total_window_overlaps) {
        coverage_map[path_overlaps.first] = (double)path_overlaps.second / (double)total_window_counts[path_overlaps.first];
    }

    return coverage_map;
}
