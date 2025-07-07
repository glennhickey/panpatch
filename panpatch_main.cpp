#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <fstream>
#include <memory>
#include <unistd.h>
#include <getopt.h>
#include <omp.h>

#include "handlegraph/path_handle_graph.hpp"
#include "bdsg/packed_graph.hpp"
#include "bdsg/hash_graph.hpp"
#include "bdsg/snarl_distance_index.hpp"
#include "bdsg/overlays/overlay_helper.hpp"
#include "panpatch.hpp"

//#define debug

using namespace std;
using namespace handlegraph;
using namespace bdsg;


// from hal2vg/clip-vg.cpp
static unique_ptr<PathHandleGraph> load_graph(istream& graph_stream);

// if the patched contig isn't at least this much as big as the set of input
// contigs to patch, then consider the patch failed
static const double fail_threshold = 0.95;

static const size_t fasta_width = 80;

void help(char** argv) {
  cerr << "usage: " << argv[0] << " [options] <graph> " << endl
       << "Use a pangenome alignment (of single sample and reference) to make patched assembly" << endl
       << endl
       << "options: " << endl
       << "    -p, --progress               Print progress" << endl
       << "    -r, --reference STRING       Reference sample" << endl
       << "    -s, --sample STRING          Input sample. Multiple allowed. Order specifies priority" << endl
       << "    -f, --fasta FILE             Output the patched assembly to this fasta file" << endl
       << "    -w, --window SIZE            Size of window used for computing identity for haplotype matching [1000]" << endl
       << "    -e, --default-sample STRING  If unable to patch, use contig from this sample (if diploid, haplotypes must be consistent with first sample!)" << endl
       << "    -t, --threads N              Number of threads to use [default: all available]" << endl      
       << endl;
}    

int main(int argc, char** argv) {

    string ref_sample;
    vector<string> sample_names;
    bool progress = false;
    string out_fasta_filename;
    string default_sample;
    int c;
    int64_t window_size = 1000;
    bool ref_default = false;
    optind = 1; 
    while (true) {

        static const struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"progress", no_argument, 0, 'p'},
            {"reference", required_argument, 0, 'r'},
            {"sample", required_argument, 0, 's'},
            {"fasta", required_argument, 0, 'f'},
            {"window", required_argument, 0, 'w'},
            {"default-sample", required_argument, 0, 'e'},
            {"threads", required_argument, 0, 't'},            
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hpr:s:f:w:e:t:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'p':
            progress = true;
            break;
        case 'r':
            ref_sample = optarg;
            break;
        case 's':
            sample_names.push_back(optarg);
            break;
        case 'f':
            out_fasta_filename = optarg;
            break;
        case 'w':
            window_size = atoi(optarg);
            break;
        case 'e':
            default_sample = optarg;
            break;
        case 't':
        {
            int num_threads = stoi(optarg);
            if (num_threads <= 0) {
                cerr << "[vg2maf] error: Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                exit(1);
            }
            omp_set_num_threads(num_threads);
            break;
        }                                    
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 1) {
        help(argv);
        return 1;
    }
    if (sample_names.empty()) {
        cerr << "[panpatch] error: -s must be used to specify at least one sample name to prioritize" << endl;
        return 1;
    }
    if (ref_sample.empty()) {
        cerr << "[panpatch] error: -r must be used to specify a reference sample" << endl;
        return 1;
    }    
    string graph_filename = argv[optind++];
    ifstream graph_stream(graph_filename);
    if (!graph_stream) {
        cerr << "[panpatch] error: Unable to open input graph " << graph_filename << endl;
        return 1;
    }
    ofstream out_fasta_file;
    if (!out_fasta_filename.empty()) {
        out_fasta_file.open(out_fasta_filename);
        if (!out_fasta_file) {
            cerr << "[panpatch] error: Unable to open fasta file for writing: " << out_fasta_filename << endl;
            return 1;
        }
    }
    
    if (progress) {
        cerr << "[panpatch]: Using " << get_thread_count() << (get_thread_count() > 1 ? " threads" : " thread") << endl;
    }
    unique_ptr<PathHandleGraph> base_graph = load_graph(graph_stream);
    graph_stream.close();
    if (progress) {
        cerr << "[panpatch]: Loaded graph" << endl;
    }
    bdsg::ReferencePathOverlayHelper overlay_helper;
    PathPositionHandleGraph* graph = overlay_helper.apply(base_graph.get());
    if (progress && dynamic_cast<PathPositionHandleGraph*>(base_graph.get()) == nullptr) {
        cerr << "[panpatch]: Applied position overlay" << endl;
    }

    // pull out the (one and only) reference path
    vector<path_handle_t> ref_paths;
    graph->for_each_path_of_sample(ref_sample, [&](path_handle_t ref_path) {
        ref_paths.push_back(ref_path);
    });
    if (ref_paths.size() != 1) {
        cerr << "[panpatch]: Exactly 1 path for reference sample " << ref_sample << " expected. " << ref_paths.size()
             << " found. panpatch only works when there's 1 for now..." << endl;
        return 1;        
    }
    path_handle_t ref_path = ref_paths.front();
    if (progress) {
        cerr << "[panpatch]: Selected reference path " << graph->get_path_name(ref_path) << endl;
    }

    // pull out all other paths selected by -s
    vector<path_handle_t> other_paths;
    for (const string& sample : sample_names) {
        graph->for_each_path_of_sample(sample, [&](path_handle_t path_handle) {
            other_paths.push_back(path_handle);
        });
    }
    if (progress) {
        cerr << "[panpatch]: Selected " << other_paths.size() << " (non-reference) paths" << endl;
    }

    // pull out the target paths (the ones we want to patch) and sort them by haplotype
    map<int64_t, vector<path_handle_t>> target_paths;
    graph->for_each_path_of_sample(sample_names.front(), [&](path_handle_t path_handle) {
        target_paths[graph->get_haplotype(path_handle)].push_back(path_handle);
    });
    if (progress) {
        cerr << "[panpatch]: Target sample " << sample_names.front() << " has " << target_paths.size() << " paths" << endl;
    }
    if (target_paths.empty()) {
        cerr << "[panpatch]: Error: No paths found for target sample " << sample_names.front()
             << " in graph containing reference path " << graph->get_path_name(ref_path) << endl;
        return 1;
    }

    // we patch each target haplotype independently, greedily selecting other haplotypes
    // up front using this simple coverage calculation
    for (const auto& hap_tgts : target_paths) {
        // break out the best-covering haplotype of each other sample
        unordered_map<path_handle_t, double> coverage_map = compute_overlap_identity(graph, hap_tgts.second, other_paths, window_size);
        if (progress) {
            cerr << "[panpatch]: Computed coverage for hap " << hap_tgts.first << " paths:";
            for (const auto& tgt_path : hap_tgts.second) {
                cerr << " " << graph->get_path_name(tgt_path);
            }
            cerr << ":" << endl;

            for (const auto& cov : coverage_map) {
                if (graph->get_sample_name(cov.first) != graph->get_sample_name(target_paths.begin()->second.front())) {
                    cerr << "[panpatch]:    " << graph->get_path_name(cov.first) << " " << cov.second << endl;
                }
            }
        }

        unordered_map<string, vector<path_handle_t>> sample_covers = select_sample_covers(graph, coverage_map);

        if (progress) {
            cerr << "[panpatch]: Sample cover selection:" << endl;
            for (const auto& sc : sample_covers) {
                cerr << "[panpatch]:    " << sc.first;
                for (const auto& p : sc.second) {
                    cerr << " " << graph->get_path_name(p);
                }
                cerr << endl;
            }
        }
                
        // run the patching on the given target haplotype, using seleted haplotypes of the
        // other paths
        if (progress) {
            cerr << "[panpatch]: Running greedy patch selection" << endl;
        }
        vector<tuple<step_handle_t, step_handle_t, bool>> patched_intervals = greedy_patch(
            graph, ref_path, hap_tgts.second, sample_names, sample_covers);

        vector<tuple<step_handle_t, step_handle_t, bool>> input_intervals;
        bool reverted = revert_bad_patch(graph, ref_path, hap_tgts.second, sample_names,
                                         patched_intervals, input_intervals, 
                                         default_sample, fail_threshold);
        if (reverted) {
            patched_intervals = input_intervals;
        }
        

        // print the intervals to cout
        cout << "#Patched assembly on " << graph->get_locus_name(ref_path) << " for "
             << graph->get_sample_name(hap_tgts.second.front()) << "#"
             << graph->get_haplotype(hap_tgts.second.front()) << ":" << endl;
        print_intervals(graph, patched_intervals);
        cout << endl;

        // save the intervals to the fasta
        if (!out_fasta_filename.empty()) {
            if (reverted) {
                // note: we have two modes, either we've reverted to the original contigs and
                // we just write them out one by one.  Or we made a single t2t patch.
                // todo: what isn't supported (and probably should be!!!) is partial patching
                // where we patch a few contigs but output is not single t2t contig.
                if (progress) {
                    cerr << "[panpatch]: Writing input contig(s) to FASTA" << endl;
                }
                for (const auto& interval : patched_intervals) {
                    path_handle_t interval_path = graph->get_path_handle_of_step(get<0>(interval));
                    string contig_name = graph->get_path_name(interval_path);
                    string sequence = intervals_to_sequence(graph, {interval});
                    out_fasta_file << ">" << contig_name << endl;
                    for (size_t written = 0; written < sequence.length(); written += fasta_width) {
                        out_fasta_file << sequence.substr(written, min(fasta_width, sequence.length() - written)) << "\n";
                    }
                }
            } else {
                if (progress) {
                    cerr << "[panpatch]: Writing patched contig to FASTA" << endl;
                }         
                string contig_name = graph->get_locus_name(ref_path) + "_hap_" + std::to_string(hap_tgts.first);
                string sequence = intervals_to_sequence(graph, patched_intervals);
                out_fasta_file << ">" << contig_name << endl;
                for (size_t written = 0; written < sequence.length(); written += fasta_width) {
                    out_fasta_file << sequence.substr(written, min(fasta_width, sequence.length() - written)) << "\n";
                }
            }
        }
    }
    return 0;
}

unique_ptr<PathHandleGraph> load_graph(istream& graph_stream) {

    char magic_bytes[4];
    graph_stream.read(magic_bytes, 4);
    uint32_t magic_number = ntohl(*((uint32_t*) magic_bytes));
    graph_stream.clear();
    graph_stream.seekg(0, ios::beg);

    PathHandleGraph* graph;
    if (magic_number == PackedGraph().get_magic_number()) {
        graph = new PackedGraph();
    } else if (magic_number == HashGraph().get_magic_number()) {
        graph = new HashGraph();
    }  else {
        cerr << "Unable to parse input graph with magic number " << magic_number << endl;
        exit(1);
    }
    dynamic_cast<SerializableHandleGraph*>(graph)->deserialize(graph_stream);

    return unique_ptr<PathHandleGraph>(graph);
}
