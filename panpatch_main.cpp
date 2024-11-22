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


void help(char** argv) {
  cerr << "usage: " << argv[0] << " [options] <graph> " << endl
       << "Use a pangenome alignment (of single sample and reference) to make patched assembly" << endl
       << endl
       << "options: " << endl
       << "    -p, --progress           Print progress" << endl
       << "    -r, --reference STRING   Reference sample" << endl
       << "    -s, --sample STRING      Input sample. Multiple allowed. Order specifies priority" << endl
       << "    -f, --fasta FILE         Output the patched assembly to this fasta file" << endl
       << "    -t, --threads N          Number of threads to use [default: all available]" << endl      
       << endl;
}    

int main(int argc, char** argv) {

    string ref_sample;
    vector<string> sample_names;
    bool progress = false;
    string out_fasta_filename;
    int c;
    optind = 1; 
    while (true) {

        static const struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"progress", no_argument, 0, 'p'},
            {"reference", required_argument, 0, 'r'},
            {"sample", required_argument, 0, 's'},
            {"fasta", required_argument, 0, 'f'},
            {"threads", required_argument, 0, 't'},            
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hpr:s:f:t:",
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

    // pull out all other paths selected by -s
    vector<path_handle_t> other_paths;
    for (const string& sample : sample_names) {
        graph->for_each_path_of_sample(sample, [&](path_handle_t path_handle) {
            other_paths.push_back(path_handle);
        });
    }

    // pull out the target paths (the ones we want to patch) and sort them by haplotype
    unordered_map<int64_t, vector<path_handle_t>> target_paths;
    graph->for_each_path_of_sample(sample_names.front(), [&](path_handle_t path_handle) {
        target_paths[graph->get_haplotype(path_handle)].push_back(path_handle);
    });

    // we patch each target haplotype independently, greedily selecting other haplotypes
    // up front using this simple coverage calculation
    for (const auto& hap_tgts : target_paths) {
        // break out the best-covering haplotype of each other sample
        unordered_map<path_handle_t, double> coverage_map = compute_overlap_identity(graph, hap_tgts.second, other_paths);
        unordered_map<string, vector<path_handle_t>> sample_covers = select_sample_covers(graph, coverage_map);

        // run the patching on the given target haplotype, using seleted haplotypes of the
        // other paths
        vector<tuple<step_handle_t, step_handle_t, bool>> patched_intervals = greedy_patch(
            graph, ref_path, hap_tgts.second, sample_names, sample_covers);

        // print the intervals to cout
        cout << "Patched assembly on " << graph->get_locus_name(ref_path) << " for "
             << graph->get_sample_name(hap_tgts.second.front()) << "#"
             << graph->get_haplotype(hap_tgts.second.front()) << ":" << endl;
        print_intervals(graph, patched_intervals);
        cout << endl;

        // save the intervals to the fasta
        if (!out_fasta_filename.empty()) {
            string contig_name = graph->get_locus_name(ref_path) + "_hap_" + std::to_string(hap_tgts.first);
            string sequence = intervals_to_sequence(graph, patched_intervals);
            out_fasta_file << ">" << contig_name << endl;
            static const size_t width = 80;
            for (size_t written = 0; written < sequence.length(); written += width) {
                out_fasta_file << sequence.substr(written, min(width, sequence.length() - written)) << "\n";
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
