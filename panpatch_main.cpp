#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <fstream>
#include <memory>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <sstream>
#include <iomanip>
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

// the per-patch report table (TSV) printed to stdout
static const char* PATCH_TABLE_HEADER =
    "chrom\thap\ttype\ttarget\ttarget_bp\tdonor\tdonor_bp\treplaced_bp\tkmer%\tflankL%\tflankR%\tdecision\treason";
static string fmt_pct(double v) { if (v < 0) return "."; ostringstream s; s << fixed << setprecision(1) << v; return s.str(); }
static string fmt_bp(int64_t v) { return v < 0 ? string(".") : to_string(v); }
// per-haplotype FASTA filename: insert ".hap<N>" before the extension of the -f base
static string fasta_name(const string& base, int64_t hap) {
    string tag = ".hap" + to_string(hap);
    size_t slash = base.find_last_of('/');
    size_t dot = base.find_last_of('.');
    return (dot == string::npos || (slash != string::npos && dot < slash))
           ? base + tag : base.substr(0, dot) + tag + base.substr(dot);
}
static void print_patch_row(ostream& o, const PatchRecord& pr) {
    o << pr.chrom << '\t' << pr.hap << '\t' << pr.type << '\t'
      << pr.target << '\t' << pr.target_bp << '\t'
      << pr.donor << '\t' << pr.donor_bp << '\t'
      << fmt_bp(pr.replaced_bp) << '\t'
      << fmt_pct(pr.kmer) << '\t' << fmt_pct(pr.flankL) << '\t' << fmt_pct(pr.flankR) << '\t'
      << (pr.accepted ? "accepted" : "rejected") << '\t'
      << (pr.reason.empty() ? "." : pr.reason) << '\n';
}

static const size_t fasta_width = 80;

void help(char** argv) {
  cerr << "usage: " << argv[0] << " [options] <graph> [graph2 ...]" << endl
       << "Use a pangenome alignment (of single sample and reference) to make patched assembly" << endl
       << endl
       << "options: " << endl
       << "    -p, --progress               Print progress" << endl
       << "    -r, --reference STRING       Reference sample" << endl
       << "    -s, --sample STRING          Input sample. Multiple allowed. Order specifies priority" << endl
       << "    -f, --fasta FILE             Output the patched assembly to FASTA, one file per haplotype" << endl
       << "                                 (FILE.hap1.fa, FILE.hap2.fa, ...); written only on full success" << endl
       << "        --bed FILE               Write the patched-assembly intervals (BED) to FILE (default: none);" << endl
       << "                                 written only on full success. The report still goes to stdout" << endl
       << "    -w, --window SIZE            Size of window used for computing identity for haplotype matching [1000]" << endl
       << "    -e, --default-sample STRING  If unable to patch, use contig from this sample (if diploid, haplotypes must be consistent with first sample!)" << endl
       << "    -t, --threads N              Number of threads to use [default: all available]" << endl
       << "    -T, --require-telomeres      Require telomeres at both ends (no internal): patch a missing terminal telomere from another assembly when possible, else revert" << endl
       << "    -M, --max-telomere-patch N   Max bp of target sequence a -T telomere patch may replace at a contig end [500000]" << endl
       << "    -b, --exclude-bed FILE       BED file of target regions to exclude from patching" << endl
       << "        --min-cover FLOAT        Revert a patch covering less than this fraction of the input length [0.95]" << endl
       << "        --telomere-threshold F   Min telomere hexamer density to call a telomere (with -T) [0.8]" << endl
       << "        --graft-recovery FLOAT   Revert a foreign interior graft sharing less than this % of the replaced k-mers [50]" << endl
       << "        --graft-min-bp N         Apply --graft-recovery only when at least this many non-N bp are replaced [10000]" << endl
       << "        --min-flank FLOAT        Revert a foreign interior graft anchored to less than this % of the target flank [50]" << endl
       << "        --flank-window N         Window (bp) each side of a graft over which flank anchoring is measured [500000]" << endl
       << endl;
}    

int main(int argc, char** argv) {

    string ref_sample;
    vector<string> sample_names;
    bool progress = false;
    string out_fasta_filename;
    string out_bed_filename;
    string default_sample;
    string exclude_bed_filename;
    int c;
    int64_t window_size = 1000;
    bool require_telomeres = false;
    int64_t max_telomere_patch = 500000;
    double fail_threshold = 0.95;
    double telo_threshold = 0.8;
    double graft_recovery = 50.0;
    int64_t graft_min_bp = 10000;
    double min_flank = 50.0;
    int64_t flank_window = 500000;
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
            {"require-telomeres", no_argument, 0, 'T'},
            {"max-telomere-patch", required_argument, 0, 'M'},
            {"exclude-bed", required_argument, 0, 'b'},
            {"min-cover", required_argument, 0, 1001},
            {"telomere-threshold", required_argument, 0, 1002},
            {"graft-recovery", required_argument, 0, 1003},
            {"graft-min-bp", required_argument, 0, 1004},
            {"min-flank", required_argument, 0, 1005},
            {"flank-window", required_argument, 0, 1006},
            {"bed", required_argument, 0, 1007},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hpr:s:f:w:e:t:TM:b:",
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
                cerr << "[panpatch] error: Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                exit(1);
            }
            omp_set_num_threads(num_threads);
            break;
        }
        case 'T':
            require_telomeres = true;
            break;
        case 'M':
        {
            char* m_end = nullptr;
            max_telomere_patch = strtol(optarg, &m_end, 10);
            if (m_end == optarg || *m_end != '\0' || max_telomere_patch < 0) {
                cerr << "[panpatch] error: --max-telomere-patch (-M) must be a non-negative integer" << endl;
                return 1;
            }
            break;
        }
        case 'b':
            exclude_bed_filename = optarg;
            break;
        case 1001:
            fail_threshold = atof(optarg);
            break;
        case 1002:
            telo_threshold = atof(optarg);
            break;
        case 1003:
            graft_recovery = atof(optarg);
            break;
        case 1004:
            graft_min_bp = strtoll(optarg, nullptr, 10);
            break;
        case 1005:
            min_flank = atof(optarg);
            break;
        case 1006:
            flank_window = strtoll(optarg, nullptr, 10);
            break;
        case 1007:
            out_bed_filename = optarg;
            break;
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
    // validate numeric thresholds (a nonsensical value silently disables or inverts a guard)
    if (fail_threshold < 0 || fail_threshold > 1) { cerr << "[panpatch] error: --min-cover must be in [0,1]" << endl; return 1; }
    if (telo_threshold < 0 || telo_threshold > 1) { cerr << "[panpatch] error: --telomere-threshold must be in [0,1]" << endl; return 1; }
    if (graft_recovery < 0 || graft_recovery > 100) { cerr << "[panpatch] error: --graft-recovery must be in [0,100]" << endl; return 1; }
    if (min_flank < 0 || min_flank > 100) { cerr << "[panpatch] error: --min-flank must be in [0,100]" << endl; return 1; }
    if (graft_min_bp < 0) { cerr << "[panpatch] error: --graft-min-bp must be >= 0" << endl; return 1; }
    if (flank_window <= 0) { cerr << "[panpatch] error: --flank-window must be > 0" << endl; return 1; }
    if (max_telomere_patch < 0) { cerr << "[panpatch] error: --max-telomere-patch must be >= 0" << endl; return 1; }
    // one or more input graphs, processed in lexicographic order (report/BED/FASTA are concatenated)
    vector<string> graph_filenames;
    while (optind < argc) graph_filenames.push_back(argv[optind++]);
    if (graph_filenames.empty()) {
        cerr << "[panpatch] error: at least one input graph is required" << endl;
        return 1;
    }
    sort(graph_filenames.begin(), graph_filenames.end());
    graph_filenames.erase(unique(graph_filenames.begin(), graph_filenames.end()), graph_filenames.end());

    BedRegions bed_regions;
    if (!exclude_bed_filename.empty()) {
        bed_regions = parse_bed_file(exclude_bed_filename);
        if (progress) {
            int64_t total_regions = 0;
            for (const auto& br : bed_regions) total_regions += br.second.size();
            cerr << "[panpatch]: Loaded " << total_regions << " regions from " << bed_regions.size()
                 << " contigs in BED exclusion file" << endl;
        }
    }

    if (progress) {
        cerr << "[panpatch]: Using " << get_thread_count() << (get_thread_count() > 1 ? " threads" : " thread") << endl;
    }

    // pre-scan all inputs: a misspelled -r/-s sample fails here, before any patching (so no partial output)
    {
        set<string> present_samples;
        for (const string& gf : graph_filenames) {
            ifstream gs(gf);
            if (!gs) { cerr << "[panpatch] error: Unable to open input graph " << gf << endl; return 1; }
            unique_ptr<PathHandleGraph> g = load_graph(gs);
            g->for_each_path_handle([&](path_handle_t p) { present_samples.insert(g->get_sample_name(p)); });
        }
        vector<string> missing;
        if (!present_samples.count(ref_sample)) missing.push_back(ref_sample);
        for (const string& s : sample_names) if (!present_samples.count(s)) missing.push_back(s);
        if (!missing.empty()) {
            cerr << "[panpatch] error: sample(s) not found in any input graph:";
            for (const string& m : missing) cerr << " " << m;
            cerr << endl;
            return 1;
        }
    }

    // haplotype -> open temp FASTA stream (<FILE>.hap<N>.fa.tmp), streamed during the run and renamed to
    // the final name only after every graph succeeds, so a failure never leaves a partial FASTA and the
    // whole genome is not held in RAM
    map<int64_t, ofstream> fasta_tmp;
    // accumulated BED intervals (small); written to a temp file and renamed at the end alongside the FASTAs
    ostringstream bed_ss;
    set<string> seen_loci;   // to warn on duplicate reference loci across graphs (would collide in the output)

    // the patch report streams to stdout: a TSV table, one row per candidate patch, per contig
    cout << PATCH_TABLE_HEADER << "\n";

    for (const string& graph_filename : graph_filenames) {
        ifstream graph_stream(graph_filename);
        if (!graph_stream) {
            cerr << "[panpatch] error: Unable to open input graph " << graph_filename << endl;
            return 1;
        }
        if (progress) {
            cerr << "[panpatch]: Processing " << graph_filename << endl;
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
        cerr << "[panpatch]: skipping " << graph_filename << ": expected exactly 1 reference path for "
             << ref_sample << ", found " << ref_paths.size() << endl;
        continue;
    }
    path_handle_t ref_path = ref_paths.front();
    if (!seen_loci.insert(graph->get_locus_name(ref_path)).second) {
        cerr << "[panpatch]: warning: reference locus " << graph->get_locus_name(ref_path)
             << " appears in more than one input graph; output records will share names" << endl;
    }
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
        cerr << "[panpatch]: skipping " << graph_filename << ": no paths for target sample "
             << sample_names.front() << " (reference " << graph->get_path_name(ref_path) << ")" << endl;
        continue;
    }

    // we patch each target haplotype independently, greedily selecting other haplotypes
    // up front using this simple coverage calculation
    for (const auto& hap_tgts : target_paths) {
        // collect this haplotype's patch records (the PatchRecord table below is the report)
        g_patch_records.clear();

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
            graph, ref_path, hap_tgts.second, sample_names, sample_covers, bed_regions,
            require_telomeres, telo_threshold, max_telomere_patch, progress);

        // Partial-patch cleanup: drop repeat-region-misjoin foreign interior grafts (low k-mer recovery
        // or low flank anchoring), restoring the target's own sequence, while keeping good sub-patches.
        // excised_nonN carries the non-N bp removed per contig so revert_bad_patch's backstop threshold
        // still accounts for them.
        unordered_map<path_handle_t, int64_t> excised_nonN;
        excise_bad_interior_grafts(graph, patched_intervals, sample_names[0], graft_recovery, graft_min_bp,
                                   min_flank, flank_window, excised_nonN);

        // record scaffold joins (patch spans >1 target contig) for the report, and guard foreign-bridged
        // joins: a bridge whose donor doesn't anchor to both contigs' flanks is a wrong-locus misjoin.
        string bridge_detail;
        bool bad_bridge = record_scaffolds(graph, patched_intervals, sample_names[0], min_flank, flank_window, bridge_detail);

        // Check telomere validation if required
        bool telomere_validation_failed = false;
        if (require_telomeres && !patched_intervals.empty()) {
            if (progress) {
                cerr << "[panpatch]: Validating telomeres" << endl;
            }
            bool telomeres_valid = validate_telomeres(graph, patched_intervals, telo_threshold, progress);
            if (!telomeres_valid) {
                telomere_validation_failed = true;
            }
        }

        vector<tuple<step_handle_t, step_handle_t, bool>> input_intervals;
        string revert_reason;
        bool reverted = revert_bad_patch(graph, ref_path, hap_tgts.second, sample_names,
                                         patched_intervals, input_intervals,
                                         default_sample, fail_threshold, graft_recovery, graft_min_bp, telo_threshold, excised_nonN, revert_reason);

        // Also revert if telomere validation failed or a foreign-bridged scaffold is a wrong-locus misjoin
        if (!reverted && (telomere_validation_failed || bad_bridge)) {
            reverted = true;
            revert_reason = bad_bridge ? bridge_detail : "telomere validation failed";
            // Generate input intervals if not already done
            if (input_intervals.empty()) {
                if (!default_sample.empty()) {
                    vector<path_handle_t> default_paths;
                    graph->for_each_path_of_sample(default_sample, [&](path_handle_t path_handle) {
                        size_t hap = graph->get_haplotype(path_handle);
                        if (hap == 0 || hap == PathMetadata::NO_HAPLOTYPE ||
                            hap == graph->get_haplotype(hap_tgts.second.front())) {
                            default_paths.push_back(path_handle);
                        }
                    });
                    if (!default_paths.empty()) {
                        const path_handle_t& def_path = default_paths[0];
                        input_intervals.push_back(make_tuple(graph->path_begin(def_path),
                                                            graph->path_back(def_path), false));
                    }
                } else {
                    for (const path_handle_t& tgt_path : hap_tgts.second) {
                        input_intervals.push_back(make_tuple(graph->path_begin(tgt_path),
                                                            graph->path_back(tgt_path), false));
                    }
                }
            }
        }

        if (reverted) {
            patched_intervals = input_intervals;
        }

        // emit this haplotype's report rows.  Finalize each decision first: a whole-contig revert
        // rejects every patch that had been provisionally accepted on it.
        for (auto& pr : g_patch_records) {
            pr.chrom = graph->get_locus_name(ref_path);
            pr.hap = hap_tgts.first;
            if (reverted && pr.accepted) {
                pr.accepted = false;
                pr.reason = revert_reason.empty() ? "contig reverted to input" : ("contig reverted: " + revert_reason);
            }
            print_patch_row(cout, pr);
        }

        // log telomere information for contigs
        log_contig_telomeres(graph, patched_intervals, telo_threshold);

        // emit the patched-assembly intervals to the BED buffer (the report comments stay on stdout)
        bed_ss << "#Patched assembly on " << graph->get_locus_name(ref_path) << " for "
             << graph->get_sample_name(hap_tgts.second.front()) << "#"
             << graph->get_haplotype(hap_tgts.second.front()) << ":" << endl;
        if (reverted) {
            // reverted output is a set of separate input contigs (each its own FASTA record below),
            // so print each on its own: print_intervals closes the last interval of each call, which
            // for a single-contig list means the full contig (BED then matches the FASTA / #Contig
            // lengths instead of dropping each non-last contig's final node)
            for (const auto& interval : patched_intervals) {
                print_intervals(graph, {interval}, bed_ss);
            }
        } else {
            print_intervals(graph, patched_intervals, bed_ss);
        }
        bed_ss << endl;

        // stream the FASTA to this haplotype's temp file (renamed to the final name only after every
        // graph succeeds, so a failure leaves no partial FASTA -- without holding the genome in RAM)
        if (!out_fasta_filename.empty()) {
            ofstream& fo = fasta_tmp[hap_tgts.first];
            if (!fo.is_open()) {
                string tmp = fasta_name(out_fasta_filename, hap_tgts.first) + ".tmp";
                fo.open(tmp);
                if (!fo) { cerr << "[panpatch] error: Unable to open fasta file for writing: " << tmp << endl; return 1; }
            }
            auto emit = [&](const string& name, const string& seq) {
                fo << ">" << name << "\n";
                for (size_t w = 0; w < seq.length(); w += fasta_width)
                    fo << seq.substr(w, min(fasta_width, seq.length() - w)) << "\n";
            };
            if (reverted) {
                // either we reverted to the original contigs (write each out), or we made a single t2t patch
                for (const auto& interval : patched_intervals) {
                    path_handle_t interval_path = graph->get_path_handle_of_step(get<0>(interval));
                    emit(graph->get_path_name(interval_path), intervals_to_sequence(graph, {interval}));
                }
            } else {
                emit(graph->get_locus_name(ref_path) + "_hap_" + std::to_string(hap_tgts.first),
                     intervals_to_sequence(graph, patched_intervals));
            }
        }
    }    // end haplotype loop
    }    // end per-graph loop

    // Every graph succeeded.  Finalize all outputs atomically: flush/close the temp files, then rename
    // them to their final names at the very end (so a failure during the run leaves only *.tmp, never a
    // partial final BED/FASTA, and the set of finals appears together).
    if (!out_bed_filename.empty()) {
        string bed_tmp = out_bed_filename + ".tmp";
        ofstream ob(bed_tmp);
        if (!ob) { cerr << "[panpatch] error: Unable to open bed file for writing: " << bed_tmp << endl; return 1; }
        ob << bed_ss.str();
        ob.flush();
        if (!ob) { cerr << "[panpatch] error: failed to write bed file (disk full?): " << bed_tmp << endl; return 1; }
    }
    for (auto& kv : fasta_tmp) {
        kv.second.flush();
        if (!kv.second) { cerr << "[panpatch] error: failed to write fasta for hap " << kv.first << " (disk full?)" << endl; return 1; }
        kv.second.close();
    }
    // rename temps -> finals
    if (!out_bed_filename.empty()) {
        if (rename((out_bed_filename + ".tmp").c_str(), out_bed_filename.c_str()) != 0) {
            cerr << "[panpatch] error: failed to finalize bed file: " << out_bed_filename << endl; return 1;
        }
        if (progress) cerr << "[panpatch]: Wrote " << out_bed_filename << endl;
    }
    for (const auto& kv : fasta_tmp) {
        string fn = fasta_name(out_fasta_filename, kv.first);
        if (rename((fn + ".tmp").c_str(), fn.c_str()) != 0) {
            cerr << "[panpatch] error: failed to finalize fasta file: " << fn << endl; return 1;
        }
        if (progress) cerr << "[panpatch]: Wrote " << fn << endl;
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
