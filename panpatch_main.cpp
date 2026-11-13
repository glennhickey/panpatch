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
#include "bdsg/overlays/overlay_helper.hpp"
#include "panpatch.hpp"
#include "rdna.hpp"

//#define debug

using namespace std;
using namespace handlegraph;
using namespace bdsg;


// from hal2vg/clip-vg.cpp
static unique_ptr<PathHandleGraph> load_graph(istream& graph_stream);

// the per-patch report table (TSV) printed to stdout
static const char* PATCH_TABLE_HEADER =
    "chrom\thap\ttype\ttarget\ttarget_bp\tdonor\tdonor_bp\treplaced_bp\tkmer%\tflankL%\tflankR%\tdecision\treason\ttarget_start\ttarget_end";
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
      << (pr.reason.empty() ? "." : pr.reason) << '\t'
      << fmt_bp(pr.target_start) << '\t' << fmt_bp(pr.target_end) << '\n';
}

static const size_t fasta_width = 80;

void help(char** argv) {
  cerr << "usage: " << argv[0] << " [options] <graph> [graph2 ...]" << endl
       << "Patch an assembly toward telomere-to-telomere chromosomes with the other assemblies in a pangenome graph" << endl
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
       << "        --patch-types LIST       Which of {gap,telomere} to attempt, comma-separated [default: both]; scaffolding is always attempted" << endl
       << "    -T, --require-telomeres      Require a telomere at both ends (revert a contig that isn't); implies telomere patching" << endl
       << "    -M, --max-telomere-patch N   Max bp of target sequence a telomere patch, or a tip graft (a donor carrying" << endl
       << "                                 the patch past the target's outermost piece), may replace at a contig end (0:" << endl
       << "                                 none: a graft may only add sequence at the tip) [20000]" << endl
       << "        --max-subtelomere-graft N  Max bp of donor sequence a telomere patch may graft ahead of the donor's" << endl
       << "                                 telomere (a subtelomere the target lacks, which reads rarely vouch for), and a" << endl
       << "                                 tip graft may add beyond the target bp it replaces.  0 disables [2000]" << endl
       << "    -b, --exclude-bed FILE       BED file of target regions to exclude from patching" << endl
       << "        --error-bed FILE         BED of assembly-error intervals on any input assembly, contigs named as the graph" << endl
       << "                                 names its paths (SAMPLE#HAP#CONTIG; HAP 0 for a haploid sample).  The target's" << endl
       << "                                 are replaced by donor sequence where a donor can thread across them (and kept" << endl
       << "                                 where none can); a donor's never go into a patch.  The assemblies are aligned" << endl
       << "                                 unmasked, so a patch differs from the target only where the donor does" << endl
       << "        --min-cover FLOAT        Revert a patch shorter than this fraction of the target contigs it is made from [0.95]" << endl
       << "        --telomere-threshold F   Min telomere hexamer density to call a telomere end (used by the" << endl
       << "                                 -T validation/patch and the always-on discarded-telomere guard) [0.8]" << endl
       << "        --graft-recovery FLOAT   Excise a foreign interior graft sharing less than this % of the replaced k-mers [50]" << endl
       << "        --graft-min-bp N         Apply --graft-recovery only when at least this many non-N bp are replaced [10000]" << endl
       << "        --min-flank FLOAT        Excise a foreign interior graft replacing fewer than --graft-min-bp non-N bp (an" << endl
       << "                                 N-gap fill), or revert a scaffold bridge, anchored to less than this % of the" << endl
       << "                                 target flank on either side [50]" << endl
       << "        --flank-window N         Window (bp) each side of a graft over which flank anchoring is measured [500000]" << endl
       << "        --max-graft-size-change N  Excise a foreign interior graft that makes its contig more than N bp longer" << endl
       << "                                 or shorter: too large a change to verify (satellite arrays, rDNA gaps).  0 disables [100000]" << endl
       << "        --min-graft-complexity F  Excise a foreign interior graft that replaces or inserts >= --graft-complexity-min-bp of" << endl
       << "                                 low-complexity sequence (distinct k-mer fraction < F over all of it, or in" << endl
       << "                                 some --graft-complexity-window), which k-mer recovery cannot check (tandem" << endl
       << "                                 satellites).  0 disables [0.5]" << endl
       << "        --graft-complexity-window N  Also measure --min-graft-complexity in windows of N bp, so a satellite" << endl
       << "                                 array inside a larger graft is not diluted by the sequence around it.  0" << endl
       << "                                 measures the whole graft only [10000]" << endl
       << "        --graft-complexity-min-bp N  Apply --min-graft-complexity only when at least this many non-N bp are" << endl
       << "                                 replaced or inserted [2000]" << endl
       << "        --max-swap-size-change N  Excise a swap (a graft over sequence with no N gap, or an insertion) that" << endl
       << "                                 makes its contig more than N bp longer or shorter: a copy-number change k-mer" << endl
       << "                                 recovery cannot check.  Gap fills are held to --max-graft-size-change only.  0 disables [10000]" << endl
       << "        --no-patch-fasta FILE    More sequences no graft may replace or put in, besides the built-in rDNA" << endl
       << "                                 unit.  A graft whose replaced or inserted sequence, as a whole or in any 10kb" << endl
       << "                                 window, shares >= --no-patch-kmer-fraction of its (distinct, sampled) k-mers" << endl
       << "                                 with them is not made: an interior graft is excised, a tip graft dropped (the" << endl
       << "                                 target keeps its own end), and a scaffold bridge's scaffold reverted.  Telomere" << endl
       << "                                 grafts are not checked" << endl
       << "        --no-patch-kmer-fraction F  [0.3]" << endl
       << "        --patch-rdna             Let grafts take rDNA: the human rDNA unit (KY962518.1, built in) is otherwise" << endl
       << "                                 a no-patch sequence, its near-identical copies being impossible to verify" << endl
       << "        --max-graft-pieces N     Excise a graft (interior, tip or scaffold bridge) spliced from more than N" << endl
       << "                                 pieces of donor sequence: two donor contigs, or one on both strands or out of" << endl
       << "                                 its order, which no donor has in one piece.  0 disables [1]" << endl
       << "        --telomere-report FILE   Standalone: print per-contig telomere cap status (#Contig lines) for" << endl
       << "                                 a finished FASTA and exit. No graph or -r/-s needed. Use this to" << endl
       << "                                 measure an assembly that was post-processed after panpatch ran" << endl
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
    string error_bed_filename;
    string telomere_report_filename;
    int c;
    int64_t window_size = 1000;
    bool require_telomeres = false;
    bool patch_gap = true;
    bool patch_telomere = true;
    int64_t max_telomere_patch = 20000;
    int64_t max_subtelomere_graft = 2000;
    double fail_threshold = 0.95;
    double telo_threshold = 0.8;
    double graft_recovery = 50.0;
    int64_t graft_min_bp = 10000;
    double min_flank = 50.0;
    int64_t flank_window = 500000;
    int64_t max_graft_size_change = 100000;
    double min_graft_complexity = 0.5;
    int64_t graft_complexity_window = 10000;
    int64_t graft_complexity_min_bp = 2000;
    int64_t max_graft_pieces = 1;
    int64_t max_swap_size_change = 10000;
    string no_patch_fasta;
    double no_patch_kmer_fraction = 0.3;
    bool patch_rdna = false;
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
            {"patch-types", required_argument, 0, 1008},
            {"max-telomere-patch", required_argument, 0, 'M'},
            {"exclude-bed", required_argument, 0, 'b'},
            {"min-cover", required_argument, 0, 1001},
            {"telomere-threshold", required_argument, 0, 1002},
            {"graft-recovery", required_argument, 0, 1003},
            {"graft-min-bp", required_argument, 0, 1004},
            {"min-flank", required_argument, 0, 1005},
            {"flank-window", required_argument, 0, 1006},
            {"bed", required_argument, 0, 1007},
            {"telomere-report", required_argument, 0, 1009},
            {"max-graft-size-change", required_argument, 0, 1010},
            {"max-subtelomere-graft", required_argument, 0, 1012},
            {"min-graft-complexity", required_argument, 0, 1011},
            {"error-bed", required_argument, 0, 1013},
            {"graft-complexity-window", required_argument, 0, 1014},
            {"graft-complexity-min-bp", required_argument, 0, 1015},
            {"max-graft-pieces", required_argument, 0, 1016},
            {"max-swap-size-change", required_argument, 0, 1017},
            {"no-patch-fasta", required_argument, 0, 1018},
            {"no-patch-kmer-fraction", required_argument, 0, 1019},
            {"patch-rdna", no_argument, 0, 1020},
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
            if (window_size <= 0) {
                cerr << "[panpatch] error: --window (-w) must be a positive integer" << endl;
                return 1;
            }
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
        case 1008:
        {
            patch_gap = patch_telomere = false;
            stringstream pss(optarg); string ptok;
            while (getline(pss, ptok, ',')) {
                if (ptok.empty()) continue;
                else if (ptok == "gap") patch_gap = true;
                else if (ptok == "telomere") patch_telomere = true;
                else { cerr << "[panpatch] error: unknown --patch-types value '" << ptok
                            << "' (expected a comma-separated subset of gap,telomere)" << endl; return 1; }
            }
            if (!patch_gap && !patch_telomere) {
                cerr << "[panpatch] error: --patch-types selected no patch types" << endl; return 1;
            }
            break;
        }
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
        case 1013:
            error_bed_filename = optarg;
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
        case 1010:
            max_graft_size_change = strtoll(optarg, nullptr, 10);
            break;
        case 1012:
            max_subtelomere_graft = strtoll(optarg, nullptr, 10);
            break;
        case 1011:
            min_graft_complexity = atof(optarg);
            break;
        case 1014:
            graft_complexity_window = strtoll(optarg, nullptr, 10);
            break;
        case 1015:
            graft_complexity_min_bp = strtoll(optarg, nullptr, 10);
            break;
        case 1016:
            max_graft_pieces = strtoll(optarg, nullptr, 10);
            break;
        case 1017:
            max_swap_size_change = strtoll(optarg, nullptr, 10);
            break;
        case 1018:
            no_patch_fasta = optarg;
            break;
        case 1019:
            no_patch_kmer_fraction = atof(optarg);
            break;
        case 1020:
            patch_rdna = true;
            break;
        case 1007:
            out_bed_filename = optarg;
            break;
        case 1009:
            telomere_report_filename = optarg;
            break;
        case 'h':
            help(argv);
            return 0;
        case '?':
            /* getopt_long already printed an error message. */
            help(argv);
            return 1;
        default:
            abort ();
        }
    }

    if (argc <= 1) {
        help(argv);
        return 1;
    }
    // --telomere-report: standalone mode.  Scan a finished FASTA, print its per-contig telomere cap
    // status ("#Contig" lines, the same ones a normal run prints), and exit.  This is how a caller
    // measures the assembly it actually delivered when that differs from what panpatch emitted (after
    // post-processing it, say), which the in-graph measurement cannot see.  No graph or samples needed.
    if (!telomere_report_filename.empty()) {
        if (telo_threshold < 0 || telo_threshold > 1) {
            cerr << "[panpatch] error: --telomere-threshold must be in [0,1]" << endl;
            return 1;
        }
        return telomere_report_fasta(telomere_report_filename, telo_threshold);
    }
    if (sample_names.empty()) {
        cerr << "[panpatch] error: -s must be used to specify at least one sample name to prioritize" << endl;
        return 1;
    }
    if (ref_sample.empty()) {
        cerr << "[panpatch] error: -r must be used to specify a reference sample" << endl;
        return 1;
    }
    // -T requires the T2T end result, which is meaningless if telomere patching is turned off
    if (require_telomeres && !patch_telomere) {
        cerr << "[panpatch] error: -T/--require-telomeres cannot be used with --patch-types that excludes telomere" << endl;
        return 1;
    }
    // validate numeric thresholds (a nonsensical value silently disables or inverts a guard)
    if (fail_threshold < 0 || fail_threshold > 1) { cerr << "[panpatch] error: --min-cover must be in [0,1]" << endl; return 1; }
    if (telo_threshold < 0 || telo_threshold > 1) { cerr << "[panpatch] error: --telomere-threshold must be in [0,1]" << endl; return 1; }
    if (graft_recovery < 0 || graft_recovery > 100) { cerr << "[panpatch] error: --graft-recovery must be in [0,100]" << endl; return 1; }
    if (min_flank < 0 || min_flank > 100) { cerr << "[panpatch] error: --min-flank must be in [0,100]" << endl; return 1; }
    if (graft_min_bp < 0) { cerr << "[panpatch] error: --graft-min-bp must be >= 0" << endl; return 1; }
    if (flank_window <= 0) { cerr << "[panpatch] error: --flank-window must be > 0" << endl; return 1; }
    if (max_graft_size_change < 0) { cerr << "[panpatch] error: --max-graft-size-change must be >= 0" << endl; return 1; }
    if (min_graft_complexity < 0 || min_graft_complexity > 1) { cerr << "[panpatch] error: --min-graft-complexity must be in [0,1]" << endl; return 1; }
    if (graft_complexity_window < 0) { cerr << "[panpatch] error: --graft-complexity-window must be >= 0" << endl; return 1; }
    if (graft_complexity_min_bp < 0) { cerr << "[panpatch] error: --graft-complexity-min-bp must be >= 0" << endl; return 1; }
    if (max_graft_pieces < 0) { cerr << "[panpatch] error: --max-graft-pieces must be >= 0" << endl; return 1; }
    if (max_swap_size_change < 0) { cerr << "[panpatch] error: --max-swap-size-change must be >= 0" << endl; return 1; }
    if (!(no_patch_kmer_fraction > 0 && no_patch_kmer_fraction <= 1)) { cerr << "[panpatch] error: --no-patch-kmer-fraction must be in (0,1]" << endl; return 1; }
    NoPatch no_patch;
    no_patch.fraction = no_patch_kmer_fraction;
    if (!patch_rdna) {
        add_no_patch_sequence(RDNA_UNIT_KY962518_1, no_patch);
    }
    if (!no_patch_fasta.empty()) {
        NoPatch more;
        int64_t n = load_no_patch_fasta(no_patch_fasta, more);
        if (n < 0) {
            cerr << "[panpatch] error: cannot open --no-patch-fasta " << no_patch_fasta << endl;
            return 1;
        }
        if (more.kmers.empty()) {
            cerr << "[panpatch] error: no sequence (FASTA record of at least 31 bases) in --no-patch-fasta " << no_patch_fasta << endl;
            return 1;
        }
        no_patch.kmers.insert(more.kmers.begin(), more.kmers.end());
    }
    if (progress && !no_patch.kmers.empty()) {
        cerr << "[panpatch]: " << no_patch.kmers.size() << " sampled k-mers of no-patch sequence ("
             << (patch_rdna ? "" : "the rDNA unit") << (!patch_rdna && !no_patch_fasta.empty() ? ", " : "")
             << no_patch_fasta << ")" << endl;
    }
    const NoPatch* no_patch_p = no_patch.kmers.empty() ? nullptr : &no_patch;
    if (max_telomere_patch < 0) { cerr << "[panpatch] error: --max-telomere-patch must be >= 0" << endl; return 1; }
    if (max_subtelomere_graft < 0) { cerr << "[panpatch] error: --max-subtelomere-graft must be >= 0" << endl; return 1; }
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

    BedRegions error_bed;
    if (!error_bed_filename.empty()) {
        error_bed = parse_bed_file(error_bed_filename);
        if (progress) {
            int64_t total_regions = 0;
            for (const auto& br : error_bed) total_regions += br.second.size();
            cerr << "[panpatch]: Loaded " << total_regions << " error intervals on " << error_bed.size()
                 << " contigs from " << error_bed_filename << endl;
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
        // an --error-bed contig whose sample is in no graph is almost certainly misnamed (its intervals
        // would silently do nothing).  a contig on a chromosome these graphs don't hold is normal
        set<string> unknown;
        for (const auto& br : error_bed) {
            string s = br.first.substr(0, br.first.find('#'));
            if (!present_samples.count(s)) unknown.insert(s);
        }
        for (const string& s : unknown) {
            cerr << "[panpatch] warning: --error-bed names sample " << s << ", which is in no input graph" << endl;
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
    // Index every path: the reference, targets and donors are all haplotype paths, and the overlay
    // leaves haplotypes out by default (they then do not exist in it at all)
    bdsg::ReferencePathOverlayHelper overlay_helper;
    PathPositionHandleGraph* graph = overlay_helper.apply(base_graph.get(), true);
    if (progress && dynamic_cast<PathPositionHandleGraph*>(base_graph.get()) == nullptr) {
        cerr << "[panpatch]: Applied position overlay" << endl;
    }

    // pull out the (one and only) reference path.  paths are always taken by name (sample_paths), never
    // in the graph's own order, which varies between graph implementations and libbdsg versions
    vector<path_handle_t> ref_paths = sample_paths(graph, ref_sample);
    if (ref_paths.size() != 1) {
        // No single reference path -- e.g. a chrOther / unplaced-contigs graph, which holds many of the
        // target's own contigs rather than one chromosome.  Rather than drop it, pass the target sample's
        // contigs through unchanged so every input contig still reaches the output.
        cerr << "[panpatch]: warning: " << graph_filename << " has " << ref_paths.size()
             << " reference paths for " << ref_sample << " (not a single chromosome); passing its "
             << sample_names.front() << " contigs through unchanged" << endl;
        string stem = graph_filename;
        { size_t s = stem.find_last_of('/'); if (s != string::npos) stem = stem.substr(s + 1);
          size_t d = stem.find(".full.vg"); if (d == string::npos) d = stem.rfind(".vg");
          if (d != string::npos) stem = stem.substr(0, d); }
        vector<pair<string, string>> unplaced_records;
        for (path_handle_t p : sample_paths(graph, sample_names.front())) {
            int64_t hap = graph->get_haplotype(p);
            vector<tuple<step_handle_t, step_handle_t, bool>> whole = { make_tuple(graph->path_begin(p), graph->path_back(p), false) };
            string seq = intervals_to_sequence(graph, whole);
            if (!out_fasta_filename.empty()) {
                ofstream& fo = fasta_tmp[hap];
                if (!fo.is_open()) {
                    string tmp = fasta_name(out_fasta_filename, hap) + ".tmp";
                    fo.open(tmp);
                    if (!fo) { cerr << "[panpatch] error: Unable to open fasta file for writing: " << tmp << endl; return 1; }
                }
                fo << ">" << graph->get_path_name(p) << "\n";
                for (size_t w = 0; w < seq.length(); w += fasta_width)
                    fo << seq.substr(w, min(fasta_width, seq.length() - w)) << "\n";
            }
            bed_ss << "#Passthrough (unplaced) for " << graph->get_sample_name(p) << "#" << hap << ":" << endl;
            print_intervals(graph, whole, bed_ss);
            bed_ss << endl;
            PatchRecord pr;
            pr.chrom = stem; pr.hap = hap; pr.type = "passthrough";
            pr.target = graph->get_path_name(p); pr.target_bp = path_bp(graph, p);
            pr.accepted = true; pr.reason = "unplaced (no reference chromosome)";
            print_patch_row(cout, pr);
            unplaced_records.push_back(make_pair(graph->get_path_name(p), std::move(seq)));
        }
        // their "#Contig" cap lines, after the rows as for every other graph
        for (const auto& rec : unplaced_records) {
            print_contig_telomere_line(rec.first, rec.second, telo_threshold);
        }
        continue;
    }
    path_handle_t ref_path = ref_paths.front();
    // the --error-bed intervals on this graph's paths
    ErrorRegions error_regions = error_bed_to_paths(graph, error_bed);
    if (progress && !error_bed.empty()) {
        int64_t n = 0, bp = 0;
        for (const auto& er : error_regions) for (const auto& iv : er.second) { ++n; bp += iv.second - iv.first; }
        cerr << "[panpatch]: " << n << " error intervals (" << bp << " bp) on " << error_regions.size()
             << " paths of " << graph_filename << endl;
    }
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
        for (const path_handle_t& path_handle : sample_paths(graph, sample)) {
            other_paths.push_back(path_handle);
        }
    }
    if (progress) {
        cerr << "[panpatch]: Selected " << other_paths.size() << " (non-reference) paths" << endl;
    }

    // pull out the target paths (the ones we want to patch) and sort them by haplotype (and then by
    // name: this is the order the target's contigs are ranked in at a fork, and written out in when a
    // patch is reverted)
    map<int64_t, vector<path_handle_t>> target_paths;
    for (const path_handle_t& path_handle : sample_paths(graph, sample_names.front())) {
        target_paths[graph->get_haplotype(path_handle)].push_back(path_handle);
    }
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

            for (const path_handle_t& p : other_paths) {
                auto cov = coverage_map.find(p);
                if (cov != coverage_map.end() && graph->get_sample_name(p) != sample_names.front()) {
                    cerr << "[panpatch]:    " << graph->get_path_name(p) << " " << cov->second << endl;
                }
            }
        }

        unordered_map<string, vector<path_handle_t>> sample_covers = select_sample_covers(graph, coverage_map);

        if (progress) {
            cerr << "[panpatch]: Sample cover selection:" << endl;
            set<string> printed;
            for (const string& sample : sample_names) {
                auto sc = sample_covers.find(sample);
                if (sc == sample_covers.end() || !printed.insert(sample).second) {
                    continue;
                }
                cerr << "[panpatch]:    " << sc->first;
                for (const auto& p : sc->second) {
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
        unordered_set<path_handle_t> threaded_tgt_paths;
        vector<tuple<step_handle_t, step_handle_t, bool>> patched_intervals = greedy_patch(
            graph, ref_path, hap_tgts.second, sample_names, sample_covers, bed_regions,
            patch_telomere, telo_threshold, max_telomere_patch, max_subtelomere_graft, progress, require_telomeres,
            &threaded_tgt_paths, &error_regions, patch_gap, max_graft_pieces, no_patch_p);

        // Partial-patch cleanup: drop repeat-region-misjoin foreign interior grafts (low k-mer recovery
        // or low flank anchoring), restoring the target's own sequence, while keeping good sub-patches.
        // excised_nonN carries the non-N bp removed per contig so revert_bad_patch's backstop threshold
        // still accounts for them.
        unordered_map<path_handle_t, int64_t> excised_nonN;
        excise_bad_interior_grafts(graph, patched_intervals, sample_names[0], graft_recovery, graft_min_bp,
                                   min_flank, flank_window, patch_gap, excised_nonN,
                                   max_graft_size_change, min_graft_complexity, graft_complexity_window,
                                   graft_complexity_min_bp, max_graft_pieces, max_swap_size_change, no_patch_p);

        // final consistency guard: a patch that would write any sequence twice or out of order (see
        // patch_reuses_sequence) reverts the contig, whatever else passed
        string reuse_detail;
        bool bad_reuse = !patched_intervals.empty() &&
            patch_reuses_sequence(graph, patched_intervals, sample_names[0], hap_tgts.second, threaded_tgt_paths, reuse_detail);

        // record scaffold joins (patch spans >1 target contig) for the report, and guard foreign-bridged
        // joins: a bridge whose donor doesn't anchor to both contigs' flanks is a wrong-locus misjoin.
        string bridge_detail;
        bool bad_bridge = record_scaffolds(graph, patched_intervals, sample_names[0], min_flank, flank_window, bridge_detail,
                                           max_graft_pieces, no_patch_p);

        // Check telomere validation.  -T requires a telomere at both ends (and no internal one); telomere
        // patching alone (on by default) only rejects an internal-telomere scaffold misjoin (missing ends are kept).
        bool telomere_validation_failed = false;
        string telo_fail_reason;
        if (patch_telomere && !patched_intervals.empty()) {
            if (progress) {
                cerr << "[panpatch]: Validating telomeres" << endl;
            }
            bool telomeres_valid = validate_telomeres(graph, patched_intervals, telo_threshold,
                                                      require_telomeres, progress, telo_fail_reason);
            if (!telomeres_valid) {
                telomere_validation_failed = true;
            }
        }

        vector<tuple<step_handle_t, step_handle_t, bool>> input_intervals;
        string revert_reason;
        bool reverted = revert_bad_patch(graph, ref_path, hap_tgts.second, sample_names,
                                         patched_intervals, input_intervals,
                                         default_sample, fail_threshold, graft_recovery, graft_min_bp, telo_threshold, excised_nonN, revert_reason);

        // Also revert if the patch reuses sequence, telomere validation failed or a foreign-bridged scaffold
        // is a wrong-locus misjoin
        if (!reverted && (bad_reuse || telomere_validation_failed || bad_bridge)) {
            reverted = true;
            revert_reason = bad_reuse ? reuse_detail : bad_bridge ? bridge_detail : telo_fail_reason;
            // Generate input intervals if not already done
            if (input_intervals.empty()) {
                vector<path_handle_t> default_paths;
                if (!default_sample.empty()) {
                    default_paths = default_sample_paths(graph, default_sample, graph->get_haplotype(hap_tgts.second.front()));
                }
                // as in revert_bad_patch: the default sample stands in only with exactly one contig here;
                // otherwise every one of the target's own contigs is kept (none may be dropped)
                if (default_paths.size() == 1) {
                    input_intervals.push_back(make_tuple(graph->path_begin(default_paths[0]),
                                                        graph->path_back(default_paths[0]), false));
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

        // a contig with gaps that no donor could patch reverts with no candidate patch to report it by;
        // say so, rather than writing the contig back with nothing in the report to show it was looked at
        if (reverted && g_patch_records.empty() && revert_reason == NO_DONOR_PATCH_REASON) {
            PatchRecord pr;
            pr.type = "gap-fill";
            pr.target = graph->get_path_name(hap_tgts.second.front());
            pr.target_bp = path_bp(graph, hap_tgts.second.front());
            pr.accepted = false;
            pr.reason = string(NO_DONOR_PATCH_REASON) + "; contig kept as it is";
            g_patch_records.push_back(pr);
        }

        // target contigs on this chromosome that an accepted patch never threaded through at all (an
        // unplaced fragment binned here, say).  they are written through as they are, each its own record,
        // so no input sequence is dropped: patched output is one record built from the patch's intervals,
        // so a contig never in them would otherwise never be written.  (one the threading used and a
        // telomere graft then replaced was patched over, and is not written again)
        vector<path_handle_t> unused_tgts;
        if (!reverted) {
            unordered_set<path_handle_t> used = threaded_tgt_paths;
            for (const auto& iv : patched_intervals) used.insert(graph->get_path_handle_of_step(get<0>(iv)));
            for (const path_handle_t& p : hap_tgts.second) {
                if (!used.count(p)) unused_tgts.push_back(p);
            }
            for (const path_handle_t& p : unused_tgts) {
                PatchRecord pr;
                pr.type = "passthrough";
                pr.target = graph->get_path_name(p); pr.target_bp = path_bp(graph, p);
                pr.accepted = true; pr.reason = "not part of the patched contig; written as it is";
                g_patch_records.push_back(pr);
            }
            // an N gap the patched contig still has, with no row about it, gets one
            report_unfilled_gaps(graph, patched_intervals, sample_names[0]);
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

        // the records this haplotype writes, in output order: the patched record built from the patch's
        // intervals (then any target contig it never threaded, written through), or, when reverted, each
        // input contig as it is.  the #Contig lines, the BED blocks and the FASTA records are all made
        // from this one list, so they describe exactly the same sequences
        struct OutRecord {
            string name;
            string header;
            vector<tuple<step_handle_t, step_handle_t, bool>> intervals;
        };
        vector<OutRecord> records;
        const string locus = graph->get_locus_name(ref_path);
        auto sample_hap = [&](path_handle_t p) {
            return graph->get_sample_name(p) + "#" + std::to_string(graph->get_haplotype(p));
        };
        if (reverted) {
            for (const auto& interval : patched_intervals) {
                path_handle_t p = graph->get_path_handle_of_step(get<0>(interval));
                string label = graph->get_sample_name(p) == sample_names.front()
                    ? "Unpatched (input contig)" : "Default sample contig (" + graph->get_sample_name(p) + ")";
                records.push_back({graph->get_path_name(p),
                                   "#" + label + " on " + locus + " for " + sample_hap(hap_tgts.second.front()) + ":",
                                   {interval}});
            }
        } else {
            records.push_back({locus + "_hap_" + std::to_string(hap_tgts.first),
                               "#Patched assembly on " + locus + " for " + sample_hap(hap_tgts.second.front()) + ":",
                               patched_intervals});
            for (const path_handle_t& p : unused_tgts) {
                records.push_back({graph->get_path_name(p),
                                   "#Passthrough (unused by the patch) on " + locus + " for " + sample_hap(p) + ":",
                                   {make_tuple(graph->path_begin(p), graph->path_back(p), false)}});
            }
        }

        // stream the FASTA to this haplotype's temp file (renamed to the final name only after every
        // graph succeeds, so a failure leaves no partial FASTA -- one chromosome is held at a time)
        ofstream* fo = nullptr;
        if (!out_fasta_filename.empty()) {
            fo = &fasta_tmp[hap_tgts.first];
            if (!fo->is_open()) {
                string tmp = fasta_name(out_fasta_filename, hap_tgts.first) + ".tmp";
                fo->open(tmp);
                if (!*fo) { cerr << "[panpatch] error: Unable to open fasta file for writing: " << tmp << endl; return 1; }
            }
        }
        for (const OutRecord& rec : records) {
            string seq = intervals_to_sequence(graph, rec.intervals);
            // the record's telomere cap status (the report's "#Contig" lines)
            print_contig_telomere_line(rec.name, seq, telo_threshold);
            // its intervals (print_intervals closes the last interval of each call, so each record's
            // block ends on its own final node)
            bed_ss << rec.header << endl;
            print_intervals(graph, rec.intervals, bed_ss);
            bed_ss << endl;
            if (fo) {
                *fo << ">" << rec.name << "\n";
                for (size_t w = 0; w < seq.length(); w += fasta_width)
                    *fo << seq.substr(w, min(fasta_width, seq.length() - w)) << "\n";
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
