#!/usr/bin/env python3
# Build gaps.csv: every N-gap in the verkko targets, classified filled (replaced by a graft) vs
# retained (still in the output), by streaming the verkko contigs and intersecting N-runs with the
# kept verkko intervals in each run bed.
import glob, os, re, csv, subprocess, sys
RUNS = "patch-runs"

# only scan chromosomes that could carry/fill a gap: a donor was used, or a revert happened
chroms = []
for bed in glob.glob(RUNS + "/*.bed"):
    t = open(bed).read()
    if ("hifiasm" in t) or ("#Reverting to input assembly" in t) or ("#Reverting patch" in t):
        chroms.append(os.path.basename(bed)[:-4])

rows = []
for sc in sorted(chroms):
    sample, chrom = sc.split(".", 1)
    graph = "%s.chroms/%s.full.vg" % (sample, chrom)
    kept = {}                                   # verkko contig -> list of kept [a,b) ranges
    for L in open("%s/%s.bed" % (RUNS, sc)):
        m = re.match(r"(verkko\S+)\t(\d+)\t(\d+)\t", L)
        if m: kept.setdefault(m.group(1), []).append((int(m.group(2)), int(m.group(3))))
    p = subprocess.Popen(["vg", "paths", "-F", "-x", graph, "-Q", "verkko-R3-%s" % sample],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    name = None; seq = []
    def flush(nm, parts):
        if not nm: return
        s = "".join(parts).upper(); ranges = kept.get(nm)
        hap = nm.split("#")[1] if "#" in nm else "?"
        for mm in re.finditer(r"N+", s):
            gs, gl = mm.start(), mm.end() - mm.start()
            if ranges is None:                                            status = "dropped"
            elif any(gs < b and a < gs + gl for a, b in ranges):          status = "retained"
            else:                                                         status = "filled"
            rows.append((sample, chrom, hap, nm, gs, gl, status))
    for L in p.stdout:
        if L.startswith(">"):
            flush(name, seq); name = L[1:].strip().split()[0]; seq = []
        else:
            seq.append(L.strip())
    flush(name, seq); p.wait()
    sys.stderr.write("scanned %s\n" % sc)

with open("gaps.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["sample", "chrom", "hap", "contig", "gap_start", "gap_len", "status"])
    for r in rows: w.writerow(r)
fs = [r for r in rows if r[6] == "filled"]; rs = [r for r in rows if r[6] == "retained"]
print("wrote gaps.csv: filled=%d (%.0f kb)  retained=%d (%.2f Mb)" % (
    len(fs), sum(r[5] for r in fs) / 1000, len(rs), sum(r[5] for r in rs) / 1e6))
