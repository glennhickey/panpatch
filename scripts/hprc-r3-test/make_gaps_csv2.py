#!/usr/bin/env python3
# Classify every verkko N-gap (CHM13 deliverable) as filled / rejected / no_donor, precisely:
#   filled    -- the gap is not in a kept verkko interval (an accepted graft replaced it)
#   rejected  -- the gap is retained but lies under a graft that was rejected/reverted
#                (using the report's gap-fill target_start/target_end coordinates)
#   no_donor  -- retained and no graft reached it
import glob, os, re, csv, subprocess, sys, collections
DEL = "deliverable/chm13"
def short(c):
    m = re.search(r"(haplotype(\d)-\d+)", c); return (m.group(1), m.group(2)) if m else (c, "?")

# kept verkko ranges + rejected-graft ranges, keyed by (sample, chrom, contig-short)
kept = collections.defaultdict(list); rej = collections.defaultdict(list); samples = set()
for bed in glob.glob(DEL + "/*.bed"):
    sample = os.path.basename(bed).split(".")[0]; samples.add(sample); chrom = None
    for L in open(bed):
        m = re.match(r"#Patched assembly on (chr\S+) for", L)
        if m: chrom = m.group(1); continue
        m = re.match(r"(verkko\S+)\t(\d+)\t(\d+)\t", L)
        if m and chrom: kept[(sample, chrom, short(m.group(1))[0])].append((int(m.group(2)), int(m.group(3))))
for rep in glob.glob(DEL + "/*.report"):
    sample = os.path.basename(rep).split(".")[0]
    for L in open(rep):
        if L.startswith("#") or L.startswith("chrom\t"): continue
        c = L.rstrip("\n").split("\t")
        if len(c) >= 15 and c[2] == "gap-fill" and c[11] == "rejected" and c[13] != "." and c[14] != ".":
            rej[(sample, c[0], short(c[3])[0])].append((int(c[13]), int(c[14])))

def overlaps(gs, gl, ranges):
    return any(gs < b and a < gs + gl for a, b in ranges)

rows = []
for sample in sorted(samples):
    for graph in sorted(glob.glob("%s.chroms/chr*.full.vg" % sample)):
        chrom = os.path.basename(graph)[:-len(".full.vg")]
        p = subprocess.Popen(["vg", "paths", "-F", "-x", graph, "-Q", "verkko-R3-%s" % sample],
                             stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
        name = None; seq = []
        def flush(nm, parts):
            if not nm: return
            s = "".join(parts).upper(); sc, hap = short(nm)
            kr = kept.get((sample, chrom, sc)); rr = rej.get((sample, chrom, sc), [])
            for mm in re.finditer(r"N+", s):
                gs, gl = mm.start(), mm.end() - mm.start()
                if kr is None:                 status = "dropped"
                elif not overlaps(gs, gl, kr): status = "filled"
                elif overlaps(gs, gl, rr):     status = "rejected"
                else:                          status = "no_donor"
                rows.append((sample, chrom, hap, nm, gs, gl, status))
        for L in p.stdout:
            if L.startswith(">"): flush(name, seq); name = L[1:].split()[0]; seq = []
            else: seq.append(L.strip())
        flush(name, seq); p.wait()
        sys.stderr.write("scanned %s %s\n" % (sample, chrom))

with open("gaps.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["sample", "chrom", "hap", "contig", "gap_start", "gap_len", "status"])
    for r in rows: w.writerow(r)
c = collections.Counter(r[6] for r in rows)
print("wrote gaps.csv (%d gaps): %s" % (len(rows), dict(c)))
