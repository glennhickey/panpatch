#!/usr/bin/env python3
# Build the slide table + the CHM13-vs-self-ref diff from the NEW deliverable TSV reports.
#   chm13_patches.md          -- the report output (accepted CHM13 patches), for the slide table
#   diff_chm13_vs_selfref.txt -- per-sample T2T comparison + the cases that differ
import glob, os, re, csv, collections

CH = "deliverable/chm13"; SR = "deliverable/self-ref"
def ckey(c):
    c = c[3:]; return {"X": 23, "Y": 24, "M": 25}.get(c, int(c) if c.isdigit() else 99)
def hap_of(contig):          # haplotypeN-... -> N
    m = re.search(r"haplotype(\d)", contig); return int(m.group(1)) if m else 0
def short(contig):           # ...#haplotypeN-NNNN(#0)? -> haplotypeN-NNNN
    m = re.search(r"(haplotype\d-\d+)", contig); return m.group(1) if m else contig

# --- c2chrom: verkko contig -> chromosome, from the CHM13 BED (#Patched headers give the chrom) ---
c2chrom = {}
for bed in glob.glob(CH + "/*.bed"):
    sample = os.path.basename(bed).split(".")[0]; chrom = None
    for L in open(bed):
        m = re.match(r"#Patched assembly on (chr\S+) for", L)
        if m: chrom = m.group(1); continue
        m = re.match(r"(verkko\S+)\t", L)
        if m and chrom: c2chrom[(sample, short(m.group(1)))] = chrom

# --- parse a report: TSV rows + #Contig cap lines ---
def parse_report(path):
    rows, caps = [], {}     # caps: contig -> (left bool, right bool)
    for L in open(path):
        if L.startswith("#Contig"):
            m = re.search(r"#Contig (\S+) len=\d+bp left=(YES|NO)\S* right=(YES|NO)", L)
            if m: caps[m.group(1)] = (m.group(2) == "YES", m.group(3) == "YES")
        elif not L.startswith("#") and not L.startswith("chrom\t"):
            c = L.rstrip("\n").split("\t")
            if len(c) >= 13: rows.append(c)
    return rows, caps

# contigs that keep an N-gap in a mode (gap is retained in a kept verkko interval, i.e. not replaced)
def mode_unfilled(mode):
    kept = collections.defaultdict(list)
    beds = glob.glob(CH + "/*.bed") if mode == "chm13" else glob.glob(SR + "/*/out.bed")
    for bed in beds:
        sample = os.path.basename(bed).split(".")[0] if mode == "chm13" else os.path.basename(os.path.dirname(bed)).split(".hap")[0]
        for L in open(bed):
            m = re.match(r"(verkko\S+)\t(\d+)\t(\d+)\t", L)
            if m: kept[(sample, short(m.group(1)))].append((int(m.group(2)), int(m.group(3))))
    def ov(gs, gl, rngs): return any(gs < b and a < gs + gl for a, b in rngs)
    unf = set()
    for r in csv.DictReader(open("gaps.csv")):
        sample, sc = r["sample"], short(r["contig"])
        kr = kept.get((sample, sc))
        if kr and ov(int(r["gap_start"]), int(r["gap_len"]), kr): unf.add((sample, sc))
    return unf

# gap-free T2T: for a (sample,chrom,hap) is there a contig with both telomeres AND no remaining N-gap?
def t2t_from_caps(sample, caps, unfilled):
    out = collections.defaultdict(bool)
    for contig, (l, r) in caps.items():
        key = (sample, c2chrom.get((sample, short(contig)), "?"), hap_of(contig))
        if (l and r) and (sample, short(contig)) not in unfilled: out[key] = True
        else: out.setdefault(key, False)
    return out

samples = ["HG00126", "HG00235", "HG01074"]

# ---- CHM13: accepted-patch table + T2T ----
chm_unf = mode_unfilled("chm13")
patch_rows = []; chm_t2t = collections.defaultdict(bool)
for s in samples:
    rows, caps = parse_report("%s/%s.report" % (CH, s))
    for k, v in t2t_from_caps(s, caps, chm_unf).items(): chm_t2t[k] = chm_t2t[k] or v
    for c in rows:
        if c[11] == "accepted":
            patch_rows.append((s, c[0], c[1], c[2], short(c[3]), short(c[5]), c[7], c[8]))  # sample chrom hap type target donor replaced kmer

patch_rows.sort(key=lambda r: (r[0], ckey(r[1]), r[2]))
with open("chm13_patches.md", "w") as f:
    f.write("| sample | chrom | hap | type | target contig | donor | replaced_bp | kmer% |\n|" + "---|" * 8 + "\n")
    for r in patch_rows:
        f.write("| " + " | ".join(str(x) for x in r) + " |\n")

# ---- self-ref: T2T ----
sr_unf = mode_unfilled("selfref")
sr_t2t = collections.defaultdict(bool)
for d in sorted(glob.glob(SR + "/*.hap*")):
    if not os.path.isdir(d): continue
    s = os.path.basename(d).split(".hap")[0]
    rows, caps = parse_report(d + "/report")
    for k, v in t2t_from_caps(s, caps, sr_unf).items(): sr_t2t[k] = sr_t2t[k] or v

# ---- diff ----
keys = sorted(set(chm_t2t) | set(sr_t2t), key=lambda k: (k[0], ckey(k[1]) if k[1].startswith("chr") else 999, k[2]))
keys = [k for k in keys if k[1] != "?"]   # drop unplaced (chrOther passthrough) from the per-chromosome comparison
with open("diff_chm13_vs_selfref.txt", "w") as f:
    f.write("CHM13 vs reference-free (self) gap-free T2T per haplotype-chromosome (both telomeres AND no N-gaps)\n\n")
    f.write("Per sample (of 46 hap-chromosomes):\n")
    tot_c = tot_s = 0
    for s in samples:
        c = sum(chm_t2t[k] for k in keys if k[0] == s); z = sum(sr_t2t[k] for k in keys if k[0] == s)
        tot_c += c; tot_s += z
        f.write("  %-9s CHM13 %2d   self-ref %2d   (delta %+d)\n" % (s, c, z, z - c))
    f.write("  %-9s CHM13 %2d   self-ref %2d   (delta %+d)\n\n" % ("TOTAL", tot_c, tot_s, tot_s - tot_c))
    f.write("self-ref T2T where CHM13 is not:\n")
    for k in keys:
        if sr_t2t[k] and not chm_t2t[k]: f.write("  %s %s hap%d\n" % k)
    f.write("\nCHM13 T2T where self-ref is not:\n")
    for k in keys:
        if chm_t2t[k] and not sr_t2t[k]: f.write("  %s %s hap%d\n" % k)
print("wrote chm13_patches.md (%d accepted patches) + diff_chm13_vs_selfref.txt" % len(patch_rows))
print(open("diff_chm13_vs_selfref.txt").read())
