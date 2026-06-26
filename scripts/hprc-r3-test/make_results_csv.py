#!/usr/bin/env python3
# Build panpatch_results.csv (per haplotype-chromosome) + panpatch_summary.csv from the CHM13 run beds.
import glob, os, re, csv
RUNS = "patch-runs"
def ckey(c):
    c = c[3:]; return {"X": 23, "Y": 24, "M": 25}.get(c, int(c) if c.isdigit() else 99)

rows = []
for bed in glob.glob(RUNS + "/*.bed"):
    sample, chrom = os.path.basename(bed)[:-4].split(".", 1)
    for block in open(bed).read().split("\n\n"):
        m = re.search(r"#Patched assembly on \S+ for (\S+)#(\d+):", block)
        if not m: continue
        hap = int(m.group(2))
        contigs = re.findall(r"#Contig (\S+) len=\d+bp left=(YES|NO)\([\d.]+\) right=(YES|NO)\([\d.]+\)", block)
        vfail    = "#Telomere validation failed" in block
        reverted = ("#Reverting" in block) or ("#No patching is required" in block)
        telo     = "#Telomere patch (" in block
        kept     = not reverted and not vfail
        after    = True if kept else (len(contigs) == 1 and contigs[0][1] == "YES" and contigs[0][2] == "YES")
        vk = [c for c in contigs if c[0].startswith("verkko")]
        # telomere patches always lacked a clean telomere at that end -> not T2T before
        before = (False if telo else (len(vk) == 1 and vk[0][1] == "YES" and vk[0][2] == "YES")) if kept else after
        pm = re.search(r"#Telomere patch \((front|back)\): donor=(\S+) replaced=(\d+)bp grafted=(\d+)bp kmer_recovery=([\d.]+)%", block)
        if   telo and kept: outcome = "telomere_patch"
        elif kept:          outcome = "scaffold_gap_patch"
        elif telo:          outcome = "telomere_attempt_reverted"
        elif after:         outcome = "already_t2t"
        else:               outcome = "incomplete"
        rows.append(dict(sample=sample, chrom=chrom, hap=hap, t2t_before=int(before), t2t_after=int(after),
            outcome=outcome, patch_end=pm.group(1) if pm else "", donor=pm.group(2) if pm else "",
            replaced_bp=int(pm.group(3)) if pm else "", grafted_bp=int(pm.group(4)) if pm else "",
            kmer_recovery_pct=float(pm.group(5)) if pm else ""))
rows.sort(key=lambda r: (r["sample"], ckey(r["chrom"]), r["hap"]))

cols = ["sample", "chrom", "hap", "t2t_before", "t2t_after", "outcome", "patch_end", "replaced_bp", "grafted_bp", "kmer_recovery_pct", "donor"]
with open("panpatch_results.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
    for r in rows: w.writerow(r)

samples = sorted({r["sample"] for r in rows})
with open("panpatch_summary.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["sample", "hap_chromosomes", "t2t_before", "t2t_after", "gained", "telomere_patches", "scaffold_gap_patches", "incomplete"])
    for s in samples:
        r = [x for x in rows if x["sample"] == s]
        w.writerow([s, len(r), sum(x["t2t_before"] for x in r), sum(x["t2t_after"] for x in r),
                    sum(x["t2t_after"] - x["t2t_before"] for x in r),
                    sum(x["outcome"] == "telomere_patch" for x in r),
                    sum(x["outcome"] == "scaffold_gap_patch" for x in r),
                    sum(x["t2t_after"] == 0 for x in r)])
print("wrote panpatch_results.csv (%d rows) + panpatch_summary.csv" % len(rows))
