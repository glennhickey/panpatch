#!/usr/bin/env python3
# Build scaffolds.csv (multi-verkko-contig haplotype-chromosomes: joined vs left fragmented) from the run beds.
import glob, os, re, csv
RUNS = "patch-runs"
def ckey(c):
    c = c[3:]; return {"X": 23, "Y": 24, "M": 25}.get(c, int(c) if c.isdigit() else 99)

ev = []
for bed in glob.glob(RUNS + "/*.bed"):
    sample, chrom = os.path.basename(bed)[:-4].split(".", 1)
    for block in open(bed).read().split("\n\n"):
        m = re.search(r"#Patched assembly on \S+ for \S+#(\d+):", block)
        if not m: continue
        vk = [int(L) for nm, L in re.findall(r"#Contig (verkko\S+) len=(\d+)bp", block)]
        if len(vk) < 2: continue
        reverted = ("#Reverting" in block) or ("#No patching is required" in block) or ("#Telomere validation failed" in block)
        ev.append((sample, chrom, m.group(1), "fragmented" if reverted else "joined", sorted(vk, reverse=True)))

with open("scaffolds.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["sample", "chrom", "hap", "status", "n_contigs", "backbone_bp", "fragment_bp(desc)"])
    for s, c, h, st, sz in sorted(ev, key=lambda e: (e[0], ckey(e[1]), e[2])):
        w.writerow([s, c, h, st, len(sz), sz[0], ";".join(map(str, sz[1:]))])
nj = sum(e[3] == "joined" for e in ev)
print("wrote scaffolds.csv (%d multi-contig haps: %d joined, %d fragmented)" % (len(ev), nj, len(ev) - nj))
