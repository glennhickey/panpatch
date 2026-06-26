#!/usr/bin/env python3
# List every kept patch (target + donor contigs with exact coordinate ranges) from the run beds.
# Writes patches.md (markdown, thousands-separated) and patches.csv (plain integers).
import glob, os, re, csv
def ckey(c):
    c = c[3:]; return {"X": 23, "Y": 24, "M": 25}.get(c, int(c) if c.isdigit() else 99)
def fmt(lst, commas):
    out = []
    for nm, a, b, st in lst:
        f = nm.split("#"); hap, contig = f[1], f[2]
        label = contig if nm.startswith("verkko") else "h%s %s" % (hap, contig)   # donor: tag hifiasm hap
        A = format(int(a), ",") if commas else a
        B = format(int(b), ",") if commas else b
        out.append("%s:%s-%s%s" % (label, A, B, "(-)" if st == "-" else ""))
    return "; ".join(out)

rows = []
for bed in glob.glob("patch-runs/*.bed"):
    sample, chrom = os.path.basename(bed)[:-4].split(".", 1)
    for block in open(bed).read().split("\n\n"):
        m = re.search(r"#Patched assembly on \S+ for \S+#(\d+):", block)
        if not m: continue
        hap = m.group(1)
        if re.search(r"#Reverting|#No patching is required|#Telomere validation failed", block): continue
        ivs = re.findall(r"^((?:verkko|hifiasm)\S+)\t(\d+)\t(\d+)\t([+-])", block, re.M)
        if len(ivs) < 2: continue
        typ = "telomere" if "#Telomere patch" in block else ("gap-fill" if any("hifiasm" in i[0] for i in ivs) else "scaffold")
        gg = re.search(r"grafted=(\d+)bp", block); mm = re.search(r"kmer_recovery=([\d.]+)%", block)
        graft = gg.group(1) if gg else ""; rec = mm.group(1) if mm else ""
        tgt = [i for i in ivs if i[0].startswith("verkko")]; dnr = [i for i in ivs if i[0].startswith("hifiasm")]
        if typ == "scaffold":
            tgt = sorted(tgt, key=lambda i: int(i[2]) - int(i[1]), reverse=True); dnr = tgt[1:]; tgt = tgt[:1]
        rows.append((sample, chrom, hap, typ, tgt, dnr, graft, rec))
rows.sort(key=lambda r: ({"telomere": 0, "scaffold": 1, "gap-fill": 2}[r[3]], r[0], ckey(r[1])))

with open("patches.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["sample", "chrom", "hap", "type", "target", "donor_or_joined", "graft_bp", "kmer_recovery_pct"])
    for s, c, h, t, tg, dn, g, r in rows:
        w.writerow([s, c, h, t, fmt(tg, False), fmt(dn, False), g, r])

with open("patches.md", "w") as f:
    hdr = ["sample", "chr·hap", "type", "target  (contig : range)", "donor / joined  (contig : range)", "graft / recovery"]
    f.write("| " + " | ".join(hdr) + " |\n|" + "|".join(["---"] * len(hdr)) + "|\n")
    for s, c, h, t, tg, dn, g, r in rows:
        met = (("graft %s bp" % format(int(g), ",")) if g else "") + ((" / rec %s%%" % r) if r else "")
        f.write("| %s | %s·h%s | %s | %s | %s | %s |\n" % (s, c, h, t, fmt(tg, True), fmt(dn, True), met))

print("wrote patches.csv and patches.md  (%d kept patches: %d telomere, %d scaffold, %d gap-fill)" % (
    len(rows), sum(r[3] == "telomere" for r in rows), sum(r[3] == "scaffold" for r in rows), sum(r[3] == "gap-fill" for r in rows)))
