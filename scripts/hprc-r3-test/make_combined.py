#!/usr/bin/env python3
# Combined CHM13-vs-reference-free figures: two half-height 3-way bars per chromosome
# (upper = CHM13, lower = reference-free).  panpatch_*_combined.png
import glob, os, re, csv, collections
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

GREEN = "#31a354"; ORANGE = "#e6550d"; GREY = "#cccccc"; LBLUE = "#9ecae1"; DBLUE = "#3182bd"; PURPLE = "#756bb1"
def ckey(c):
    c = c[3:]; return {"X": 23, "Y": 24, "M": 25}.get(c, int(c) if c.isdigit() else 99)
def short(c):
    m = re.search(r"(haplotype\d-\d+)", c); return m.group(1) if m else c
def hapnum(c):
    m = re.search(r"haplotype(\d)", c); return int(m.group(1)) if m else 0
SAMPLES = ["HG00126", "HG00235", "HG01074"]

c2chrom = {}
for bed in glob.glob("deliverable/chm13/*.bed"):
    sample = os.path.basename(bed).split(".")[0]; chrom = None
    for L in open(bed):
        m = re.match(r"#Patched assembly on (chr\S+) for", L)
        if m: chrom = m.group(1); continue
        m = re.match(r"(verkko\S+)\t", L)
        if m and chrom: c2chrom[(sample, short(m.group(1)))] = chrom

def compute(mode):
    if mode == "chm13":
        reps = [(os.path.basename(f).split(".")[0], f) for f in glob.glob("deliverable/chm13/*.report")]
    else:
        reps = [(os.path.basename(os.path.dirname(f)).split(".hap")[0], f) for f in glob.glob("deliverable/self-ref/*/report")]
    rows = []; caps = collections.defaultdict(dict)
    for sample, f in reps:
        for L in open(f):
            if L.startswith("#Contig"):
                m = re.search(r"#Contig (\S+) len=\d+bp left=(YES|NO)\S* right=(YES|NO)", L)
                if m and "verkko" in m.group(1):
                    ct = short(m.group(1)); ch = c2chrom.get((sample, ct), "?")
                    caps[(sample, ch, hapnum(m.group(1)))][ct] = (m.group(2) == "YES", m.group(3) == "YES")
            elif not L.startswith("#") and not L.startswith("chrom\t"):
                c = L.rstrip("\n").split("\t")
                if len(c) < 15: continue
                ch = c[0] if mode == "chm13" else c2chrom.get((sample, short(c[3])), "?")
                rows.append(dict(sample=sample, chrom=ch, hap=hapnum(c[3]), type=c[2], decision=c[11], reason=c[12]))
    # gaps
    kept = collections.defaultdict(list); rej = collections.defaultdict(list)
    beds = glob.glob("deliverable/chm13/*.bed") if mode == "chm13" else glob.glob("deliverable/self-ref/*/out.bed")
    for bed in beds:
        sample = os.path.basename(bed).split(".")[0] if mode == "chm13" else os.path.basename(os.path.dirname(bed)).split(".hap")[0]
        for L in open(bed):
            m = re.match(r"(verkko\S+)\t(\d+)\t(\d+)\t", L)
            if m: kept[(sample, short(m.group(1)))].append((int(m.group(2)), int(m.group(3))))
    for sample, f in reps:
        for L in open(f):
            if L.startswith("#") or L.startswith("chrom\t"): continue
            c = L.rstrip("\n").split("\t")
            if len(c) >= 15 and c[2] == "gap-fill" and c[11] == "rejected" and c[13] != "." and c[14] != ".":
                rej[(sample, short(c[3]))].append((int(c[13]), int(c[14])))
    def ov(gs, gl, r): return any(gs < b and a < gs + gl for a, b in r)
    gaps = collections.defaultdict(lambda: [0, 0, 0]); unfilled = set()
    for r in csv.DictReader(open("gaps.csv")):
        sample, sc, chrom = r["sample"], short(r["contig"]), r["chrom"]
        gs, gl = int(r["gap_start"]), int(r["gap_len"])
        kr = kept.get((sample, sc)); rr = rej.get((sample, sc), [])
        if kr is None: continue
        if not ov(gs, gl, kr): gaps[chrom][0] += 1
        else:
            unfilled.add((sample, sc))            # this contig keeps an N-gap -> not gap-free
            if ov(gs, gl, rr): gaps[chrom][1] += 1
            else:              gaps[chrom][2] += 1
    # scaffolds
    scaff_dec = {}
    for r in rows:
        if r["type"] == "scaffold" and r["chrom"] != "?":
            k = (r["sample"], r["chrom"], r["hap"])
            if scaff_dec.get(k) != "rejected": scaff_dec[k] = "joined" if r["decision"] == "accepted" else "rejected"
    scaff = collections.defaultdict(lambda: [0, 0, 0])
    for (sample, chrom, hap), contigs in caps.items():
        if chrom == "?" or len(contigs) < 2: continue
        frags = len(contigs) - 1; d = scaff_dec.get((sample, chrom, hap))
        scaff[chrom][0 if d == "joined" else 1 if d == "rejected" else 2] += frags
    if mode != "chm13":
        cc = collections.defaultdict(set)
        for s in SAMPLES:
            for L in open("deliverable/chm13/%s.report" % s):
                if L.startswith("#Contig") and "verkko" in L:
                    nm = L.split()[1]; ch = c2chrom.get((s, short(nm)))
                    if ch: cc[(s, ch, hapnum(nm))].add(short(nm))
        scaff = collections.defaultdict(lambda: [0, 0, 0])
        for (s, ch, h), cts in cc.items():
            if len(cts) > 1: scaff[ch][2] += len(cts) - 1
    # telomeres
    telo = collections.defaultdict(lambda: [0, 0, 0])
    for r in rows:
        if r["type"] != "telomere" or r["chrom"] == "?": continue
        if r["decision"] == "accepted": telo[r["chrom"]][0] += 1
        else:
            rs = r["reason"].lower()
            telo[r["chrom"]][2 if ("no telomere" in rs or "buried" in rs or "no donor reaches" in rs) else 1] += 1
    # gap-free T2T after (per sample): a contig with both telomere caps AND no remaining N-gap
    t2t = collections.defaultdict(int)
    for (sample, chrom, hap), contigs in caps.items():
        if chrom == "?": continue
        if any((l and r) and (sample, ct) not in unfilled for ct, (l, r) in contigs.items()):
            t2t[sample] += 1
    return gaps, scaff, telo, t2t

C = compute("chm13"); S = compute("selfref")

def combined_bar(cdat, sdat, fname, title, xlabel, lab_done):
    chroms = [c for c in sorted(set(cdat) | set(sdat), key=ckey) if c != "?"]
    h = 0.40
    fig, ax = plt.subplots(figsize=(10, max(4.5, len(chroms) * 0.64)))
    yticks = []; yticklabels = []
    for i, c in enumerate(chroms):
        for off, dat, name in [(-0.21, cdat, "CHM13"), (0.21, sdat, "R.F.")]:
            a, b, d = dat.get(c, [0, 0, 0]); tot = a + b + d
            ax.barh(i + off, a, height=h, color=GREEN)
            ax.barh(i + off, b, left=a, height=h, color=ORANGE)
            ax.barh(i + off, d, left=a + b, height=h, color=GREY)
            yticks.append(i + off); yticklabels.append(name)
            if tot:
                ax.text(tot + 0.12, i + off, "%d/%d" % (a, tot) + (" (%d rej)" % b if b else ""),
                        va="center", fontsize=7.5, color="#333")
        ax.text(-0.068, i, c, transform=ax.get_yaxis_transform(), ha="right", va="center",
                fontsize=11, fontweight="bold")
    ax.set_yticks(yticks); ax.set_yticklabels(yticklabels, fontsize=8)
    ax.tick_params(axis="y", length=0); ax.invert_yaxis()
    ax.set_xlabel(xlabel); ax.set_title(title)
    mx = max(max((sum(cdat.get(c, [0, 0, 0])) for c in chroms), default=1),
             max((sum(sdat.get(c, [0, 0, 0])) for c in chroms), default=1))
    ax.set_xlim(0, mx * 1.5); ax.grid(axis="x", ls=":", alpha=0.4)
    ax.legend(handles=[Patch(color=GREEN, label=lab_done), Patch(color=ORANGE, label="rejected by a guard"),
                       Patch(color=GREY, label="none / no donor / unplaced")],
              loc="best", fontsize=9, framealpha=0.9)
    fig.savefig(fname, dpi=140, bbox_inches="tight"); plt.close(fig)

combined_bar(C[0], S[0], "panpatch_1_gaps_combined.png",
             "N-gaps per chromosome — CHM13 vs reference-free", "number of N-gaps", "filled")
combined_bar(C[1], S[1], "panpatch_2_scaffolds_combined.png",
             "Contig fragments per chromosome — CHM13 vs reference-free", "number of contig fragments", "joined")
combined_bar(C[2], S[2], "panpatch_3_telomeres_combined.png",
             "Telomere ends per chromosome — CHM13 vs reference-free", "number of contig ends", "patched")

# ---- combined summary: mechanism rollup (grouped) + T2T after panpatch (CHM13 vs reference-free) ----
def roll(d): return [sum(v[k] for v in d.values()) for k in range(3)]
cg, cs, ct = roll(C[0]), roll(C[1]), roll(C[2]); sg, ss, st = roll(S[0]), roll(S[1]), roll(S[2])
fig, (axL, axR) = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw={"width_ratios": [1.3, 1]})
mech = ["Telomere ends", "Contig fragments", "N-gaps"]
cvals = [ct, cs, cg]; svals = [st, ss, sg]; hh = 0.40
yticks = []; yticklabels = []
for i, m in enumerate(mech):
    for off, v, name in [(-0.21, cvals[i], "CHM13"), (0.21, svals[i], "R.F.")]:
        a, b, d = v; tot = a + b + d
        axL.barh(i + off, a, height=hh, color=GREEN)
        axL.barh(i + off, b, left=a, height=hh, color=ORANGE)
        axL.barh(i + off, d, left=a + b, height=hh, color=GREY)
        yticks.append(i + off); yticklabels.append(name)
        if tot:
            axL.text(tot + 0.7, i + off, "%d done" % a + (", %d rej" % b if b else ""),
                     va="center", fontsize=8.5, color="#333")
    axL.text(-0.07, i, m, transform=axL.get_yaxis_transform(), ha="right", va="center",
             fontsize=10.5, fontweight="bold")
axL.set_yticks(yticks); axL.set_yticklabels(yticklabels, fontsize=8)
axL.tick_params(axis="y", length=0); axL.invert_yaxis()
axL.set_xlabel("count"); axL.set_xlim(0, max(sum(v) for v in cvals + svals) * 1.6)
axL.set_title("What panpatch did (3 samples)")
axL.legend(handles=[Patch(color=GREEN, label="done"), Patch(color=ORANGE, label="rejected by a guard"),
                    Patch(color=GREY, label="none / no donor")], loc="upper right", fontsize=9, framealpha=0.9)
x = np.arange(len(SAMPLES)); w = 0.38
ca = [C[3][s] for s in SAMPLES]; sa = [S[3][s] for s in SAMPLES]
axR.bar(x - w / 2, ca, w, label="CHM13", color=DBLUE)
axR.bar(x + w / 2, sa, w, label="reference-free", color=PURPLE)
for i in range(len(SAMPLES)):
    for xx, v in [(x[i] - w / 2, ca[i]), (x[i] + w / 2, sa[i])]:
        axR.text(xx, v + 0.3, v, ha="center", fontsize=8)
axR.set_xticks(x); axR.set_xticklabels(SAMPLES); axR.set_ylim(0, 50)
axR.set_ylabel("gap-free T2T (of 46)"); axR.set_title("Gap-free T2T after panpatch")
axR.legend(loc="upper left", fontsize=8.5, framealpha=0.9)
fig.suptitle("panpatch summary — gap-free T2T (both telomeres, no N-gaps): CHM13 %d · reference-free %d"
             % (sum(ca), sum(sa)), fontsize=13, y=1.03)
fig.tight_layout(); fig.savefig("panpatch_4_summary_combined.png", dpi=140, bbox_inches="tight"); plt.close(fig)

print("CHM13   gaps=%s scaff=%s telo=%s T2T=%d" % (roll(C[0]), roll(C[1]), roll(C[2]), sum(ca)))
print("selfref gaps=%s scaff=%s telo=%s T2T=%d" % (roll(S[0]), roll(S[1]), roll(S[2]), sum(sa)))
print("wrote panpatch_{1_gaps,2_scaffolds,3_telomeres,4_summary}_combined.png")
