#!/usr/bin/env python3
# Generate the 4 panpatch slide figures (3-way: done / rejected / no-donor-or-bridge) for a mode,
# from the new TSV reports in deliverable.  Usage: python3 make_allfigs.py {chm13|selfref}
import sys, glob, os, re, csv, collections
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

GREEN = "#31a354"; ORANGE = "#e6550d"; GREY = "#cccccc"; LBLUE = "#9ecae1"; DBLUE = "#3182bd"
def ckey(c):
    c = c[3:]; return {"X": 23, "Y": 24, "M": 25}.get(c, int(c) if c.isdigit() else 99)
def short(c):
    m = re.search(r"(haplotype\d-\d+)", c); return m.group(1) if m else c
def hapnum(c):
    m = re.search(r"haplotype(\d)", c); return int(m.group(1)) if m else 0

mode = sys.argv[1] if len(sys.argv) > 1 else "chm13"
TAG = "" if mode == "chm13" else "_selfref"
LABEL = "CHM13-referenced" if mode == "chm13" else "reference-free (self)"
SAMPLES = ["HG00126", "HG00235", "HG01074"]

# verkko contig -> chromosome, from the CHM13 BED
c2chrom = {}
for bed in glob.glob("deliverable/chm13/*.bed"):
    sample = os.path.basename(bed).split(".")[0]; chrom = None
    for L in open(bed):
        m = re.match(r"#Patched assembly on (chr\S+) for", L)
        if m: chrom = m.group(1); continue
        m = re.match(r"(verkko\S+)\t", L)
        if m and chrom: c2chrom[(sample, short(m.group(1)))] = chrom

def report_paths():
    if mode == "chm13":
        return [(os.path.basename(f).split(".")[0], f) for f in glob.glob("deliverable/chm13/*.report")]
    return [(os.path.basename(os.path.dirname(f)).split(".hap")[0], f) for f in glob.glob("deliverable/self-ref/*/report")]

# parse rows + verkko contig caps, mapped to (sample, real-chrom, hap)
rows = []; caps = collections.defaultdict(dict)
for sample, f in report_paths():
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
            rows.append(dict(sample=sample, chrom=ch, hap=hapnum(c[3]), type=c[2],
                             decision=c[11], reason=c[12]))

# ---- gaps: reuse the verkko N-gap list from gaps.csv, classify for this mode ----
kept = collections.defaultdict(list); rej = collections.defaultdict(list)
beds = glob.glob("deliverable/chm13/*.bed") if mode == "chm13" else glob.glob("deliverable/self-ref/*/out.bed")
for bed in beds:
    sample = os.path.basename(bed).split(".")[0] if mode == "chm13" else os.path.basename(os.path.dirname(bed)).split(".hap")[0]
    for L in open(bed):
        m = re.match(r"(verkko\S+)\t(\d+)\t(\d+)\t", L)
        if m: kept[(sample, short(m.group(1)))].append((int(m.group(2)), int(m.group(3))))
for sample, f in report_paths():
    for L in open(f):
        if L.startswith("#") or L.startswith("chrom\t"): continue
        c = L.rstrip("\n").split("\t")
        if len(c) >= 15 and c[2] == "gap-fill" and c[11] == "rejected" and c[13] != "." and c[14] != ".":
            rej[(sample, short(c[3]))].append((int(c[13]), int(c[14])))
def ov(gs, gl, rngs): return any(gs < b and a < gs + gl for a, b in rngs)
gaps = collections.defaultdict(lambda: [0, 0, 0]); unfilled = set()   # chrom -> [filled, rejected, no_donor]
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

# ---- scaffolds: fragments per chrom (joined / rejected / no-bridge) ----
scaff_dec = {}
for r in rows:
    if r["type"] == "scaffold" and r["chrom"] != "?":
        k = (r["sample"], r["chrom"], r["hap"])
        if scaff_dec.get(k) != "rejected": scaff_dec[k] = "joined" if r["decision"] == "accepted" else "rejected"
scaff = collections.defaultdict(lambda: [0, 0, 0])  # chrom -> [joined, rejected, no_bridge] (fragment counts)
for (sample, chrom, hap), contigs in caps.items():
    if chrom == "?" or len(contigs) < 2: continue
    frags = len(contigs) - 1; d = scaff_dec.get((sample, chrom, hap))
    scaff[chrom][0 if d == "joined" else 1 if d == "rejected" else 2] += frags
if mode != "chm13":
    # self-ref uses per-contig graphs -> no scaffolding; show the shared verkko fragment structure
    # (from the CHM13 #Contig grouping) as "left unplaced (chrOther)".
    cc = collections.defaultdict(set)
    for s in SAMPLES:
        for L in open("deliverable/chm13/%s.report" % s):
            if L.startswith("#Contig") and "verkko" in L:
                nm = L.split()[1]; ch = c2chrom.get((s, short(nm)))
                if ch: cc[(s, ch, hapnum(nm))].add(short(nm))
    scaff = collections.defaultdict(lambda: [0, 0, 0])
    for (s, ch, h), cts in cc.items():
        if len(cts) > 1: scaff[ch][2] += len(cts) - 1

# ---- telomeres: ends per chrom (patched / rejected / no-donor) ----
telo = collections.defaultdict(lambda: [0, 0, 0])
for r in rows:
    if r["type"] != "telomere" or r["chrom"] == "?": continue
    if r["decision"] == "accepted": telo[r["chrom"]][0] += 1
    else:
        rs = r["reason"].lower()
        telo[r["chrom"]][2 if ("no telomere" in rs or "buried" in rs or "no donor reaches" in rs) else 1] += 1

# ---- gap-free T2T after per (sample,chrom,hap): a contig with both telomere caps AND no remaining N-gap ----
t2t_after = collections.defaultdict(bool)
for (sample, chrom, hap), contigs in caps.items():
    if chrom == "?": continue
    t2t_after[(sample, chrom, hap)] = any((l and r) and (sample, ct) not in unfilled for ct, (l, r) in contigs.items())

# ================= plotting =================
def bar3(ax, data, title, xlabel, labels, annot):
    chroms = sorted(data, key=ckey); y = np.arange(len(chroms))
    for i, c in enumerate(chroms):
        a, b, d = data[c]
        ax.barh(i, a, color=GREEN); ax.barh(i, b, left=a, color=ORANGE); ax.barh(i, d, left=a + b, color=GREY)
        ax.text(a + b + d + 0.15, i, annot(a, b, d), va="center", fontsize=8, color="#333")
    ax.set_yticks(y); ax.set_yticklabels(chroms); ax.invert_yaxis()
    ax.set_xlabel(xlabel); ax.set_title(title)
    mx = max((sum(v) for v in data.values()), default=1); ax.set_xlim(0, mx * 1.5); ax.grid(axis="x", ls=":", alpha=0.4)
    ax.legend(handles=[Patch(color=GREEN, label=labels[0]), Patch(color=ORANGE, label=labels[1]),
                       Patch(color=GREY, label=labels[2])], loc="best", fontsize=9, framealpha=0.9)

# 1) gaps
fig, ax = plt.subplots(figsize=(8.5, 5.5))
bar3(ax, gaps, "N-gaps per chromosome (%s)" % LABEL, "number of N-gaps",
     ["filled by panpatch", "donor graft rejected by a guard", "no donor reached it"],
     lambda a, b, d: "%d/%d filled" % (a, a + b + d) + (("  (%d rejected)" % b) if b else ""))
fig.tight_layout(); fig.savefig("panpatch_1_gaps%s.png" % TAG, dpi=140); plt.close(fig)

# 2) scaffolds
fig, ax = plt.subplots(figsize=(8.5, 5))
if mode == "chm13":
    bar3(ax, scaff, "Contig fragments scaffolded per chromosome (%s)" % LABEL, "number of contig fragments",
         ["joined into the chromosome", "join rejected by a guard", "not joined (no bridge)"],
         lambda a, b, d: "%d/%d joined" % (a, a + b + d) + (("  (%d rejected)" % b) if b else ""))
else:
    bar3(ax, scaff, "Contig fragments per chromosome (%s) — none scaffolded" % LABEL, "number of contig fragments",
         ["joined into the chromosome", "join rejected by a guard", "left unplaced (self-ref does not scaffold)"],
         lambda a, b, d: "%d unplaced" % d)
    ax.text(0.99, 0.06, "reference-free mode patches each contig in its own graph;\nfragments stay as unplaced contigs (chrOther)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8, style="italic", color="#555")
fig.tight_layout(); fig.savefig("panpatch_2_scaffolds%s.png" % TAG, dpi=140); plt.close(fig)

# 3) telomeres
fig, ax = plt.subplots(figsize=(8.5, 4.8))
bar3(ax, telo, "Telomere ends patched per chromosome (%s)" % LABEL, "number of contig ends",
     ["telomere grafted from donor", "patch rejected by a guard", "no donor / degraded"],
     lambda a, b, d: "%d/%d patched" % (a, a + b + d) + (("  (%d rejected)" % b) if b else ""))
fig.tight_layout(); fig.savefig("panpatch_3_telomeres%s.png" % TAG, dpi=140); plt.close(fig)

# 4) summary
g = [sum(v[k] for v in gaps.values()) for k in range(3)]
s = [sum(v[k] for v in scaff.values()) for k in range(3)]
t = [sum(v[k] for v in telo.values()) for k in range(3)]
fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={"width_ratios": [1.2, 1]})
mech = ["Telomere\nends", "Contig\nfragments", "N-gaps"]
done = [t[0], s[0], g[0]]; rejc = [t[1], s[1], g[1]]; none = [t[2], s[2], g[2]]
yy = np.arange(len(mech))
axL.barh(yy, done, color=GREEN); axL.barh(yy, rejc, left=done, color=ORANGE)
axL.barh(yy, none, left=[done[i] + rejc[i] for i in range(3)], color=GREY)
for i in range(3):
    axL.text(done[i] + rejc[i] + none[i] + 0.6, i, "%d done, %d rejected, %d none" % (done[i], rejc[i], none[i]),
             va="center", fontsize=9, color="#222")
axL.set_yticks(yy); axL.set_yticklabels(mech); axL.invert_yaxis()
axL.set_xlabel("count"); axL.set_xlim(0, max(done[i] + rejc[i] + none[i] for i in range(3)) * 1.7)
axL.set_title("What panpatch did (%s, 3 samples)" % LABEL)
axL.legend(handles=[Patch(color=GREEN, label="done"), Patch(color=ORANGE, label="rejected by a guard"),
                    Patch(color=GREY, label="none / no donor")], loc="upper right", fontsize=9, framealpha=0.9)
x = np.arange(len(SAMPLES)); w = 0.5
aft = [sum(t2t_after.get((sm, c, h), False) for (s2, c, h) in t2t_after if s2 == sm) for sm in SAMPLES]
axR.bar(x, aft, w, label="after panpatch", color=DBLUE)
for i, a in enumerate(aft):
    axR.text(i, a + 0.3, a, ha="center", fontsize=9)
axR.set_xticks(x); axR.set_xticklabels(SAMPLES); axR.set_ylim(0, 50)
axR.set_ylabel("gap-free T2T (of 46)"); axR.set_title("Gap-free T2T after panpatch (%s)" % LABEL)
fig.suptitle("panpatch summary (%s) — gap-free T2T (both telomeres, no N-gaps): %d" % (LABEL, sum(aft)), fontsize=13, y=1.02)
fig.tight_layout(); fig.savefig("panpatch_4_summary%s.png" % TAG, dpi=140, bbox_inches="tight"); plt.close(fig)

print("[%s] gaps f/r/n=%s  scaff j/r/n=%s  telo p/r/n=%s  T2T after=%d" %
      (mode, g, s, t, sum(aft)))
print("wrote panpatch_{1_gaps,2_scaffolds,3_telomeres,4_summary}%s.png" % TAG)
