#!/usr/bin/env python3
# Regenerate the panpatch summary plots from the per-haplotype / gap / scaffold CSVs.
import csv
from collections import defaultdict
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

GREEN = "#31a354"; GREY = "#cccccc"; LBLUE = "#9ecae1"; DBLUE = "#3182bd"
def ckey(c):
    c = c[3:]; return {"X": 23, "Y": 24, "M": 25}.get(c, int(c) if c.isdigit() else 99)

results = list(csv.DictReader(open("panpatch_results.csv")))
gaps    = list(csv.DictReader(open("gaps.csv")))
scaf    = list(csv.DictReader(open("scaffolds.csv")))

def chrom_bar(ax, perchr, title, xlabel, done_lbl, not_lbl, annot):
    chroms = sorted(perchr, key=ckey); y = np.arange(len(chroms))
    for i, c in enumerate(chroms):
        d, n = perchr[c]
        ax.barh(i, d, color=GREEN); ax.barh(i, n, left=d, color=GREY)
        ax.text(d + n + 0.15, i, annot(c, d, n), va="center", fontsize=8, color="#333")
    ax.set_yticks(y); ax.set_yticklabels(chroms); ax.invert_yaxis()
    ax.set_xlabel(xlabel); ax.set_title(title)
    mx = max((d + n) for d, n in perchr.values())
    ax.set_xlim(0, mx * 1.45); ax.grid(axis="x", ls=":", alpha=0.4)
    ax.legend(handles=[Patch(color=GREEN, label=done_lbl), Patch(color=GREY, label=not_lbl)],
              loc="lower right", fontsize=9)

# 1) GAPS
g = defaultdict(lambda: [0, 0]); gbp = defaultdict(lambda: [0, 0])
for r in gaps:
    j = 0 if r["status"] == "filled" else 1
    g[r["chrom"]][j] += 1; gbp[r["chrom"]][j] += int(r["gap_len"])
fig, ax = plt.subplots(figsize=(8.5, 5.5))
chrom_bar(ax, g, "N-gaps filled per chromosome", "number of N-gaps",
          "filled by panpatch", "not filled (no donor bridge)",
          lambda c, d, n: "%d/%d  (%.0f kb filled)" % (d, d + n, gbp[c][0] / 1000))
fig.tight_layout(); fig.savefig("panpatch_1_gaps.png", dpi=140); plt.close(fig)

# 2) SCAFFOLDS
s = defaultdict(lambda: [0, 0]); sbp = defaultdict(lambda: [0, 0])
for r in scaf:
    frags = [int(x) for x in r["fragment_bp(desc)"].split(";") if x]
    j = 0 if r["status"] == "joined" else 1
    s[r["chrom"]][j] += len(frags); sbp[r["chrom"]][j] += sum(frags)
fig, ax = plt.subplots(figsize=(8.5, 5))
chrom_bar(ax, s, "Contig fragments scaffolded per chromosome", "number of contig fragments",
          "joined into the chromosome", "not joined (reverted separate)",
          lambda c, d, n: "%d/%d  (%.0f kb joined)" % (d, d + n, sbp[c][0] / 1000))
fig.tight_layout(); fig.savefig("panpatch_2_scaffolds.png", dpi=140); plt.close(fig)

# 3) TELOMERES
frag_haps = {(r["sample"], r["chrom"], r["hap"]) for r in scaf if r["status"] == "fragmented"}
t = defaultdict(lambda: [0, 0]); tbp = defaultdict(int)
for r in results:
    key = (r["sample"], r["chrom"], r["hap"])
    if r["outcome"] == "telomere_patch":
        t[r["chrom"]][0] += 1; tbp[r["chrom"]] += int(r["grafted_bp"])
    elif r["outcome"] == "telomere_attempt_reverted" or (r["outcome"] == "incomplete" and key not in frag_haps):
        t[r["chrom"]][1] += 1
fig, ax = plt.subplots(figsize=(8.5, 4.8))
chrom_bar(ax, t, "Telomere ends patched per chromosome", "number of contig ends",
          "telomere grafted from donor", "not patched (no donor / degraded)",
          lambda c, d, n: "%d/%d%s" % (d, d + n, ("  (+%.0f kb)" % (tbp[c] / 1000) if tbp[c] else "")))
fig.tight_layout(); fig.savefig("panpatch_3_telomeres.png", dpi=140); plt.close(fig)

# 4) SUMMARY
samples = ["HG00126", "HG00235", "HG01074"]
gf = sum(v[0] for v in g.values()); gn = sum(v[1] for v in g.values()); gfbp = sum(v[0] for v in gbp.values())
sf = sum(v[0] for v in s.values()); sn = sum(v[1] for v in s.values()); sfbp = sum(v[0] for v in sbp.values())
tf = sum(v[0] for v in t.values()); tn = sum(v[1] for v in t.values()); tfbp = sum(tbp.values())
fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={"width_ratios": [1.1, 1]})
mech = ["Telomere\ncompletions", "Contig\nscaffolds", "N-gap\nfills"]
done = [tf, sf, gf]; nots = [tn, sn, gn]
bpnote = ["+%.0f kb grafted" % (tfbp / 1000), "+%.2f Mb joined" % (sfbp / 1e6), "+%.2f Mb filled" % (gfbp / 1e6)]
yy = np.arange(len(mech))
axL.barh(yy, done, color=GREEN); axL.barh(yy, nots, left=done, color=GREY)
for i in range(len(mech)):
    axL.text(done[i] + nots[i] + 0.6, i, "%d of %d done   %s" % (done[i], done[i] + nots[i], bpnote[i]),
             va="center", fontsize=9, color="#222")
axL.set_yticks(yy); axL.set_yticklabels(mech); axL.invert_yaxis()
axL.set_xlabel("number of events"); axL.set_xlim(0, max(d + n for d, n in zip(done, nots)) * 1.7)
axL.set_title("What panpatch did (3 samples combined)")
axL.legend(handles=[Patch(color=GREEN, label="patched"), Patch(color=GREY, label="not patched")], loc="lower right", fontsize=9)
x = np.arange(len(samples)); w = 0.38
bef = [sum(int(r["t2t_before"]) for r in results if r["sample"] == sm) for sm in samples]
aft = [sum(int(r["t2t_after"])  for r in results if r["sample"] == sm) for sm in samples]
axR.bar(x - w / 2, bef, w, label="before", color=LBLUE); axR.bar(x + w / 2, aft, w, label="after", color=DBLUE)
for i, (b, a) in enumerate(zip(bef, aft)):
    axR.text(i - w / 2, b + 0.3, b, ha="center", fontsize=9); axR.text(i + w / 2, a + 0.3, a, ha="center", fontsize=9)
    axR.text(i + w / 2, a + 2.0, "+%d" % (a - b), ha="center", fontsize=9, color="#08519c", fontweight="bold")
axR.set_xticks(x); axR.set_xticklabels(samples); axR.set_ylim(0, 50)
axR.set_ylabel("single-contig T2T (of 46)"); axR.set_title("Outcome: T2T haplotype-chromosomes")
axR.legend(loc="lower right")
fig.suptitle("panpatch summary — T2T %d → %d (+%d) across 3 samples" % (sum(bef), sum(aft), sum(aft) - sum(bef)), fontsize=13, y=1.02)
fig.tight_layout(); fig.savefig("panpatch_4_summary.png", dpi=140, bbox_inches="tight"); plt.close(fig)
print("wrote panpatch_1_gaps.png panpatch_2_scaffolds.png panpatch_3_telomeres.png panpatch_4_summary.png")
