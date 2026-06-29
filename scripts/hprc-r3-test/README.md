# HPRC R3 panpatch test (HG00126, HG00235, HG01074)

Scripts to reproduce the panpatch evaluation on three HPRC release-3 samples, run two ways:

- **CHM13-referenced** — one graph per sample (`--reference CHM13`); target = `verkko-R3-<id>`, donor = `hifiasm-R2-<id>`. Supports scaffolding and gap-fills across the whole genome.
- **Self-referenced (no CHM13)** — two graphs per sample, one per verkko haplotype, each using *that haplotype as its own reference* (`--reference verkko-R3-<id>_{1,2}`), donor = `hifiasm-R2-<id>`. Better satellite/centromere alignment (same individual), but per-contig graphs so no cross-reference scaffolding (reference == target, both passed to panpatch).

`panpatch` is `~/dev/panpatch/panpatch` (override with `PANPATCH=...`).

## Pipeline

1. **Build graphs** (cluster, slurm; needs the per-sample `*.seqfile`s):
   ```
   ./make-graphs.sh            # CHM13 graphs    -> mc-<sample>/
   ./make-graphs-no-chm13.sh   # self-ref graphs -> mc-<sample>.no-chm13.{1,2}/
   ```
2. **Collect per-chromosome graphs** from cactus's `--chrom-vg full` output into:
   - `<sample>.chroms/chr*.full.vg`                              (CHM13; one graph per chromosome)
   - `<sample>.no-chm13.{1,2}.chroms/*.full.vg`                  (self-ref; one graph per verkko contig, plus `chrOther`)
3. **Patch.** Each run makes a single multi-graph panpatch call and writes the deliverable layout directly:
   ```
   ./run-patch.sh                                                  # CHM13, all three samples
   ./run-nochm13.sh HG01074.no-chm13.1.chroms verkko-R3-HG01074_1 hifiasm-R2-HG01074
   ./run-nochm13.sh HG01074.no-chm13.2.chroms verkko-R3-HG01074_2 hifiasm-R2-HG01074
   ```
   Outputs (under `deliverable/`, override with `PATCHOUT=` / the 4th arg):
   - `deliverable/chm13/<id>.{report,bed,hap1.fa,hap2.fa,stderr}`
   - `deliverable/self-ref/<id>.hap<h>/{report,out.bed,out.hap0.fa,stderr}`

   where **`report`** is the TSV report (panpatch stdout), **`.bed`/`out.bed`** the patched-assembly intervals (`--bed`), and the **FASTAs** the per-haplotype patched sequence (`-f`, diploid split automatic; self-ref's single haplotype is `out.hap0.fa`). Both runners are resumable (skip a target whose report exists) and print a per-target accept/reject tally read straight from the report's `type`/`decision` columns.
4. **Analyze** (run from the dir holding `deliverable/` and the `<sample>.chroms/` graphs):
   ```
   python3 make_gaps_csv2.py            # -> gaps.csv   (every verkko N-gap: filled/rejected/no_donor; scans the graphs, slow)
   python3 make_allfigs.py chm13        # -> panpatch_{1_gaps,2_scaffolds,3_telomeres,4_summary}.png
   python3 make_allfigs.py selfref      # -> the same four, *_selfref.png
   python3 make_combined.py             # -> *_combined.png  (CHM13 vs reference-free, two bars per chromosome)
   python3 make_slides2.py              # -> chm13_patches.md + diff_chm13_vs_selfref.txt
   ```
   The figures are 3-way (done / rejected-by-a-guard / no-donor); `make_gaps_csv2.py` writes `gaps.csv`, which the figure scripts consume.

## Output model

panpatch writes three independent sinks: the **report TSV to stdout** (one row per candidate patch — see the top-level `README.md`), the **patched intervals to `--bed FILE`**, and the **per-haplotype FASTA to `-f FILE`** (`FILE.hap<N>.fa`). The BED/FASTA are written only on full success. There is no separate assembly step — `-f` already emits one complete FASTA per haplotype (patched + reverted/passthrough contigs).

## Caveats

- **Completeness.** These FASTAs contain what survives into the per-chromosome graphs, *not* necessarily the full original verkko assembly. cactus-pangenome can drop contigs during graph construction, and the CHM13 build has no `chrOther`, so unplaced verkko contigs are lost there (e.g. HG01074 hap1 loses `haplotype1-0000031`/`-0000059`, ~0.5 + 0.7 Mb). The self-reference build keeps a `chrOther` (panpatch passes it through) and recovers those — but contigs dropped at construction are invisible to both. Pull missing contigs from the original verkko FASTAs (URLs in the `*.seqfile`s) if a strictly complete assembly is required.
- **chrM.** verkko has no mitochondrion in these graphs (chrM is CHM13-only), so the deliverables carry none.
- **Headers.** Patched records are `chrN_hap_X` (CHM13) or `<contig>_hap_0` (self-reference); unpatched/reverted/passthrough records keep their original verkko contig names. The `.bed` files give full provenance, and the `report` gives the per-patch decision and reason.

## Results (for the record)

- **CHM13, 3 samples:** 11 telomere completions + 3 N-gap fills, **0 scaffolds** (all 23 candidate fragment-joins rejected as repeat-region misjoins). **Gap-free T2T (single contig, both telomeres, no N-gaps): 91.**
- **Self-reference, 3 samples:** fills *more* N-gaps than CHM13 (10 vs 3 — hifiasm aligns directly to the verkko target rather than through CHM13) but also *rejects more* (12 vs 2, incl. 4 divergent-flank fills caught by the flank-anchoring guard), and does no scaffolding (per-contig graphs; fragments stay unplaced in `chrOther`). **Gap-free T2T: 95 — reference-free wins**, because for completeness interior gap-filling outweighs terminal-telomere completion. (On the looser both-telomeres-only metric the two are ~tied, 126 vs 125; 50 of the 60 N-gaps are acrocentric rDNA with no clean bridge in either mode.)

Patches that survive all guards are listed in `chm13_patches.md`; the guards and the messages they emit are documented in the top-level `README.md` (*Why a patch is rejected*).
