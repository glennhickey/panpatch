# HPRC R3 panpatch test (HG00126, HG00235, HG01074)

Scripts to reproduce the panpatch evaluation on three HPRC release-3 samples, run two ways:

- **CHM13-referenced** — one graph per sample (`--reference CHM13`); target = `verkko-R3-<id>`, donor = `hifiasm-R2-<id>`. Supports scaffolding and gap-fills across the whole genome.
- **Self-referenced (no CHM13)** — two graphs per sample, one per verkko haplotype, each using *that haplotype as its own reference* (`--reference verkko-R3-<id>_{1,2}`), donor = `hifiasm-R2-<id>`. Better satellite/centromere alignment (same individual), but no cross-reference scaffolding. Here reference == target (both still passed to panpatch).

## Pipeline

All paths assume a scratch work dir (we used `~/dev/work/panpatch-jun24`); adjust `ROOT`/`PATCHOUT` in the run scripts to taste. `panpatch` is `~/dev/panpatch/panpatch`.

1. **Build graphs** (cluster, slurm; needs the per-sample `*.seqfile`s):
   ```
   ./make-graphs.sh            # CHM13 graphs       -> mc-<sample>/
   ./make-graphs-no-chm13.sh   # self-ref graphs    -> mc-<sample>.no-chm13.{1,2}/
   ```
2. **Collect per-chromosome graphs** from cactus's `--chrom-vg full` output into:
   - `<sample>.chroms/chr*.full.vg`                       (CHM13; named by chromosome)
   - `<sample>.no-chm13.{1,2}.chroms/haplotype{1,2}-*.full.vg`  (self-ref; named by verkko contig; `chrOther` is skipped)
3. **Patch:**
   ```
   PATCHOUT=patch-runs ./run-patch.sh                     # CHM13, all three samples (parallel)
   ./run-nochm13.sh HG01074.no-chm13.1.chroms verkko-R3-HG01074_1 hifiasm-R2-HG01074
   ./run-nochm13.sh HG01074.no-chm13.2.chroms verkko-R3-HG01074_2 hifiasm-R2-HG01074
   ```
   Both write one `<run>.bed` (panpatch stdout) + `<run>.stderr` per chromosome and are resumable.
4. **Analyze** (run from the dir holding `patch-runs/`; override with `PATCHRUNS=<dir>`):
   ```
   python3 make_results_csv.py     # -> panpatch_results.csv, panpatch_summary.csv  (per hap-chromosome T2T + outcome)
   python3 make_scaffolds_csv.py   # -> scaffolds.csv          (multi-contig haps: joined vs fragmented)
   python3 make_gaps_csv.py        # -> gaps.csv               (every N-gap: filled vs retained; streams verkko, slow)
   python3 make_plots.py           # -> panpatch_{1_gaps,2_scaffolds,3_telomeres,4_summary}.png
   python3 make_patch_table.py     # -> patches.csv, patches.md (one row per kept patch)
   ```

## Results (for the record)

- **CHM13, 3 samples:** T2T (single contig + both telomeres) **115 → 125 (+10)** — 10 telomere completions, 3 N-gap fills, **0 scaffolds**. The only scaffold candidate (HG01074 chr14) was a CHM13-subtelomere artifact, rejected by the telomere-preservation guard; the main contig alone was already T2T.
- **Self-reference (HG01074 only so far):** all 46 hap-chromosomes come out single-contig with both telomeres. Cleaner than CHM13 in a few small cases — completes chr9·h2 that CHM13 cannot (CHM13 reverts a 335 kb deletion misjoin), and fills chr12 more faithfully (1.2 Mb @ 98% vs CHM13's 2.7 Mb @ 73%). rDNA misjoins still occur (chr22, both haps) but are caught by the k-mer guard. Memory: a 2 GB satellite graph peaks ~4.9 GB.

Patches that survive all guards are listed in `patches.md`; the guards themselves and the messages they emit are documented in the top-level `README.md` (*Why a patch is rejected*).
