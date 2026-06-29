# panpatch

Use a pangenome graph to patch (slightly) fragmented assemblies into T2T chromosomes.

Panpatch supports two types of patches:

1) [Gaps inside scaffolds (`N`s)](#patching-gaps)
2) [Scaffolding disconnected contigs](#scaffolding)

In all cases, `panpatch` requires a "reference" assembly with chromosome-scale scaffolds that it uses for orientation (`-r`), an assembly to patch (first `-s`), and one or more assemblies to use for patching (subsequent `-s`).  If the assembly being patched is chromosome-scale, it can be used as the reference. 

You must always begin by building a pangenome graph of your assemblies with [Minigraph-Cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md).  Minigraph-Cactus requires you specify a reference backbone with `--reference` and this should be used consistently with the `-r` reference option in `panpatch`.  It is important to remember that Minigraph-Cactus only supports haploid references, so in some cases this will require building two separate graphs.  

## Patching Gaps

In this scenario, let's say you have a diploid assembly, `PAN028-verkko`, with chromosome scale scaffolds. But it contains gaps (represented as runs of `N`s) that you want to patch with two other assemblies, `PAN028-hifiasm` and `PAN028-duplex`.  These assemblies are also diploid, but their haplotype annotations aren't necessarily the same (ie maybe `chr2#1` from duplex is the same as `chr2#2` from verkko).

Because the reference assembly is diploid, you need to make two graphs and patch the two haplotypes separately.  Below is an example for haplotype 1 (you'd have to repeat the process for haplotype 2).

Begin by making the graph. You specify the inputs to Cactus with a two-column file  `pan028.verkko1.seqfile`:

```
PAN028-verkko_1  PAN028.haplotype1.full.verkko2.fa
PAN028-hifiasm.1 PAN028.hifiasm.20240417.hic.hap1.fa
PAN028-hifiasm.2 PAN028.hifiasm.20240417.hic.hap2.fa
PAN028-duplex.1  PAN028.haplotype1.duplex.verkko2.0.fa
PAN028-duplex.2  PAN028.haplotype2.duplex.verkko2.0.scaff.fa
```

The `.1/.2` suffixes specify the haplotype, but because diploid references aren't supported we use an underscore.  More information [here](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md#sample-names).

Now make the graph.  The `--reference` and `--chrom-vg full` options are essential. You may need to adapt the others for your computing environment. 
```
cactus-pangenome ./js ./pan028.verkko1.seqfile --outName pan028-mc-verkko-1 --outDir pan028-mc-verkko-1 --logFile pan028-mc-verkko1.log --reference PAN028-verkko_1 --chrom-vg full --batchSystem slurm --consCores 65 --mgCores 64 --indexCores 64 --mapCores 16
```

Now run `panpatch` over all the chromosome graphs in `pan028-mc-verkko-1/pan028-mc-verkko-1.chroms/` in a single command:
```
cd pan028-mc-verkko-1/pan028-mc-verkko-1.chroms/
panpatch *.full.vg -r PAN028-verkko_1 -p -s PAN028-verkko_1 -s PAN028-hifiasm -s PAN028-duplex \
    --bed patched.bed -f patched.fa > patched.report 2> patched.stderr
```

This produces `patched.report` (the per-patch report table, on stdout), `patched.bed` (the patched intervals), and a per-haplotype FASTA (`patched.hap*.fa`).  A chromosome that couldn't be patched contributes its input contig(s) unchanged.  With `-p`, progress and the reasons anything wasn't patched go to stderr; the report table on stdout records the accept/reject decision for every candidate patch.

Note that even though the `hifiasm` and `duplex` assemblies are diploid, only the most relevant haplotype for each will be selected for each chromosome.

When running on the second haplotype, the only difference is the first line of the seqfile, as well as the `--reference / -r` options.

## Scaffolding

If the assembly you want to patch does not have chromosome-scale scaffolds, you must use a reference that does.  Here is an example of using T2T-CHM13 (aka `hs1`) as the reference to patch `PAN028-verkko`.  In this case, since we only have one reference we can do all the patching at once, with a single graph.  For example, use `pan028.hs1.seqfile`:

```
hs1              https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz
PAN028-verkko.1  PAN028.haplotype1.full.verkko2.fa
PAN028-verkko.2  PAN028.haplotype2.full.verkko2.fa
PAN028-hifiasm.1 PAN028.hifiasm.20240417.hic.hap1.fa
PAN028-hifiasm.2 PAN028.hifiasm.20240417.hic.hap2.fa
PAN028-duplex.1  PAN028.haplotype1.duplex.verkko2.0.fa
PAN028-duplex.2  PAN028.haplotype2.duplex.verkko2.0.scaff.fa
```

(note the `.1/.2` suffixes in the first column denote haplotype and are [important](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md#sample-names))

The alignment is then done with

```
cactus-pangenome ./js ./pan028.hs1.seqfile --outName pan028-mc-hs1 --outDir pan028-mc-hs1 --logFile pan028-mc-hs1.log --reference hs1 --consCores 65 --batchSystem slurm --mgCores 64 --indexCores 64 --mapCores 16 --chrom-vg full
```

Now run `panpatch` over all the chromosome graphs in `pan028-mc-hs1/pan028-mc-hs1.chroms/` in a single command:
```
cd pan028-mc-hs1/pan028-mc-hs1.chroms/
panpatch *.full.vg -r hs1 -p -s PAN028-verkko -s PAN028-hifiasm -s PAN028-duplex \
    --bed patched.bed -f patched.fa > patched.report 2> patched.stderr
```

Unlike the previous example the output here is diploid, since both verkko haplotypes are patched, so you get two FASTAs: `patched.hap1.fa` and `patched.hap2.fa` (plus the single `patched.bed` and `patched.report`).

## Building Panpatch

Clone it with submodules then `make`.  The `panpatch` binary should be built in the same directory if all went well.

```
git clone --recursive https://github.com/glennhickey/panpatch.git --branch development
cd panpatch
make
```

If you are missing system dependencies, you can try following the [Install Dependencies](https://github.com/vgteam/vg?tab=readme-ov-file#linux-install-dependencies) section from vg -- `panpatch` uses a subset of these so they will be more than sufficient. 

Linux binaries are available [here](https://github.com/glennhickey/panpatch/releases).

## PanPatch Interface

You specify the graph(s) and sample names in order of priority (first column of the above input file, excluding `.1/2` suffixes)
```
panpatch <graph.vg> [graph2.vg ...] -r <reference sample> -s <sample to patch> -s <first sample to patch with> -s <second sample> etc.
```

For example

```
panpatch chr20.full.vg -r hs1 -s PAN028-verkko -s PAN028-hifiasm -s PAN028-duplex
```

will patch the `PAN028-verkko` assembly, using `PAN028-hifiasm` where possible, then `PAN028-duplex` as a backup.

**Multiple graphs.**  You can pass more than one graph (e.g. `panpatch chr*.full.vg ...`) to patch a whole assembly in one run.  The inputs are processed in lexicographic order and their results concatenated.  Before any patching, panpatch scans every input to confirm the `-r`/`-s` samples exist (a typo fails immediately); a graph that legitimately lacks the target sample (e.g. `chrY` for a female) is skipped with a note.  A graph with no single reference path — e.g. a `chrOther` graph holding many unplaced contigs rather than one chromosome — is **passed through**: its target contigs are emitted unchanged (as `passthrough` rows), so every input contig still reaches the output.

**Three outputs.**

- **the report** — a TSV table, one row per candidate patch, streamed to **stdout** as each contig is processed.  Columns: `chrom hap type target target_bp donor donor_bp replaced_bp kmer% flankL% flankR% decision reason target_start target_end`.  `type` is `telomere`, `gap-fill`, or `scaffold`; a 2-sided graft shows both its k-mer recovery (content) and flank anchoring (locus), a telomere shows its k-mer only (`.` = n/a); `decision` is `accepted` or `rejected`, with the `reason` column giving the rejection cause; `target_start`/`target_end` give the target-contig region (forward coords) a graft replaced (`.` for telomere, scaffold, and passthrough rows).  Per-contig telomere cap status (`#Contig ...`) is printed after the rows.  These numbers let you scrutinise borderline calls and retune the guards (`--graft-recovery`, `--min-flank`, ...) for your own assemblies.
- **`--bed FILE`** writes the patched-assembly contig intervals (BED format, spanning each reference chromosome telomere-to-telomere).
- **`-f/--fasta FILE`** writes the patched sequence as **one FASTA per haplotype** (`FILE.hap1.fa`, `FILE.hap2.fa`, ...), giving the diploid split automatically.

The BED and FASTA are written only **after every input graph has been processed successfully**, so an error never leaves a partial BED/FASTA behind.  Use `--exclude-bed FILE` to prevent patching in specific regions (see [Excluding Regions from Patching](#excluding-regions-from-patching)).

### Options

| option | description | default |
|---|---|---|
| `-r, --reference STR` | reference sample (required) | |
| `-s, --sample STR` | sample to patch (first), then donors in priority order (required, repeatable) | |
| `-f, --fasta FILE` | write the patched assembly to FASTA, one file per haplotype (`FILE.hap1.fa`, ...); written only on full success | |
| `--bed FILE` | write the patched-assembly intervals (BED) to FILE; written only on full success | |
| `-e, --default-sample STR` | use this sample's contig when a patch is rejected | |
| `-b, --exclude-bed FILE` | target regions to leave untouched (input filter; distinct from `--bed`) | |
| `-w, --window N` | window size for haplotype-identity binning | 1000 |
| `-t, --threads N` | threads | all |
| `-T, --require-telomeres` | require a telomere at both ends and none internal | off |
| `-M, --max-telomere-patch N` | max bp a telomere graft may replace | 500000 |
| `--telomere-threshold F` | min hexamer density to call a telomere | 0.8 |
| `--min-cover F` | revert a patch covering less than this fraction of the input length | 0.95 |
| `--graft-recovery F` | revert a foreign interior graft sharing less than this % of the replaced k-mers | 50 |
| `--graft-min-bp N` | apply `--graft-recovery` only when at least this many non-N bp are replaced | 10000 |
| `--min-flank F` | revert an N-gap fill anchored to less than this % of the target flank | 50 |
| `--flank-window N` | window (bp) each side of an N-gap fill over which flank anchoring is measured | 500000 |
| `-p, --progress` | print progress to stderr | |

The report table on **stdout** has one row per candidate patch, e.g.:

```
chrom  hap  type      target              target_bp  donor       donor_bp   replaced_bp  kmer%  flankL%  flankR%  decision  reason                                               target_start  target_end
chr12  2    gap-fill  haplotype2-0000064  132285855  CM088792.1  133100000  1175239      98.4   100.0    100.0    accepted  .                                                    63000000      64175239
chr14  1    gap-fill  haplotype1-0000004  101799395  CM090131.1  101948476  197297       0.4    98.1     100.0    rejected  k-mer recovery 0.4% < 50.0% (repeat-region misjoin)  101051458     101248755
```

The optional **`--bed`** file lists the contig intervals of the patched assembly (the path taken through the graph for each haplotype):

```
#Patched assembly on chr8 for PAN028-verkko#1:
PAN028-hifiasm#2#h2tg000032l#0	2	27	+
PAN028-verkko#1#haplotype1-0000008#0	0	64821530	+
PAN028-verkko#1#haplotype1-0000043#0	11	673283	+
```

(Small intervals should probably be filtered out; there's no such logic yet in `panpatch`.)

## Excluding Regions from Patching

Use the `-b/--exclude-bed` option to provide a BED file of regions that panpatch should not attempt to patch.  Coordinates are in the target assembly being patched (the first `-s` sample).  In excluded regions, panpatch will always keep the original target assembly sequence and never substitute sequence from another assembly.

This is primarily useful for acrocentric chromosomes, where rDNA-adjacent gaps can get overfilled by patches, or where a centromere gap patch should be skipped.  It can also be used to prevent patching of any region where you want to preserve the original assembly sequence.

For example, to prevent panpatch from patching the centromere region of chr22 (leaving the rDNA gap patch intact):

```
echo -e "PAN027-verkko-1#0#PAN027.chr22.paternal\t11000000\t13000000" > exclude.bed
panpatch chr22.full.vg -r PAN027-verkko-1 -s PAN027-verkko-1 -s PAN027-herro -s PAN011-verkko -b exclude.bed
```

The BED file uses standard 0-based half-open coordinates.  The first column must match the full path name of the target assembly contig in the graph (as shown in panpatch output).  Multiple regions can be specified, across multiple contigs, one per line.  Comment lines beginning with `#` are ignored.

## Telomeres

Use the `-T` option to require that each patched haplotype begins and ends with a telomere (and has none internally).  With `-T`, panpatch will:

1. **Patch a missing terminal telomere.**  If the target assembly stops short of a telomere at one end, panpatch looks for a higher-priority assembly (one of the other `-s` samples) that *does* reach a telomere there, and splices it on.  Because subtelomeric sequence is typically divergent (and absent from the reference path), the handoff is made at the nearest graph node shared by the target and the donor: the target's capless tip beyond that node is replaced by the donor's run to its telomere.  This is graph-coherent — the junction is a node both assemblies actually traverse.
2. **Reject bad patches.**  If, after any such patching, the assembly still does not begin and end with a telomere (or has an internal one — e.g. an alignment placed a telomere mid-assembly), the patch is reverted to the input contigs.

Without `-T`, telomeres are neither patched nor checked.

Because the replaced tip is the *divergent* subtelomere (everything past the last shared node), a telomere patch is not a swap of equivalent sequence — it completes the arm with the donor's version.  Two controls make this safe and auditable:

- `-M, --max-telomere-patch N` caps how much target sequence a single telomere patch may replace (default 500000).  The nearest shared handoff can be far in when the subtelomere is large (e.g. acrocentric arms); a patch that would replace more than `N` bp is skipped and the assembly is left to the normal revert.  The skipped patch appears as a rejected `telomere` row whose `reason` reports the distance, so you can opt in with a larger `-M`.
- Each applied patch is a `telomere` row in the report table; its `kmer%` column is the fraction of the replaced target sequence's k-mers also found in the graft — a low value flags that the donor's subtelomere differs substantially from the target's (telomere patches are accepted regardless, since the subtelomere is expected to be divergent).

When an end lacks a telomere and panpatch cannot lift one over, it appears as a rejected `telomere` row whose `reason` distinguishes a simple gap (no telomere here and no donor reaches one) from an assembly issue beyond panpatch's scope (telomeric repeats present near the tip but not as a clean terminal telomere — i.e. sequence extending past the telomere, or a degraded/fragmented one).

## Why a patch is rejected (quality control)

panpatch only emits a patched sequence when it passes a series of checks; otherwise it reverts to the target's input contig(s) — or to `--default-sample`'s contig if that option is given.  Every candidate patch appears in the stdout report table with `decision` = `accepted`/`rejected` and, when rejected, the cause in the `reason` column, so the outcome is auditable.

In the report table each candidate patch is a row with `decision` = `accepted` or `rejected`; the rejection cause appears in the `reason` column (shown below for each check).  A whole-contig revert prefixes its rows with `contig reverted: `.

**Checks applied to every run**

- **Patch too short.**  If the patched sequence covers less than `--min-cover` (default 0.95) of the combined length of the target sample's input contigs, the whole contig is reverted.
  reason: `contig reverted: patch covers only <frac> of target (< --min-cover)`
- **Repeat-region misjoin — same-sample splice.**  If a contig is used non-contiguously with another contig *of the same sample* spliced into the interior between its pieces, the patch is reverted.  This catches the failure mode where non-unique satellite / segmental-duplication anchors cause the target assembly's own spare fragments to be stitched into the middle of a contig that already spans the region — observed on acrocentric short arms and pericentromeres, where it collapses or scrambles megabase satellite arrays.  Legitimate operations are unaffected: a genuine gap-fill is bridged by a *foreign* donor (different sample), and an end-to-end scaffold uses each contig as a single contiguous block.
  reason: `contig reverted: contig <C> was used non-contiguously with same-sample fragment <F> spliced into its interior (repeat-region misjoin)`
- **Repeat-region misjoin — foreign graft.**  When a *foreign* donor is spliced into a target contig's interior, the graft is **excised** — dropped, with the target's own sequence spliced back in, while the rest of the patch (telomere completions, faithful fills) is *kept* — if it fails the test appropriate to its kind:
  - **content** (a replacement of real, non-N sequence, ≥ `--graft-min-bp` default 10000): it must recapitulate ≥ `--graft-recovery` (default 50%) of the replaced sequence's k-mers, else the donor came from a different locus/array;
  - **anchoring** (an N-gap fill, where there is no replaced sequence to recapitulate): the donor must stay homologous to the target's *own* sequence over `--flank-window` (default 500000) bp on **both** flanks, sharing ≥ `--min-flank` (default 50%); otherwise the fill sits on a structural disagreement (a wrong-locus join).
  Each graft is judged independently, so a faithful graft is kept even when a sibling graft on the same contig is excised.  The excised graft shows as a rejected `gap-fill` row; the kept ones stay `accepted`.  (Anything that doesn't fit the clean excision shape is still backstopped by a full revert.)
  reason: `k-mer recovery <r>% < 50% (repeat-region misjoin)`  or  `flank anchoring <l>%/<r>% < 50% (wrong-locus join)`
- **Repeat-region misjoin — foreign bridge.**  When a foreign donor bridges *two* of the target's own contigs into a scaffold (`A | donor | B`), the bridging donor must anchor to both contigs' own sequence over `--flank-window` bp (≥ `--min-flank`); otherwise it is a wrong-locus join and the whole scaffold is reverted.
  reason: `contig reverted: contig <A> scaffolded to <B> by foreign <D>: foreign bridge anchoring <l>%/<r>% < 50% (wrong-locus join)`
- **Discarded telomere.**  A patch must not trim off a target contig's telomere-bearing end.  If a contig is capped at a natural end (telomere density ≥ `--telomere-threshold`) but the patch uses it starting — or ending — ≥20 kb past that cap, the contig was already complete there, so the join is redundant/erroneous and is reverted.  (A divergent reference's subtelomere can otherwise make panpatch trim a real telomere to bolt on an overlapping same-haplotype fragment, when the main contig alone was already T2T.)  Telomere patches (which trim a *capless* tip), gap-fills (interior only), and end-to-end scaffolds (each contig used in full) are unaffected.
  reason: `contig reverted: patch trimmed the telomere-bearing 5'|3' end of <C> ... (target was already capped)`
- **Nothing to patch.**  If no foreign sequence was used and fewer than two of the target's own contigs were joined, there was no patch — so the contig produces *no* report row.  Its `#Contig <name> ... left=/right=` cap-status line (printed after the table) shows whether it is already complete (no `N` gaps) or simply had no donor that covered its gaps.

**Additional checks with `-T` (telomere requirements)**

- **Telomere validation failed.**  After any telomere patching, the assembly must begin and end with a telomere (hexamer density ≥ `--telomere-threshold`, default 0.8), contain none internally, and be at least 2 kb long; otherwise the contig is reverted.
  reason: `contig reverted: telomere validation failed`
- **Telomere patch over the cap.**  A missing terminal telomere is grafted from a donor only if the handoff replaces at most `--max-telomere-patch` (`-M`, default 500000) bp of target sequence; a larger graft is skipped, shown as a rejected `telomere` row.
  reason: `handoff would replace <N>bp, over --max-telomere-patch <M>`
- **No telomere to lift over.**  If an end has no telomere and no donor can supply one (or the telomere is degraded/buried — beyond panpatch's scope), it is a rejected `telomere` row.
  reason: `no telomere at this end and no donor reaches one`  (or a buried/degraded variant)

**Options that limit what is patched**

- `-b, --exclude-bed FILE` — keep listed target regions untouched (see *Excluding Regions from Patching*).
- `-e, --default-sample STRING` — when a patch is rejected, output this sample's contig for the chromosome instead of the target's input contigs.

### Running time

The above examples take about 2 hours on the cluster to run `cactus-pangenome`.  The single `panpatch` command over all the chromosome graphs takes about 2 minutes on my desktop. 

### Algorithm

All contigs are first binned by haplotype.  Since the input not necessarily trio-phased, this determines whether, for example, haplotype 1 from verkko corresponds to haplotype 1 or haplotype 2 from hifasm, etc.

This is accomplished by looking at the average alignment identity in the graph between pairs of haplotypes, over windows of `1000bp`.

<img src="panpatch-1.png" height=60% width=60%>

Next, the reference path of the chromosome (ie CHM13) is scanned left to right for potential anchors.  An anchor is a node on the reference path where one or more paths either starts, ends or branches off from another.

Finally a path through the anchors is searched in the graph that connects the first and last anchors (tips of the reference path), giving a T2T assembly of the contig (if possible).  The path stays on the highest priority path at every junction (verkko, then hifiasm, then duplex in our example). 

<img src="panpatch-2.png" height=60% width=60%>

### Limitations and future work

* Needs better checking for obviously bad patches:
     * Are the supporting alignments sketchy?
     * Are sequences being patched in unreasonably long?
     * Etc.
* Entirely reference-based.  If graph doesn't align contigs to reference, then no anchors will be found.  This could happen in acrocentric short arms, for example.
* Left-to-right reference-based graph search is very simplistic, and some cases could probably be improved with more general search.



