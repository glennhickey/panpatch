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

Begin by making the graph. You specify the inputs to Cactus with a two-column file  `pan28.hs1.seqfile`:

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
cactus-pangenome ./js ./pan028.hs1.seqfile --outName pan028-mc-verkko-1 --outDir pan028-mc-verkko-1 --logFile pan028-mc-hs1.log --reference PAN028-verkko_1 --chrom-vg full --batchSystem slurm --consCores 65 --mgCores 64 --indexCores 64 --mapCores 16
```

Now you can run `panpatch` individually on each `.vg` file in `pan028-mc-verkko-1/pan028-mc-verkko-1.chroms/` :
```
cd pan028-mc-verkko-1/pan028-mc-verkko-1.chroms/
for CHR in *.vg; do \
    panpatch $CHR -r PAN028-verkko_1 -p -s PAN028-verkko_1 -s PAN028-hifiasm -s PAN028-duplex -f ${CHR::-3}-patched.fa > ${CHR::-3}-patched.bed 2>${CHR::-3}-patched.stderr ; \
done 
```

This will produce a FASTA file for each chromosome, as well as a BED file listing the patched regions.  IF a chromosome couldn't be patched, the FASTA output will be the same as the input.  The `.stderr` files will contain additional information about what wasn't patched and why.

Note that even though the `hifiasm` and `duplex` assemblies are diploid, only the most relevant haplotype for each will be selected for each chromosome.

When running on the second haplotype, the only difference is the first line of the seqfile, as well as the `--reference / -r` options.

## Scaffolding

If the assembly you want to patch does not have chromosome-scale scaffolds, you must use a reference that does.  Here is an example of using T2T-CHM13 (aka `hs1`) as the reference to patch `PAN028-verkko`.  In this case, since we only have one reference we can do all the patching at once, with a single graph.  For example, use `pan28.hs1.seqfile`:

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

Now you can run `panpatch` individually on each `.vg` file in `pan028-mc-hs1/pan028-mc-hs1.chroms/` :
```
cd pan028-mc-hs1/pan028-mc-hs1.chroms/
for CHR in *.vg; do \
    panpatch $CHR -r hs1 -p -s PAN028-verkko -s PAN028-hifiasm -s PAN028-duplex -f ${CHR::-3}-patched.fa > ${CHR::-3}-patched.bed 2>${CHR::-3}-patched.stderr ; \
done 
```

Note that unlike the above example, the output here will be diploid since both verkko haplotypes are getting patched. 

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

You specify the graph and sample names in order of priority (first column of the above input file, excluding `.1/2` suffixes)
```
panpatch <graph.vg> -r <reference sample> -s <sample to patch> -s <first sample to patch with> -s <second sample> etc.
```

For example

```
panpatch chr20.full.vg -r hs1 -s PAN028-verkko -s PAN028-hifiasm -s PAN028-duplex
```

will patch the `PAN028-verkko` assembly, using `PAN028-hifiasm` where possible, then `PAN028-duplex` as a backup.

The output will be a list of contig intervals (BED format), for each haplotype, that span the reference chromosome from telomere to telomere, which form the patched T2T assembly.

You can write a FASTA file for the patched contigs with `--fasta FILE`. 

Note: small intervals should probably be filtered out, there's not such logic yet in `panpatch`.  The output of the above is

```
Patched assembly for PAN028-verkko#1:
PAN028-hifiasm#2#h2tg000032l#0	2	27
PAN028-verkko#1#haplotype1-0000008#0	0	64821530
PAN028-verkko#1#haplotype1-0000043#0	11	673283
PAN028-hifiasm#2#h2tg000052l#0	649569	683780
PAN028-verkko#1#haplotype1-0000046#0	0	682357

Patched assembly for PAN028-verkko#2:
PAN028-verkko#2#haplotype2-0000073#0	1	65910828
```

## Check Telomeres

In some cases, the best alignment / patch can place a telomere inside the output assembly (as opposed to at the tip).  Use the `-T` option to check for these cases and reject patches that do not begin and end with telomeres. 

### Running time

The above examples takes about 2 hours on the cluster to run `cactus-pangenome`.  Running `panpatch` on each chromsome in series takes about 2 minutes total on my desktop. 

### Algorithm

All contigs are first binned by haplotype.  Since the input not necessarily trio-phased, this determines whether, for example, haplotype 1 from verkko corresponds to haplotype 1 or haplotype 2 from hifasm, etc.

This is accomplished by looking at the average alignment identity in the graph between pairs of haplotypes, over windows of `1000bp`.

<img src="panpatch-1.png" height=60% width=60%>

Next, the reference path of the chromosome (ie CHM13) is scanned left to right for potential anchors.  An anchor is a node in on the reference path that where one more paths either starts, ends or branches off from another.

Finally a path through the anchors is searched in the graph that connects the first and last anchors (tips of the reference path), giving a T2T assembly of the contig (if possible).  The path stays on the highest priority path at every junction (verkko, then hifiasm, then duplex in our example). 

<img src="panpatch-2.png" height=60% width=60%>

### Limitations and future work

* Needs better checking for obviously bad patches:
     * Are the supporting alignments sketchy?
     * Are sequences being patched in unreasonably long?
     * Etc.
* Entirely reference-based.  If graph doesn't align contigs to reference, then no anchors will be found.  This could happen in acrocentric short arms, for example.
* Left-to-right reference-based graph search is very simplistic, and some cases could probably be improved with more general search.



