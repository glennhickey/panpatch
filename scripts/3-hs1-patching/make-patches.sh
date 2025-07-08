#!/usr/bin/bash
set -ex

PANPATCH=/private/home/ghickey/dev/work/panpatch/panpatch
VERSION=v0.1

mkdir -p patches.${VERSION}
for SAMPLE in PAN010 PAN011 PAN027 PAN028; do
    for VG in ${SAMPLE}-mc-hs1/${SAMPLE}-mc-hs1.chroms/chr*.full.vg; do
	$PANPATCH $VG -t 4 -r hs1 -s ${SAMPLE}-verkko -s ${SAMPLE}-hifiasm -s ${SAMPLE}-duplex -f patches.${VERSION}/${SAMPLE}.$(basename $VG).patch.fa -p > patches.${VERSION}/${SAMPLE}.$(basename $VG).patch.bed 2> patches.${VERSION}/${SAMPLE}.$(basename $VG).patch.fa.stderr && bgzip -f patches.${VERSION}/${SAMPLE}.$(basename $VG).patch.fa &
    done
    wait
    rm -f patches.${VERSION}/${SAMPLE}.patch.fa.gz patches.${VERSION}/${SAMPLE}.patch.bed patches.${VERSION}/${SAMPLE}.patch.stderr
    cat patches.${VERSION}/${SAMPLE}*.fa.gz > patches.${VERSION}/${SAMPLE}.patch.fa.gz
    cat patches.${VERSION}/${SAMPLE}*.bed > patches.${VERSION}/${SAMPLE}.patch.bed
    cat patches.${VERSION}/${SAMPLE}*.stderr > patches.${VERSION}/${SAMPLE}.patch.stderr
    rm -f patches.${VERSION}/${SAMPLE}*chr*full*
done
