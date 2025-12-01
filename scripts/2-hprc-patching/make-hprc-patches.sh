#!/usr/bin/bash
set -ex

PANPATCH=/private/home/ghickey/dev/work/panpatch/panpatch
SAMPLE=$1
VERSION=$2

mkdir -p hs1-patches.${VERSION}
for VG in ${SAMPLE}-mc-hs1/${SAMPLE}-mc-hs1.chroms/chr*.full.vg; do
    $PANPATCH $VG -t 4 -T -r hs1 -s ${SAMPLE}-verkko -s ${SAMPLE} -f hs1-patches.${VERSION}/${SAMPLE}.$(basename $VG).patch.fa -p > hs1-patches.${VERSION}/${SAMPLE}.$(basename $VG).patch.bed 2> hs1-patches.${VERSION}/${SAMPLE}.$(basename $VG).patch.fa.stderr && bgzip -f hs1-patches.${VERSION}/${SAMPLE}.$(basename $VG).patch.fa &
done
wait
rm -f hs1-patches.${VERSION}/${SAMPLE}.patch.fa.gz hs1-patches.${VERSION}/${SAMPLE}.patch.bed hs1-patches.${VERSION}/${SAMPLE}.patch.stderr
cat hs1-patches.${VERSION}/${SAMPLE}*.fa.gz > hs1-patches.${VERSION}/${SAMPLE}.patch.fa.gz
cat hs1-patches.${VERSION}/${SAMPLE}*.bed > hs1-patches.${VERSION}/${SAMPLE}.patch.bed
cat hs1-patches.${VERSION}/${SAMPLE}*.stderr > hs1-patches.${VERSION}/${SAMPLE}.patch.stderr
rm -f hs1-patches.${VERSION}/${SAMPLE}*chr*full*

#mkdir -p verkko-patches.${VERSION}
#for HAP in 1 2; do
#    for VG in ${SAMPLE}-mc-${HAP}/${SAMPLE}-mc-${HAP}.chroms/hap*.full.vg; do
#	$PANPATCH $VG -t 4 -r ${SAMPLE}-verkko-${HAP} -s ${SAMPLE}-verkko-${HAP} -s ${SAMPLE} -f verkko-patches.${VERSION}/${SAMPLE}.$(basename $VG).patch.${HAP}.fa -pe > verkko-patches.${VERSION}/${SAMPLE}.$(basename $VG).patch.${HAP}.bed 2> verkko-patches.${VERSION}/${SAMPLE}.$(basename $VG).patch.${HAP}.stderr && bgzip -f verkko-patches.${VERSION}/${SAMPLE}.$(basename $VG).patch.${HAP}.fa &
#    done
#    wait
#    rm -f verkko-patches.${VERSION}/${SAMPLE}.patch.${HAP}.fa.gz verkko-patches.${VERSION}/${SAMPLE}.patch.${HAP}.bed verkko-patches.${VERSION}/${SAMPLE}.patch.${HAP}.stderr
#    cat verkko-patches.${VERSION}/${SAMPLE}*.${HAP}.fa.gz > verkko-patches.${VERSION}/${SAMPLE}.patch.${HAP}.fa.gz
#    cat verkko-patches.${VERSION}/${SAMPLE}*.${HAP}.bed > verkko-patches.${VERSION}/${SAMPLE}.patch.${HAP}.bed
#    cat verkko-patches.${VERSION}/${SAMPLE}*.${HAP}.stderr > verkko-patches.${VERSION}/${SAMPLE}.patch.${HAP}.stderr
#    rm -f verkko-patches.${VERSION}/${SAMPLE}*hap*full*
#done
	   

