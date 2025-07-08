#!/usr/bin/bash
set -ex

PANPATCH=/private/home/ghickey/dev/work/panpatch/panpatch
VERSION=v0.1

SAMPLE=PAN027
mkdir -p family-patches.${VERSION}
for HAP in 1 2; do
    for VG in ${SAMPLE}-family-mc-${HAP}/${SAMPLE}-family-mc-${HAP}.chroms/*.full.vg; do
	$PANPATCH $VG -t 4 -r PAN027-verkko-${HAP} -s ${SAMPLE}-verkko-${HAP} -s ${SAMPLE}-herro -s PAN011-verkko -f family-patches.${VERSION}/${SAMPLE}.$(basename $VG).${HAP}.patch.fa -p  > family-patches.${VERSION}/${SAMPLE}.$(basename $VG).${HAP}.patch.bed 2> family-patches.${VERSION}/${SAMPLE}.$(basename $VG).${HAP}.patch.stderr && bgzip -f family-patches.${VERSION}/${SAMPLE}.$(basename $VG).${HAP}.patch.fa &
    done
    wait
    rm -f family-patches.${VERSION}/${SAMPLE}.${HAP}.patch.fa.gz family-patches.${VERSION}/${SAMPLE}.${HAP}.patch.bed family-patches.${VERSION}/${SAMPLE}.${HAP}.patch.stderr
    for VG in ${SAMPLE}-family-mc-${HAP}/${SAMPLE}-family-mc-${HAP}.chroms/*.full.vg; do
	cat family-patches.${VERSION}/${SAMPLE}.$(basename $VG).${HAP}.patch.fa.gz >> family-patches.${VERSION}/${SAMPLE}.${HAP}.patch.fa.gz
	cat family-patches.${VERSION}/${SAMPLE}.$(basename $VG).${HAP}.patch.bed >> family-patches.${VERSION}/${SAMPLE}.${HAP}.patch.bed
	cat family-patches.${VERSION}/${SAMPLE}.$(basename $VG).${HAP}.patch.stderr >> family-patches.${VERSION}/${SAMPLE}.${HAP}.patch.stderr
	rm -f family-patches.${VERSION}/${SAMPLE}.$(basename $VG).${HAP}.patch.fa.gz family-patches.${VERSION}/${SAMPLE}.$(basename $VG).${HAP}.patch.bed family-patches.${VERSION}/${SAMPLE}.$(basename $VG).${HAP}.patch.stderr
    done
done
