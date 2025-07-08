#!/usr/bin/bash
set -ex

SAMPLE=$1

VERKKO_URL_BASE="https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/6807247E-4F71-45D8-AECE-9E5813BA1D9F--verkko-v2.2.1-release2_asms"
HPRC_URL_BASE="https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION"

for HAP in 1 2; do
    aria2c -x8 -j8 -s8 ${VERKKO_URL_BASE}/${SAMPLE}/verkko-hi-c/${SAMPLE}.assembly.haplotype${HAP}.fasta.gz &
done
wait

for HAP in 1 2; do
    seqtk cutN ${SAMPLE}.assembly.haplotype${HAP}.fasta.gz | bgzip > ${SAMPLE}.assembly.haplotype${HAP}.cut.fasta.gz &
done
wait

for HAP in 1 2; do
    rm -f ${SAMPLE}.assembly.haplotype${HAP}.fasta.gz
done

printf "${SAMPLE}-verkko-cut.1  ${SAMPLE}.assembly.haplotype1.cut.fasta.gz\n" > ${SAMPLE}.seqfile
printf "${SAMPLE}-verkko-cut.2  ${SAMPLE}.assembly.haplotype2.cut.fasta.gz\n" >> ${SAMPLE}.seqfile

printf "${SAMPLE}-verkko.1  ${VERKKO_URL_BASE}/${SAMPLE}/verkko-hi-c/${SAMPLE}.assembly.haplotype1.fasta.gz\n" >> ${SAMPLE}.seqfile
printf "${SAMPLE}-verkko.2  ${VERKKO_URL_BASE}/${SAMPLE}/verkko-hi-c/${SAMPLE}.assembly.haplotype2.fasta.gz\n" >> ${SAMPLE}.seqfile

grep  "${SAMPLE}_" /private/groups/cgl/hprc-graphs/hprc-v2.0-feb28/hprc-v2.0-mc-chm13.seqfile >> ${SAMPLE}.seqfile

printf "hs1  https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz\n" > ${SAMPLE}-hs1.seqfile
cat ${SAMPLE}.seqfile >> ${SAMPLE}-hs1.seqfile

sed -e "s/${SAMPLE}-verkko.1/${SAMPLE}-verkko-1/g" ${SAMPLE}.seqfile | grep -v verkko.2 | grep -v cut.2 > ${SAMPLE}-1.seqfile
sed -e "s/${SAMPLE}-verkko.2/${SAMPLE}-verkko-2/g" ${SAMPLE}.seqfile | grep -v verkko.1 | grep -v cut.1 > ${SAMPLE}-2.seqfile
