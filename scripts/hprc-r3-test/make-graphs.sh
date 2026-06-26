#!/usr/bin/bash

set -ex

for SAMPLE in HG00235 HG00126 HG01074
do
    cactus-pangenome js-${SAMPLE} ${SAMPLE}.seqfile --outDir mc-${SAMPLE} --outName ${SAMPLE} --reference CHM13 --chrom-vg full --gfa full --logFile mc-${SAMPLE}.log --disableProgress --batchSystem slurm --slurmTime 10:00:00 --doubleMem true --mgCores 32 --mapCores 16 --consCores 16 --indexCores 32 &
done
wait

