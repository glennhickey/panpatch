#!/usr/bin/bash

set -ex

for SAMPLE in HG00235 HG00126 HG01074
do
    cactus-pangenome js-${SAMPLE}-nc1 ${SAMPLE}.no-chm13.1.seqfile --outDir mc-${SAMPLE}.no-chm13.1 --outName ${SAMPLE}.no-chm13.1 --reference verkko-R3-${SAMPLE}_1 --chrom-vg full --gfa full --logFile mc-${SAMPLE}.no-chm13.1.log --disableProgress --batchSystem slurm --slurmTime 10:00:00 --doubleMem true --mgCores 32 --mapCores 16 --consCores 16 --indexCores 32 --retryCount 10 --maxMemory 1.5T &

    cactus-pangenome js-${SAMPLE}-nc2 ${SAMPLE}.no-chm13.2.seqfile --outDir mc-${SAMPLE}.no-chm13.2 --outName ${SAMPLE}.no-chm13.2 --reference verkko-R3-${SAMPLE}_2 --chrom-vg full --gfa full --logFile mc-${SAMPLE}.no-chm13.2.log --disableProgress --batchSystem slurm --slurmTime 10:00:00 --doubleMem true --mgCores 32 --mapCores 16 --consCores 16 --indexCores 32 --retryCount 10 --maxMemory 1.5T &
    
done
wait

