#!/usr/bin/bash

for SAMPLE in PAN010 PAN011 PAN027 PAN028; do
    cactus-pangenome ./js-${SAMPLE} ./${SAMPLE}-hs1.seqfile --outName ${SAMPLE}-mc-hs1 --outDir ${SAMPLE}-mc-hs1 --logFile ${SAMPLE}-mc-hs1.log --batchSystem slurm --gbz --reference hs1 --consCores 65 --mgCores 64 --indexCores 64 --mapCores 16 --xg full --chrom-vg full --slurmTime 100:00:00 --lastTrain --maxLen 100000 2> ${SAMPLE}-mc-hs1.stderr &
done

wait

    
