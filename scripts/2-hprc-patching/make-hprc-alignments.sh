#!/usr/bin/bash
set -ex
SAMPLE=$1

cactus-pangenome ./js-${SAMPLE}-hs1 ./${SAMPLE}-hs1.seqfile --outName ${SAMPLE}-mc-hs1 --outDir ${SAMPLE}-mc-hs1 --logFile ${SAMPLE}-mc-hs1.log --batchSystem slurm --gbz --reference hs1 --consCores 8 --mgCores 26 --indexCores 16 --mapCores 16 --xg full --chrom-vg full --slurmTime 100:00:00 --lastTrain --maxLen 100000 --doubleMem true --retryCount 10 2> ${SAMPLE}-mc-hs1.stderr &

#cactus-pangenome ./js-${SAMPLE}-1 ./${SAMPLE}-1.seqfile --outName ${SAMPLE}-mc-1 --outDir ${SAMPLE}-mc-1 --logFile ${SAMPLE}-mc-1.log --batchSystem slurm --gbz --reference ${SAMPLE}-verkko-1 --consCores 8 --mgCores 16 --indexCores 16 --mapCores 16 --xg full --chrom-vg full --slurmTime 100:00:00 --lastTrain --maxLen 100000 --doubleMem true --retryCount 10 2> ${SAMPLE}-mc-1.stderr &

#cactus-pangenome ./js-${SAMPLE}-2 ./${SAMPLE}-2.seqfile --outName ${SAMPLE}-mc-2 --outDir ${SAMPLE}-mc-2 --logFile ${SAMPLE}-mc-2.log --batchSystem slurm --gbz --reference ${SAMPLE}-verkko-2 --consCores 8 --mgCores 16 --indexCores 16 --mapCores 16 --xg full --chrom-vg full --slurmTime 100:00:00 --lastTrain --maxLen 100000 --doubleMem true --retryCount 10 2> ${SAMPLE}-mc-2.stderr &

wait

    
