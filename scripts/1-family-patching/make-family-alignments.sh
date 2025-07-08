#!/usr/bin/bash

for HAP in 1 2; do
    cactus-pangenome ./js-${HAP} ./PAN027-family-${HAP}.seqfile --outName PAN027-family-mc-${HAP} --outDir PAN027-family-mc-${HAP} --logFile PAN027-family-mc-${HAP}.log --batchSystem slurm --gbz --reference PAN027-verkko-${HAP} --consCores 8 --mgCores 16 --indexCores 16 --mapCores 8 --xg full --chrom-vg full --slurmTime 100:00:00 --lastTrain --maxLen 100000 2> PAN027-family-mc-${HAP}.stderr &
done

wait

    
