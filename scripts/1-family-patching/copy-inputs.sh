#!/usr/bin/bash
set -ex

./split-diploid.sh  /private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/diploid/assembly.v1.0.PAN027.diploid.fa maternal paternal &
./split-diploid.sh  /private/groups/migalab/mcechova/pedigree_assemblies/assembly_comparisons/PAN027.verkko2.1.Q28_k50000.diploid.20250425.fasta maternal paternal &
./split-diploid.sh /private/groups/migalab/mcechova/pedigree_assemblies/verkko2.0/assembly_release/diploid/assembly.v1.0.PAN011.diploid.fa haplotype1 haplotype2 &
wait

