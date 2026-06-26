#!/usr/bin/python3

import os, sys
import subprocess

temp_files = []

fail_count = 0

def run(c):
    global fail_count
    try:
        return subprocess.check_call(c, shell=True)
    except:
        print('failed cmd: {}'.format(c))
        fail_count += 1
        
run('vg mod -X2 tiny.gfa | vg convert - > tiny.x2.vg')
run('vg convert tiny.gfa > tiny.vg')
run('echo ">chrX_hap_1" > tiny.truth.fa')
run('vg paths -x tiny.gfa -Q hifi -F | tail -1 >> tiny.truth.fa')
run('echo ">chrX_hap_1" > tiny.longo.truth.fa')
run('vg paths -x tiny.gfa -Q plongo -F | tail -1 >> tiny.longo.truth.fa')
temp_files += ['tiny.x2.vg', 'tiny.vg', 'tiny.truth.fa', 'tiny.longo.truth.fa']

run('panpatch tiny.vg -r x -s verkko -s hifi -f tiny.hifi.fa > tiny.hifi.bed')
run('diff tiny.hifi.bed tiny.hifi.bed.truth')
run('diff tiny.hifi.fa tiny.truth.fa')
temp_files += ['tiny.hifi.bed', 'tiny.hifi.fa']

run('panpatch tiny.x2.vg -r x -s verkko -s hifi -f tiny.x2.hifi.fa > tiny.x2.hifi.bed')
run('diff tiny.x2.hifi.bed tiny.x2.hifi.bed.truth')
run('diff tiny.x2.hifi.fa tiny.truth.fa')
temp_files += ['tiny.x2.hifi.bed', 'tiny.x2.hifi.fa']

run('panpatch tiny.vg -r x -s verkko -s duplex -f tiny.duplex.fa > tiny.duplex.bed')
run('diff tiny.duplex.bed tiny.duplex.bed.truth')
run('diff tiny.duplex.fa tiny.truth.fa')
temp_files += ['tiny.duplex.bed', 'tiny.duplex.fa']

run('panpatch tiny.x2.vg -r x -s verkko -s duplex -f tiny.x2.duplex.fa > tiny.x2.duplex.bed')
run('diff tiny.x2.duplex.bed tiny.x2.duplex.bed.truth')
run('diff tiny.x2.duplex.fa tiny.truth.fa')
temp_files += ['tiny.x2.duplex.bed', 'tiny.x2.duplex.fa']

run('panpatch tiny.vg -r x -s backo -s hifi -f tiny.backo.hifi.fa > tiny.backo.hifi.bed')
run('diff tiny.backo.hifi.bed tiny.backo.hifi.bed.truth')
run('diff tiny.backo.hifi.fa tiny.truth.fa')
temp_files += ['tiny.backo.hifi.bed', 'tiny.backo.hifi.fa']

run('panpatch tiny.vg -r x -s backo -s duplex -f tiny.backo.duplex.fa > tiny.backo.duplex.bed')
run('diff tiny.backo.duplex.bed tiny.backo.duplex.bed.truth')
run('diff tiny.backo.duplex.fa tiny.truth.fa')
temp_files += ['tiny.backo.duplex.bed', 'tiny.backo.duplex.fa']

run('panpatch tiny.vg -r x -s longo -s hifi -f tiny.longo.hifi.fa > tiny.longo.hifi.bed')
run('diff tiny.longo.hifi.bed tiny.longo.hifi.bed.truth')
run('diff tiny.longo.hifi.fa tiny.longo.truth.fa')
temp_files += ['tiny.longo.hifi.bed', 'tiny.longo.hifi.fa']

run('panpatch tiny.vg -r x -s longo -s duplex -f tiny.longo.duplex.fa > tiny.longo.duplex.bed')
run('diff tiny.longo.duplex.bed tiny.longo.duplex.bed.truth')
run('diff tiny.longo.duplex.fa tiny.longo.truth.fa')
temp_files += ['tiny.longo.duplex.bed', 'tiny.longo.duplex.fa']

# Telomere tests
run('vg convert telomere.gfa > telomere.vg')
temp_files += ['telomere.vg']

run('panpatch telomere.vg -r x -s verkko -s hifi > telomere.hifi.bed')
run('diff telomere.hifi.bed telomere.hifi.bed.truth')
temp_files += ['telomere.hifi.bed']

run('panpatch telomere.vg -r x -s verkko -s noTelo -T > telomere.noTelo.T.bed')
run('diff telomere.noTelo.T.bed telomere.noTelo.T.bed.truth')
temp_files += ['telomere.noTelo.T.bed']

# No-gaps message test: hifi has no N bases, so should get friendly message
run('panpatch tiny.vg -r x -s hifi -s verkko > tiny.nogaps.bed')
run('diff tiny.nogaps.bed tiny.nogaps.bed.truth')
temp_files += ['tiny.nogaps.bed']

# BED exclusion test: exclude the gap region so it doesn't get patched
run('printf "verkko#1#chrX#0\\t12\\t36\\n" > tiny.exclude.bed')
run('panpatch tiny.vg -r x -s verkko -s hifi -b tiny.exclude.bed > tiny.exclude.bed.out')
run('diff tiny.exclude.bed.out tiny.exclude.bed.truth')
temp_files += ['tiny.exclude.bed', 'tiny.exclude.bed.out']

# Scaffolding test: two contigs of the SAME target sample (frag) span the chromosome
# and meet at one shared anchor node, with no foreign sample bridging them. This
# exercises the target-only scaffold join that used to be incorrectly reverted as
# "no patching required".
run('vg convert scaffold.gfa > scaffold.vg')
temp_files += ['scaffold.vg']

# without -T: ctgA and ctgB are stitched into a single chrX_hap_1 record
run('panpatch scaffold.vg -r x -s frag -f scaffold.frag.fa > scaffold.frag.bed')
run('diff scaffold.frag.bed scaffold.frag.bed.truth')
run('diff scaffold.frag.fa scaffold.frag.truth.fa')
temp_files += ['scaffold.frag.bed', 'scaffold.frag.fa']

# with -T: the join has no telomeres, so it must fail validation and revert to the inputs
run('panpatch scaffold.vg -r x -s frag -T > scaffold.frag.T.bed')
run('diff scaffold.frag.T.bed scaffold.frag.T.bed.truth')
temp_files += ['scaffold.frag.T.bed']

# Interior-splice guard: the main contig (frag#1#main) spans the whole chromosome but diverges
# in the middle (node 8); a same-sample fragment (frag#1#alt) matches the reference there. Greedy
# threading splices alt into main's interior (main used non-contiguously with alt between), which
# is the repeat-region misjoin signature. The guard must reject this and revert to the input
# contigs. The big flanks keep the length check satisfied so the guard is the sole reason to revert
# (contrast with the legit cases above: end-to-end scaffolds and foreign-donor gap-fills are kept).
run('vg convert interior_splice.gfa > interior_splice.vg')
run('panpatch interior_splice.vg -r x -s frag > interior_splice.bed')
run('diff interior_splice.bed interior_splice.bed.truth')
temp_files += ['interior_splice.vg', 'interior_splice.bed']

# Telomere-preservation guard: frag#1#main spans the chromosome but its 5' telomere is on an
# off-reference node, while the same-sample fragment frag#1#tip carries the reference's 5' node.
# Greedy threading therefore trims main's telomere-bearing 5' end to use tip there (the chr14
# scaffold pattern). The guard must reject discarding a real telomere and revert to the inputs.
# Big interior nodes keep the length check satisfied so the guard is the sole reason to revert.
run('vg convert telo_preserve.gfa > telo_preserve.vg')
run('panpatch telo_preserve.vg -r x -s frag -T > telo_preserve.bed')
run('diff telo_preserve.bed telo_preserve.bed.truth')
temp_files += ['telo_preserve.vg', 'telo_preserve.bed']

# Telomere patching (folded into -T): when the target is missing a telomere, graft one in
# from a foreign cover (donor) by handing off at the nearest shared node. Small graphs
# exercise the splice paths. The non-deterministic-order #Contig lines are dropped, and the
# #Telomere patch line is kept (pinning the replaced=/grafted= byte counts) with only the
# kmer_recovery percentage masked.
#   telopatch       - forward target, telomere appended at the back
#   telopatch_rev   - reverse-oriented target (as in real assemblies), back patch
#   telopatch_front - forward target, telomere prepended at the front
#   telopatch_multi - scaffolded 2-interval target, front patch with handoff at the inner
#                     boundary (exercises the multi-interval rebuild + empty-interval trim)
#   telopatch_buried- target's back telomere is buried under terminal junk (verkko over-extension
#                     pattern): must read as capless and get patched, where the old lenient
#                     "any 500bp window" tip check would call it present and skip -> revert.
telo_filter = 'grep -v "^#Contig" | sed -E "s/kmer_recovery=[0-9.]+%/kmer_recovery=NA/"'
for base in ['telopatch', 'telopatch_rev', 'telopatch_front', 'telopatch_multi', 'telopatch_buried']:
    run('vg convert {0}.gfa > {0}.vg'.format(base))
    run('panpatch {0}.vg -r x -s frag -s donor -T -f {0}.fa 2>/dev/null | {1} > {0}.patch.bed'.format(base, telo_filter))
    run('diff {0}.patch.bed {0}.bed.truth'.format(base))
    run('diff {0}.fa {0}.truth.fa'.format(base))
    temp_files += ['{0}.vg'.format(base), '{0}.patch.bed'.format(base), '{0}.fa'.format(base)]


for f in temp_files:
    if os.path.isfile(f):
        os.remove(f)

if fail_count:
    print('{} tests failed'.format(fail_count))
else:
    print('all tests passed')
    
sys.exit(fail_count)
