#!/usr/bin/python3

import os, sys
import random
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

run('panpatch tiny.vg -r x -s verkko -s hifi -f tiny.hifi.fa --bed tiny.hifi.bed.int > tiny.hifi.bed.rep && cat tiny.hifi.bed.rep tiny.hifi.bed.int > tiny.hifi.bed')
run('diff tiny.hifi.bed tiny.hifi.bed.truth')
run('diff tiny.hifi.hap1.fa tiny.truth.fa')
temp_files += ['tiny.hifi.bed', 'tiny.hifi.hap1.fa']

run('panpatch tiny.x2.vg -r x -s verkko -s hifi -f tiny.x2.hifi.fa --bed tiny.x2.hifi.bed.int > tiny.x2.hifi.bed.rep && cat tiny.x2.hifi.bed.rep tiny.x2.hifi.bed.int > tiny.x2.hifi.bed')
run('diff tiny.x2.hifi.bed tiny.x2.hifi.bed.truth')
run('diff tiny.x2.hifi.hap1.fa tiny.truth.fa')
temp_files += ['tiny.x2.hifi.bed', 'tiny.x2.hifi.hap1.fa']

run('panpatch tiny.vg -r x -s verkko -s duplex -f tiny.duplex.fa --bed tiny.duplex.bed.int > tiny.duplex.bed.rep && cat tiny.duplex.bed.rep tiny.duplex.bed.int > tiny.duplex.bed')
run('diff tiny.duplex.bed tiny.duplex.bed.truth')
run('diff tiny.duplex.hap1.fa tiny.truth.fa')
temp_files += ['tiny.duplex.bed', 'tiny.duplex.hap1.fa']

run('panpatch tiny.x2.vg -r x -s verkko -s duplex -f tiny.x2.duplex.fa --bed tiny.x2.duplex.bed.int > tiny.x2.duplex.bed.rep && cat tiny.x2.duplex.bed.rep tiny.x2.duplex.bed.int > tiny.x2.duplex.bed')
run('diff tiny.x2.duplex.bed tiny.x2.duplex.bed.truth')
run('diff tiny.x2.duplex.hap1.fa tiny.truth.fa')
temp_files += ['tiny.x2.duplex.bed', 'tiny.x2.duplex.hap1.fa']

run('panpatch tiny.vg -r x -s backo -s hifi -f tiny.backo.hifi.fa --bed tiny.backo.hifi.bed.int > tiny.backo.hifi.bed.rep && cat tiny.backo.hifi.bed.rep tiny.backo.hifi.bed.int > tiny.backo.hifi.bed')
run('diff tiny.backo.hifi.bed tiny.backo.hifi.bed.truth')
run('diff tiny.backo.hifi.hap1.fa tiny.truth.fa')
temp_files += ['tiny.backo.hifi.bed', 'tiny.backo.hifi.hap1.fa']

run('panpatch tiny.vg -r x -s backo -s duplex -f tiny.backo.duplex.fa --bed tiny.backo.duplex.bed.int > tiny.backo.duplex.bed.rep && cat tiny.backo.duplex.bed.rep tiny.backo.duplex.bed.int > tiny.backo.duplex.bed')
run('diff tiny.backo.duplex.bed tiny.backo.duplex.bed.truth')
run('diff tiny.backo.duplex.hap1.fa tiny.truth.fa')
temp_files += ['tiny.backo.duplex.bed', 'tiny.backo.duplex.hap1.fa']

run('panpatch tiny.vg -r x -s longo -s hifi -f tiny.longo.hifi.fa --bed tiny.longo.hifi.bed.int > tiny.longo.hifi.bed.rep && cat tiny.longo.hifi.bed.rep tiny.longo.hifi.bed.int > tiny.longo.hifi.bed')
run('diff tiny.longo.hifi.bed tiny.longo.hifi.bed.truth')
run('diff tiny.longo.hifi.hap1.fa tiny.longo.truth.fa')
temp_files += ['tiny.longo.hifi.bed', 'tiny.longo.hifi.hap1.fa']

run('panpatch tiny.vg -r x -s longo -s duplex -f tiny.longo.duplex.fa --bed tiny.longo.duplex.bed.int > tiny.longo.duplex.bed.rep && cat tiny.longo.duplex.bed.rep tiny.longo.duplex.bed.int > tiny.longo.duplex.bed')
run('diff tiny.longo.duplex.bed tiny.longo.duplex.bed.truth')
run('diff tiny.longo.duplex.hap1.fa tiny.longo.truth.fa')
temp_files += ['tiny.longo.duplex.bed', 'tiny.longo.duplex.hap1.fa']

# Path sense: vg makes every PanSN path of a plain GFA (as in the tests above) a haplotype, but
# cactus-pangenome graphs carry the reference as a REFERENCE-sense path. Tag x as the reference
# sample (RS:Z:x) and check that nothing changes.
with open('tiny.gfa') as gfa_in, open('tiny.ref.gfa', 'w') as gfa_out:
    gfa_lines = gfa_in.readlines()
    assert gfa_lines[0].startswith('H\t')
    gfa_out.write('H\tVN:Z:1.1\tRS:Z:x\n')
    gfa_out.writelines(gfa_lines[1:])
run('vg convert tiny.ref.gfa > tiny.ref.vg')
run('vg paths -M -x tiny.ref.vg | grep -q "^x#0#chrX\tREFERENCE\t"')
run('panpatch tiny.ref.vg -r x -s verkko -s hifi -f tiny.ref.hifi.fa --bed tiny.ref.hifi.bed.int > tiny.ref.hifi.bed.rep && cat tiny.ref.hifi.bed.rep tiny.ref.hifi.bed.int > tiny.ref.hifi.bed')
run('diff tiny.ref.hifi.bed tiny.hifi.bed.truth')
run('diff tiny.ref.hifi.hap1.fa tiny.truth.fa')
temp_files += ['tiny.ref.gfa', 'tiny.ref.vg', 'tiny.ref.hifi.bed', 'tiny.ref.hifi.hap1.fa']

# Telomere tests
run('vg convert telomere.gfa > telomere.vg')
temp_files += ['telomere.vg']

run('panpatch telomere.vg -r x -s verkko -s hifi --bed telomere.hifi.bed.int > telomere.hifi.bed.rep && cat telomere.hifi.bed.rep telomere.hifi.bed.int > telomere.hifi.bed')
run('diff telomere.hifi.bed telomere.hifi.bed.truth')
temp_files += ['telomere.hifi.bed']

run('panpatch telomere.vg -r x -s verkko -s noTelo -T --bed telomere.noTelo.T.bed.int > telomere.noTelo.T.bed.rep && cat telomere.noTelo.T.bed.rep telomere.noTelo.T.bed.int > telomere.noTelo.T.bed')
run('diff telomere.noTelo.T.bed telomere.noTelo.T.bed.truth')
temp_files += ['telomere.noTelo.T.bed']

# --- --patch-types: gap/telomere toggles (both on by default) ---
# default: telomere patching is on, but stays silent when there is nothing to complete (no telomere rows)
run('! (panpatch tiny.vg -r x -s verkko -s hifi 2>/dev/null | grep -q telomere)')
# -T additionally reports the per-end telomere failures (and requires a telomere at both ends)
run('panpatch tiny.vg -r x -s verkko -s hifi -T 2>/dev/null | grep -q telomere')
# --patch-types telomere turns gap-fill off: the tiny gap-fill is excised (marked excluded in the report)
run('panpatch tiny.vg -r x -s verkko -s hifi --patch-types telomere 2>/dev/null | grep -q "excluded by --patch-types"')
# -T with telomere excluded from --patch-types is contradictory and must error
run('! panpatch tiny.vg -r x -s verkko -s hifi -T --patch-types gap 2>/dev/null')

# No-gaps message test: hifi has no N bases, so should get friendly message
run('panpatch tiny.vg -r x -s hifi -s verkko --bed tiny.nogaps.bed.int > tiny.nogaps.bed.rep && cat tiny.nogaps.bed.rep tiny.nogaps.bed.int > tiny.nogaps.bed')
run('diff tiny.nogaps.bed tiny.nogaps.bed.truth')
temp_files += ['tiny.nogaps.bed']

# BED exclusion test: exclude the gap region so it doesn't get patched
run('printf "verkko#1#chrX#0\\t12\\t36\\n" > tiny.exclude.bed')
run('panpatch tiny.vg -r x -s verkko -s hifi -b tiny.exclude.bed --bed tiny.exclude.bed.out.int > tiny.exclude.bed.out.rep && cat tiny.exclude.bed.out.rep tiny.exclude.bed.out.int > tiny.exclude.bed.out')
run('diff tiny.exclude.bed.out tiny.exclude.bed.truth')
temp_files += ['tiny.exclude.bed', 'tiny.exclude.bed.out']

# Scaffolding test: two contigs of the SAME target sample (frag) span the chromosome
# and meet at one shared anchor node, with no foreign sample bridging them. This
# exercises the target-only scaffold join that used to be incorrectly reverted as
# "no patching required".
run('vg convert scaffold.gfa > scaffold.vg')
temp_files += ['scaffold.vg']

# without -T: ctgA and ctgB are stitched into a single chrX_hap_1 record
run('panpatch scaffold.vg -r x -s frag -f scaffold.frag.fa --bed scaffold.frag.bed.int > scaffold.frag.bed.rep && cat scaffold.frag.bed.rep scaffold.frag.bed.int > scaffold.frag.bed')
run('diff scaffold.frag.bed scaffold.frag.bed.truth')
run('diff scaffold.frag.hap1.fa scaffold.frag.truth.fa')
temp_files += ['scaffold.frag.bed', 'scaffold.frag.hap1.fa']

# with -T: the join has no telomeres, so it must fail validation and revert to the inputs
run('panpatch scaffold.vg -r x -s frag -T --bed scaffold.frag.T.bed.int > scaffold.frag.T.bed.rep && cat scaffold.frag.T.bed.rep scaffold.frag.T.bed.int > scaffold.frag.T.bed')
run('diff scaffold.frag.T.bed scaffold.frag.T.bed.truth')
temp_files += ['scaffold.frag.T.bed']

# Pass-through: a graph whose -r sample has more than one path (no single reference chromosome -- e.g. a
# chrOther graph of unplaced contigs) is not skipped; its target contigs are emitted unchanged (a warning
# on stderr, a "passthrough" row per contig, the contig in the FASTA). Here frag itself (ctgA + ctgB) is
# the reference, so both contigs pass through.
run('panpatch scaffold.vg -r frag -s frag -f scaffold.pass.fa --bed scaffold.pass.bed.int 2>/dev/null > scaffold.pass.bed.rep && cat scaffold.pass.bed.rep scaffold.pass.bed.int > scaffold.pass.bed')
run('diff scaffold.pass.bed scaffold.pass.bed.truth')
temp_files += ['scaffold.pass.bed', 'scaffold.pass.hap1.fa']

# Interior-splice guard: the main contig (frag#1#main) spans the whole chromosome but leaves the
# reference in the middle (node 8, an assembly gap of Ns); a same-sample fragment (frag#1#alt)
# matches the reference there. Threading can't follow main across the gap on its first pass, so it
# splices alt into main's interior (main used non-contiguously with alt between), which is the
# repeat-region misjoin signature. The gap is what makes this happen whichever of the two contigs
# ranks first (if main could be followed through node 8, the splice would depend on that order).
# The guard must reject this and revert to the input contigs. The big flanks keep the length check
# satisfied so the guard is the sole reason to revert (contrast with the legit cases above:
# end-to-end scaffolds and foreign-donor gap-fills are kept).
run('vg convert interior_splice.gfa > interior_splice.vg')
run('panpatch interior_splice.vg -r x -s frag --bed interior_splice.bed.int > interior_splice.bed.rep && cat interior_splice.bed.rep interior_splice.bed.int > interior_splice.bed')
run('diff interior_splice.bed interior_splice.bed.truth')
temp_files += ['interior_splice.vg', 'interior_splice.bed']

# Telomere-preservation guard: frag#1#main spans the chromosome but its 5' telomere is on an
# off-reference node, while the same-sample fragment frag#1#tip carries the reference's 5' node.
# Greedy threading therefore trims main's telomere-bearing 5' end to use tip there (the chr14
# scaffold pattern). The guard must reject discarding a real telomere and revert to the inputs.
# Big interior nodes keep the length check satisfied so the guard is the sole reason to revert.
run('vg convert telo_preserve.gfa > telo_preserve.vg')
run('panpatch telo_preserve.vg -r x -s frag -T --bed telo_preserve.bed.int > telo_preserve.bed.rep && cat telo_preserve.bed.rep telo_preserve.bed.int > telo_preserve.bed')
run('diff telo_preserve.bed telo_preserve.bed.truth')
temp_files += ['telo_preserve.vg', 'telo_preserve.bed']

# Telomere patching (folded into -T): when the target is missing a telomere, graft one in
# from a foreign cover (donor) by handing off at the nearest shared node. Small graphs
# exercise the splice paths. The #Contig lines are dropped; the telomere rows of the report
# are kept (pinning the replaced_bp/donor_bp byte counts and the k-mer recovery).
#   telopatch       - forward target, telomere appended at the back
#   telopatch_rev   - reverse-oriented target (as in real assemblies), back patch
#   telopatch_front - forward target, telomere prepended at the front
#   telopatch_multi - scaffolded 2-interval target, front patch with handoff at the inner
#                     boundary (exercises the multi-interval rebuild + empty-interval trim)
#   telopatch_buried- target's back telomere is buried under terminal junk (verkko over-extension
#                     pattern): must read as capless and get patched, where the old lenient
#                     "any 500bp window" tip check would call it present and skip -> revert.
# These grafts carry kb of donor subtelomere ahead of the telomere, so --max-subtelomere-graft is
# disabled to exercise the splices (the guard has its own cases, subtelo_* below).
telo_filter = 'grep -v "^#Contig"'
for base in ['telopatch', 'telopatch_rev', 'telopatch_front', 'telopatch_multi', 'telopatch_buried']:
    run('vg convert {0}.gfa > {0}.vg'.format(base))
    run('panpatch {0}.vg -r x -s frag -s donor -T --max-subtelomere-graft 0 --bed {0}.int -f {0}.fa 2>/dev/null | {1} > {0}.rep && cat {0}.rep {0}.int > {0}.patch.bed'.format(base, telo_filter))
    run('diff {0}.patch.bed {0}.bed.truth'.format(base))
    run('diff {0}.hap1.fa {0}.truth.fa'.format(base))
    temp_files += ['{0}.vg'.format(base), '{0}.patch.bed'.format(base), '{0}.hap1.fa'.format(base)]

# Telomere handoff must leave a capped end.  When the target's own short telomere array is aligned to
# the donor's (as cactus telomere anchors align them), the nearest shared node is inside the arrays:
# handing off there grafts a few bp, the end stays capless and -T reverts the contig.  The handoff is
# taken further in instead, where the graft carries the donor's array.  The graphs are written here:
#   node 1 (6bp) + node 2 (234bp) - frag's short p array, shared with the outer end of donor's
#   node 3 - frag-only junk between its array and the subtelomere
#   node 4 - the rest of donor's array (donor only)     node 8 - donor2's own array
#   nodes 5, 6, 7 - subtelomere, body and q telomere, the reference x
#   handoff_short   - frag 1 2 3 5 6 7, donor 1 2 4 5 6 7: the 6bp graft at node 2 is skipped, and the
#                     handoff is at node 5, grafting donor's whole array (1 2 4) over frag's 1 2 3
#   handoff_noarray - frag 5 6 7: no frag array, donor's array is grafted whole at the tip, as before
#   handoff_donor2  - handoff_short plus a lower-priority donor2 8 2 4 5 6 7: the nearest node still
#                     comes first, so once donor's graft at node 2 is rejected, donor2's (node 8) is taken
def handoff_gfa(path, walks):
    rng = random.Random(7)
    rnd = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
    nodes = {1: 'CCCTAA', 2: 'CCCTAA' * 39, 3: rnd(1000), 4: 'CCCTAA' * 210, 5: rnd(3000), 6: rnd(3000),
             7: 'TTAGGG' * 150, 8: 'CCCTAA' * 250}
    edges = sorted(set(e for w in walks.values() for e in zip(w, w[1:])))
    with open(path, 'w') as f:
        f.write('H\tVN:Z:1.1\n')
        for n in sorted(set(n for w in walks.values() for n in w)):
            f.write('S\t{}\t{}\n'.format(n, nodes[n]))
        for a, b in edges:
            f.write('L\t{}\t+\t{}\t+\t0M\n'.format(a, b))
        for sample, w in walks.items():
            f.write('W\t{}\t{}\t{}\t0\t{}\t{}\n'.format(sample, 0 if sample == 'x' else 1, 'chrX' if sample == 'x' else 'ctg',
                                                       sum(len(nodes[n]) for n in w), ''.join('>{}'.format(n) for n in w)))

handoff_cases = {
    'handoff_short'   : ({'x' : [5, 6, 7], 'frag' : [1, 2, 3, 5, 6, 7], 'donor' : [1, 2, 4, 5, 6, 7]},
                         'chrX\t1\ttelomere\tfrag#1#ctg#0\t8140\tdonor#1#ctg#0\t8400\t1240\t3.1\t.\t.\taccepted\t.\t0\t1240\n'
                         '#Patched assembly on chrX for frag#1:\n'
                         'donor#1#ctg#0\t0\t1500\t+\n'
                         'frag#1#ctg#0\t1240\t8140\t+\n\n'),
    'handoff_noarray' : ({'x' : [5, 6, 7], 'frag' : [5, 6, 7], 'donor' : [1, 2, 4, 5, 6, 7]},
                         'chrX\t1\ttelomere\tfrag#1#ctg#0\t6900\tdonor#1#ctg#0\t8400\t0\t0.0\t.\t.\taccepted\t.\t0\t0\n'
                         '#Patched assembly on chrX for frag#1:\n'
                         'donor#1#ctg#0\t0\t1500\t+\n'
                         'frag#1#ctg#0\t0\t6900\t+\n\n'),
    'handoff_donor2'  : ({'x' : [5, 6, 7], 'frag' : [1, 2, 3, 5, 6, 7], 'donor' : [1, 2, 4, 5, 6, 7],
                          'donor2' : [8, 2, 4, 5, 6, 7]},
                         'chrX\t1\ttelomere\tfrag#1#ctg#0\t8140\tdonor2#1#ctg#0\t9894\t6\t0.0\t.\t.\taccepted\t.\t0\t6\n'
                         '#Patched assembly on chrX for frag#1:\n'
                         'donor2#1#ctg#0\t0\t1500\t+\n'
                         'frag#1#ctg#0\t6\t8140\t+\n\n'),
}
for base, (walks, truth) in handoff_cases.items():
    handoff_gfa(base + '.gfa', walks)
    with open(base + '.truth', 'w') as f:
        f.write(truth)
    donors = ' '.join('-s ' + s for s in ['donor', 'donor2'] if s in walks)
    run('vg convert {0}.gfa > {0}.vg'.format(base))
    run('panpatch {0}.vg -r x -s frag {1} -T --bed {0}.int 2>/dev/null | grep -v "^chrom" | {2} > {0}.rep && cat {0}.rep {0}.int > {0}.out'.format(base, donors, telo_filter))
    run('diff {0}.out {0}.truth'.format(base))
    temp_files += ['{0}.gfa'.format(base), '{0}.vg'.format(base), '{0}.truth'.format(base), '{0}.out'.format(base)]

# A patch never takes a donor's N's (an assembly gap, or a region cactus-panpatch masked as an error).
# Donor N's used to come through three ways, one case each here:
#   ntelo_skip    - ndonor (higher priority) reaches the q telomere through an N node (7): its graft is
#                   skipped and donor's clean one is taken
#   ntelo_only    - only ndonor: the end is left as it is, with a rejected row saying why
#   ntelo_handoff - frag's tip is the N node, shared with ndonor: the handoff node is on the donor's side of
#                   the splice, so it counts too; the handoff moves in to node 3 and donor's run is grafted
#   ngap_shared   - reference-free: frag's N gap is two nodes, the first (13) shared with donor's N run,
#                   past which donor has sequence (15) where frag is still N (14).  no interval may start
#                   on the shared N anchor, so frag keeps its gap rather than donor being threaded from it
# (the ntelo grafts carry donor's 10kb node 4 ahead of its telomere, so --max-subtelomere-graft is off there)
# A telomere patch grafts the donor's telomere array, not a subtelomere the target lacks.  In n_telo, donor's
# graft is node 4 (10kb, no telomere) then its array (5); the telomeric run is found from 9904bp in:
#   subtelo_graft    - the default --max-subtelomere-graft (2000) rejects it, and the end is left as it is
#   subtelo_graft_ok - at --max-subtelomere-graft 9904 (the graft's exact subtelomere) it is taken
def n_gfa(path, nodes, walks):
    edges = sorted(set(e for w in walks.values() for e in zip(w, w[1:])))
    with open(path, 'w') as f:
        f.write('H\tVN:Z:1.1\n')
        for n in sorted(set(n for w in walks.values() for n in w)):
            f.write('S\t{}\t{}\n'.format(n, nodes[n]))
        for a, b in edges:
            f.write('L\t{}\t+\t{}\t+\t0M\n'.format(a, b))
        for sample, w in walks.items():
            f.write('W\t{}\t{}\t{}\t0\t{}\t{}\n'.format(sample, 0 if sample == 'x' else 1, 'chrX' if sample == 'x' else 'ctg',
                                                       sum(len(nodes[n]) for n in w), ''.join('>{}'.format(n) for n in w)))
n_rng = random.Random(11)
n_rnd = lambda n: ''.join(n_rng.choice('ACGT') for _ in range(n))
n_nodes = {1: 'CCCTAA' * 334, 2: n_rnd(10000), 3: n_rnd(10000), 4: n_rnd(10000), 5: 'TTAGGG' * 334, 6: n_rnd(3000),
           7: 'N' * 1000, 11: n_rnd(10000), 12: n_rnd(10000), 13: 'N' * 500, 14: 'N' * 700, 15: n_rnd(700), 16: n_rnd(10000)}
n_telo = {'x' : [1, 2, 3, 4, 5], 'frag' : [1, 2, 3, 6], 'ndonor' : [2, 3, 7, 4, 5], 'donor' : [2, 3, 4, 5]}
n_cases = {
    'ntelo_skip'    : (n_telo, '-r x -s frag -s ndonor -s donor --max-subtelomere-graft 0', (
                       'chrX\t1\ttelomere\tfrag#1#ctg#0\t25004\tdonor#1#ctg#0\t32004\t3000\t0.0\t.\t.\taccepted\t.\t22004\t25004\n'
                       '#Patched assembly on chrX for frag#1:\n'
                       'frag#1#ctg#0\t0\t12004\t+\n'
                       'donor#1#ctg#0\t10000\t32004\t+\n'
                       '\n')),
    'ntelo_only'    : (n_telo, '-r x -s frag -s ndonor', (
                       "chrX\t1\ttelomere\tfrag#1#ctg#0\t25004\tndonor#1#ctg#0\t33004\t3000\t.\t.\t.\trejected\tdonor's run to its telomere contains N (an assembly gap)\t22004\t25004\n"
                       '#Unpatched (input contig) on chrX for frag#1:\n'
                       'frag#1#ctg#0\t0\t25004\t+\n'
                       '\n')),
    'ntelo_handoff' : (dict(n_telo, frag=[1, 2, 3, 7]), '-r x -s frag -s ndonor -s donor --max-subtelomere-graft 0', (
                       'chrX\t1\ttelomere\tfrag#1#ctg#0\t23004\tdonor#1#ctg#0\t32004\t1000\t0.0\t.\t.\taccepted\t.\t22004\t23004\n'
                       '#Patched assembly on chrX for frag#1:\n'
                       'frag#1#ctg#0\t0\t12004\t+\n'
                       'donor#1#ctg#0\t10000\t32004\t+\n'
                       '\n')),
    'subtelo_graft' : (n_telo, '-r x -s frag -s donor', (
                       "chrX\t1\ttelomere\tfrag#1#ctg#0\t25004\tdonor#1#ctg#0\t32004\t3000\t.\t.\t.\trejected\tdonor's graft carries 9904bp ahead of its telomere, over --max-subtelomere-graft 2000\t22004\t25004\n"
                       '#Unpatched (input contig) on chrX for frag#1:\n'
                       'frag#1#ctg#0\t0\t25004\t+\n'
                       '\n')),
    'subtelo_graft_ok' : (n_telo, '-r x -s frag -s donor --max-subtelomere-graft 9904', (
                       'chrX\t1\ttelomere\tfrag#1#ctg#0\t25004\tdonor#1#ctg#0\t32004\t3000\t0.0\t.\t.\taccepted\t.\t22004\t25004\n'
                       '#Patched assembly on chrX for frag#1:\n'
                       'frag#1#ctg#0\t0\t12004\t+\n'
                       'donor#1#ctg#0\t10000\t32004\t+\n'
                       '\n')),
    'ngap_shared'   : ({'frag' : [11, 12, 13, 14, 16], 'donor' : [11, 12, 13, 15, 16]}, '-r frag -s frag -s donor', (
                       'ctg\t1\tgap-fill\tfrag#1#ctg#0\t31200\t.\t0\t.\t.\t.\t.\trejected\tno patches from other assemblies were found; contig kept as it is\t.\t.\n'
                       '#Unpatched (input contig) on ctg for frag#1:\n'
                       'frag#1#ctg#0\t0\t31200\t+\n'
                       '\n')),
}
for base, (walks, args, truth) in n_cases.items():
    n_gfa(base + '.gfa', n_nodes, walks)
    with open(base + '.truth', 'w') as f:
        f.write(truth)
    run('vg convert {0}.gfa > {0}.vg'.format(base))
    # (';', not '&&': ngap_shared patches nothing, so its report is empty and grep exits 1.  a panpatch
    #  failure still fails the case, as it leaves no .int for cat)
    run('panpatch {0}.vg {1} --bed {0}.int -f {0}.fa 2>/dev/null | grep -v "^chrom" | {2} > {0}.rep ; cat {0}.rep {0}.int > {0}.out'.format(base, args, telo_filter))
    run('diff {0}.out {0}.truth'.format(base))
    # and the patched assembly carries none of the donor's N's (ngap_shared keeps frag's own gap)
    if base != 'ngap_shared':
        run('! grep -v "^>" {0}.hap1.fa | grep -q N'.format(base))
    temp_files += ['{0}.gfa'.format(base), '{0}.vg'.format(base), '{0}.truth'.format(base), '{0}.out'.format(base),
                   '{0}.hap1.fa'.format(base)]



# A target contig that an accepted patch never threads through (an unplaced fragment binned to this
# chromosome) is written through as its own record, with a passthrough row and BED block.  Before, the
# patched chromosome was written and the fragment was silently dropped.  telopatch's layout (frag's q
# telomere is grafted from donor, with --max-subtelomere-graft off), plus frag's second contig ctg2 on a
# node of its own (17)
n_nodes[17] = n_rnd(5000)
with open('unthreaded_frag.gfa', 'w') as f:
    uf_walks = [('x', 0, 'chrX', [1, 2, 3, 4, 5]), ('frag', 1, 'ctg', [1, 2, 3, 6]), ('frag', 1, 'ctg2', [17]),
                ('donor', 1, 'ctg', [2, 3, 4, 5])]
    f.write('H\tVN:Z:1.1\n')
    for n in sorted(set(n for w in uf_walks for n in w[3])):
        f.write('S\t{}\t{}\n'.format(n, n_nodes[n]))
    for a, b in sorted(set(e for w in uf_walks for e in zip(w[3], w[3][1:]))):
        f.write('L\t{}\t+\t{}\t+\t0M\n'.format(a, b))
    for sample, hap, contig, w in uf_walks:
        f.write('W\t{}\t{}\t{}\t0\t{}\t{}\n'.format(sample, hap, contig, sum(len(n_nodes[n]) for n in w), ''.join('>{}'.format(n) for n in w)))
with open('unthreaded_frag.truth', 'w') as f:
    f.write(
        'chrX\t1\ttelomere\tfrag#1#ctg#0\t25004\tdonor#1#ctg#0\t32004\t3000\t0.0\t.\t.\taccepted\t.\t22004\t25004\n'
        'chrX\t1\tpassthrough\tfrag#1#ctg2#0\t5000\t.\t0\t.\t.\t.\t.\taccepted\tnot part of the patched contig; written as it is\t.\t.\n'
        '#Patched assembly on chrX for frag#1:\n'
        'frag#1#ctg#0\t0\t12004\t+\n'
        'donor#1#ctg#0\t10000\t32004\t+\n'
        '\n'
        '#Passthrough (unused by the patch) on chrX for frag#1:\n'
        'frag#1#ctg2#0\t0\t5000\t+\n'
        '\n')
run('vg convert unthreaded_frag.gfa > unthreaded_frag.vg')
run('panpatch unthreaded_frag.vg -r x -s frag -s donor --max-subtelomere-graft 0 --bed unthreaded_frag.int -f unthreaded_frag.fa 2>/dev/null | grep -v "^chrom" | {0} > unthreaded_frag.rep ; cat unthreaded_frag.rep unthreaded_frag.int > unthreaded_frag.out'.format(telo_filter))
run('diff unthreaded_frag.out unthreaded_frag.truth')
# and the fragment is in the FASTA, as it was input
run('grep -q "^>frag#1#ctg2" unthreaded_frag.hap1.fa')
temp_files += ['unthreaded_frag.gfa', 'unthreaded_frag.vg', 'unthreaded_frag.truth', 'unthreaded_frag.out', 'unthreaded_frag.hap1.fa']

# The interior-graft backstop judges what is left of a contig's foreign grafts by k-mer content once the
# partial-patch pass has excised the bad ones.  When what is left replaced only N's, there is nothing to
# judge, and it must not revert the contig.  Reference-free (frag is the reference):
#   frag  0 1 2 3 9 4 5 6 10 - an N gap (2), then 12kb of real sequence (4) and a short N run (5)
#   donor 0 1 7 3 9 8 6 10   - fills the gap with real sequence (7), and replaces 4 5 with unrelated sequence (8)
# frag is threaded on 0-1, 3-9 and 6-10, between two separate grafts.  The 8 graft shares none of 4's
# k-mers, so it is excised (12kb of real sequence counted toward the backstop's threshold); the 7 fill of
# the N gap passes flank anchoring and is kept.  Before, the backstop then found only N's left, scored them
# 0.0% and reverted the whole contig, the good fill with it.
bs_rng = random.Random(13)
bs_rnd = lambda n: ''.join(bs_rng.choice('ACGT') for _ in range(n))
bs_nodes = {0: bs_rnd(10000), 1: bs_rnd(10000), 2: 'N' * 1000, 3: bs_rnd(10000), 9: bs_rnd(10000), 4: bs_rnd(12000),
            5: 'N' * 100, 6: bs_rnd(10000), 10: bs_rnd(10000), 7: bs_rnd(1000), 8: bs_rnd(12100)}
bs_walks = {'frag' : [0, 1, 2, 3, 9, 4, 5, 6, 10], 'donor' : [0, 1, 7, 3, 9, 8, 6, 10]}
with open('backstop_nfill.gfa', 'w') as f:
    f.write('H\tVN:Z:1.1\n')
    for n in sorted(bs_nodes):
        f.write('S\t{}\t{}\n'.format(n, bs_nodes[n]))
    for a, b in sorted(set(e for w in bs_walks.values() for e in zip(w, w[1:]))):
        f.write('L\t{}\t+\t{}\t+\t0M\n'.format(a, b))
    for sample, w in bs_walks.items():
        f.write('W\t{}\t1\tctg\t0\t{}\t{}\n'.format(sample, sum(len(bs_nodes[n]) for n in w), ''.join('>{}'.format(n) for n in w)))
with open('backstop_nfill.truth', 'w') as f:
    f.write(
        'ctg\t1\tgap-fill\tfrag#1#ctg#0\t73100\tdonor#1#ctg#0\t73100\t12100\t0.0\t96.8\t100.0\trejected\tk-mer recovery 0.0% < 50.0% (repeat-region misjoin)\t41000\t53100\n'
        'ctg\t1\tgap-fill\tfrag#1#ctg#0\t73100\tdonor#1#ctg#0\t73100\t1000\t.\t100.0\t71.3\taccepted\t.\t20000\t21000\n'
        '#Patched assembly on ctg for frag#1:\n'
        'frag#1#ctg#0\t0\t10000\t+\n'
        'donor#1#ctg#0\t10000\t21000\t+\n'
        'frag#1#ctg#0\t21000\t73100\t+\n'
        '\n')
run('vg convert backstop_nfill.gfa > backstop_nfill.vg')
run('panpatch backstop_nfill.vg -r frag -s frag -s donor --bed backstop_nfill.int 2>/dev/null | grep -v "^chrom" | {0} > backstop_nfill.rep ; cat backstop_nfill.rep backstop_nfill.int > backstop_nfill.out'.format(telo_filter))
run('diff backstop_nfill.out backstop_nfill.truth')
temp_files += ['backstop_nfill.gfa', 'backstop_nfill.vg', 'backstop_nfill.truth', 'backstop_nfill.out']

# A patch that ends on a donor carries the donor's tail out to its end, replacing the target's tip beyond
# its outermost threaded piece: a tip graft, reported as a gap-fill row.  Nothing anchors its far end, so it
# is held to a telomere graft's limits (-M on the target bp it replaces, --max-subtelomere-graft on the donor bp
# it adds beyond that).  Reference-free:
#   frag  1 2 3 5    - an N gap (3) near its end, so the threading cannot get from 2 to 5 on frag
#   donor 1 2 4 5 6  - bridges it (4) and runs on past frag's end (6)
# frag 1-2, then donor 2-6 would replace frag's 3 and 5 (11kb) with 15kb: 4kb more than the 2kb default allows,
# so frag keeps its own end; with --max-subtelomere-graft 5000 the graft goes in
te_rng = random.Random(19)
te_rnd = lambda n: ''.join(te_rng.choice('ACGT') for _ in range(n))
te_nodes = {1: te_rnd(10000), 2: te_rnd(10000), 3: 'N' * 1000, 4: te_rnd(2000), 5: te_rnd(10000), 6: te_rnd(3000)}
te_walks = {'frag' : [1, 2, 3, 5], 'donor' : [1, 2, 4, 5, 6]}
with open('tip_extension.gfa', 'w') as f:
    f.write('H\tVN:Z:1.1\n')
    for n in sorted(te_nodes):
        f.write('S\t{}\t{}\n'.format(n, te_nodes[n]))
    for a, b in sorted(set(e for w in te_walks.values() for e in zip(w, w[1:]))):
        f.write('L\t{}\t+\t{}\t+\t0M\n'.format(a, b))
    for sample, w in te_walks.items():
        f.write('W\t{}\t1\tctg\t0\t{}\t{}\n'.format(sample, sum(len(te_nodes[n]) for n in w), ''.join('>{}'.format(n) for n in w)))
te_truth = {
    '' : ('ctg\t1\tgap-fill\tfrag#1#ctg#0\t31000\tdonor#1#ctg#0\t35000\t11000\t.\t.\t.\trejected\ttip graft adds 15000bp '
          "of donor sequence for 11000bp of the target's end, more than --max-subtelomere-graft 2000 beyond it; the target "
          'keeps its own end\t20000\t31000\n'
          '#Unpatched (input contig) on ctg for frag#1:\n'
          'frag#1#ctg#0\t0\t31000\t+\n'
          '\n'),
    '--max-subtelomere-graft 5000' : (
          'ctg\t1\tgap-fill\tfrag#1#ctg#0\t31000\tdonor#1#ctg#0\t35000\t11000\t.\t.\t.\taccepted\t.\t20000\t31000\n'
          '#Patched assembly on ctg for frag#1:\n'
          'frag#1#ctg#0\t0\t10000\t+\n'
          'donor#1#ctg#0\t10000\t35000\t+\n'
          '\n'),
    # -M 0 lets a tip graft replace nothing, as it does a telomere graft
    '--max-subtelomere-graft 5000 -M 0' : (
          'ctg\t1\tgap-fill\tfrag#1#ctg#0\t31000\tdonor#1#ctg#0\t35000\t11000\t.\t.\t.\trejected\ttip graft replaces 11000bp '
          "of the target's end, over --max-telomere-patch 0; the target keeps its own end\t20000\t31000\n"
          '#Unpatched (input contig) on ctg for frag#1:\n'
          'frag#1#ctg#0\t0\t31000\t+\n'
          '\n'),
    # a tip graft is gap filling, so --patch-types telomere leaves it out
    '--max-subtelomere-graft 5000 --patch-types telomere' : (
          'ctg\t1\tgap-fill\tfrag#1#ctg#0\t31000\tdonor#1#ctg#0\t35000\t11000\t.\t.\t.\trejected\ttip graft excluded by '
          "--patch-types (gap filling is off); the target keeps its own end\t20000\t31000\n"
          '#Unpatched (input contig) on ctg for frag#1:\n'
          'frag#1#ctg#0\t0\t31000\t+\n'
          '\n'),
}
run('vg convert tip_extension.gfa > tip_extension.vg')
for te_args, te_expect in te_truth.items():
    with open('tip_extension.truth', 'w') as f:
        f.write(te_expect)
    run('panpatch tip_extension.vg -r frag -s frag -s donor {1} --bed tip_extension.int 2>/dev/null | grep -v "^chrom" | {0} > tip_extension.rep ; cat tip_extension.rep tip_extension.int > tip_extension.out'.format(telo_filter, te_args))
    run('diff tip_extension.out tip_extension.truth')
temp_files += ['tip_extension.gfa', 'tip_extension.vg', 'tip_extension.truth', 'tip_extension.out']

# Interior grafts that the content and locus tests cannot vouch for are excised (the rest of the patch kept):
#   graft_growth - frag's 1kb N gap (2) is filled with 5kb of donor (7): accepted by default, but over
#                  --max-graft-size-change 1000 the +4kb change is too large to verify
#   graft_satellite - frag's 12kb tandem array (4, a 171bp unit, then a short N run 5 that sends the threading
#                  to donor) is replaced by donor's copy of the same array (8): k-mer recovery passes, as it
#                  would for any array of the family, but the run is low-complexity and cannot be verified
#   graft_satellite_window - the same array inside a larger swap: frag's 20kb either side of it (11, 12) differ
#                  from donor's (13, 14) by a SNP every 300bp, so the graft replaces 52kb whose distinct k-mer
#                  fraction is well over 0.5 as a whole.  measured in 10kb windows too (the default) the array
#                  still shows, and the graft is excised; with --graft-complexity-window 0 (the whole graft only)
#                  it is kept
gv_rng = random.Random(17)
gv_rnd = lambda n: ''.join(gv_rng.choice('ACGT') for _ in range(n))
gv_unit = gv_rnd(171)
gv_nodes = {0: gv_rnd(10000), 1: gv_rnd(10000), 2: 'N' * 1000, 3: gv_rnd(10000), 9: gv_rnd(10000), 10: gv_rnd(10000),
            7: gv_rnd(5000), 4: (gv_unit * 71)[:12000], 5: 'N' * 100, 8: (gv_unit * 72)[:12200],
            11: gv_rnd(20000), 12: gv_rnd(20000)}
gv_snp = lambda s: ''.join({'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[c] if i % 300 == 150 else c for i, c in enumerate(s))
gv_nodes[13], gv_nodes[14] = gv_snp(gv_nodes[11]), gv_snp(gv_nodes[12])
gv_win_walks = {'frag' : [0, 1, 3, 11, 4, 5, 12, 9, 10], 'donor' : [0, 1, 3, 13, 8, 14, 9, 10]}
gv_cases = {
    'graft_growth'    : ({'frag' : [0, 1, 2, 3, 9, 10], 'donor' : [0, 1, 7, 3, 9, 10]}, '', (
                          'ctg\t1\tgap-fill\tfrag#1#ctg#0\t51000\tdonor#1#ctg#0\t55000\t1000\t.\t100.0\t100.0\taccepted\t.\t20000\t21000\n'
                          '#Patched assembly on ctg for frag#1:\n'
                          'frag#1#ctg#0\t0\t10000\t+\n'
                          'donor#1#ctg#0\t10000\t25000\t+\n'
                          'frag#1#ctg#0\t21000\t51000\t+\n'
                          '\n')),
    # --max-swap-size-change does not apply to a gap fill: an N gap's length is only a placeholder
    'graft_growth_swapcap': ({'frag' : [0, 1, 2, 3, 9, 10], 'donor' : [0, 1, 7, 3, 9, 10]}, '--max-swap-size-change 1000', (
                          'ctg\t1\tgap-fill\tfrag#1#ctg#0\t51000\tdonor#1#ctg#0\t55000\t1000\t.\t100.0\t100.0\taccepted\t.\t20000\t21000\n'
                          '#Patched assembly on ctg for frag#1:\n'
                          'frag#1#ctg#0\t0\t10000\t+\n'
                          'donor#1#ctg#0\t10000\t25000\t+\n'
                          'frag#1#ctg#0\t21000\t51000\t+\n'
                          '\n')),
    'graft_growth_cap': ({'frag' : [0, 1, 2, 3, 9, 10], 'donor' : [0, 1, 7, 3, 9, 10]}, '--max-graft-size-change 1000', (
                          "ctg\t1\tgap-fill\tfrag#1#ctg#0\t51000\tdonor#1#ctg#0\t55000\t1000\t.\t100.0\t100.0\trejected\tchanges the contig's length by +4000bp, over --max-graft-size-change 1000 (unverifiable)\t20000\t21000\n"
                          '#Unpatched (input contig) on ctg for frag#1:\n'
                          'frag#1#ctg#0\t0\t51000\t+\n'
                          '\n')),
    'graft_satellite' : ({'frag' : [0, 1, 3, 4, 5, 9, 10], 'donor' : [0, 1, 3, 8, 9, 10]}, '', (
                          'ctg\t1\tgap-fill\tfrag#1#ctg#0\t62100\tdonor#1#ctg#0\t62200\t12100\t100.0\t100.0\t100.0\trejected\tlow-complexity sequence (distinct k-mer fraction 0.01 replaced/inserted < --min-graft-complexity 0.50) (unverifiable)\t30000\t42100\n'
                          '#Unpatched (input contig) on ctg for frag#1:\n'
                          'frag#1#ctg#0\t0\t62100\t+\n'
                          '\n')),
    'graft_satellite_window' : (gv_win_walks, '', (
                          'ctg\t1\tgap-fill\tfrag#1#ctg#0\t102100\tdonor#1#ctg#0\t102200\t52100\t89.6\t100.0\t100.0\trejected\tlow-complexity sequence (distinct k-mer fraction 0.02 in a 10000bp window replaced/inserted < --min-graft-complexity 0.50) (unverifiable)\t30000\t82100\n'
                          '#Unpatched (input contig) on ctg for frag#1:\n'
                          'frag#1#ctg#0\t0\t102100\t+\n'
                          '\n')),
    'graft_satellite_whole' : (gv_win_walks, '--graft-complexity-window 0', (
                          'ctg\t1\tgap-fill\tfrag#1#ctg#0\t102100\tdonor#1#ctg#0\t102200\t52100\t89.6\t100.0\t100.0\taccepted\t.\t30000\t82100\n'
                          '#Patched assembly on ctg for frag#1:\n'
                          'frag#1#ctg#0\t0\t20000\t+\n'
                          'donor#1#ctg#0\t20000\t82200\t+\n'
                          'frag#1#ctg#0\t82100\t102100\t+\n'
                          '\n')),
}
for base, (walks, args, truth) in gv_cases.items():
    with open(base + '.gfa', 'w') as f:
        f.write('H\tVN:Z:1.1\n')
        for n in sorted(set(n for w in walks.values() for n in w)):
            f.write('S\t{}\t{}\n'.format(n, gv_nodes[n]))
        for a, b in sorted(set(e for w in walks.values() for e in zip(w, w[1:]))):
            f.write('L\t{}\t+\t{}\t+\t0M\n'.format(a, b))
        for sample, w in walks.items():
            f.write('W\t{}\t1\tctg\t0\t{}\t{}\n'.format(sample, sum(len(gv_nodes[n]) for n in w), ''.join('>{}'.format(n) for n in w)))
    with open(base + '.truth', 'w') as f:
        f.write(truth)
    run('vg convert {0}.gfa > {0}.vg'.format(base))
    run('panpatch {0}.vg -r frag -s frag -s donor {1} --bed {0}.int 2>/dev/null | grep -v "^chrom" | {2} > {0}.rep ; cat {0}.rep {0}.int > {0}.out'.format(base, args, telo_filter))
    run('diff {0}.out {0}.truth'.format(base))
    temp_files += ['{0}.gfa'.format(base), '{0}.vg'.format(base), '{0}.truth'.format(base), '{0}.out'.format(base)]

# --error-bed: assembly-error intervals block the steps of their own path, not the node, so where the target's
# flagged sequence is aligned to a donor's clean copy, the donor can thread across it and the patch differs
# from the target only where the donor does.  Reference-free, every node 10kb unless noted:
#   frag  10 0 1 2 3 5 6 8 9 - 3 is frag's 1kb allele of a bubble
#   donor 10 0 1 2 4 5 6 8 9 - 4 is the donor's; the donor is ctg[1000-82000], so its BED is shifted by 1000
#   errbed_swap     - frag's 2 3 5 flagged: the donor's 4 replaces frag's 3, and only that (the target takes
#                     over again on 5, the last shared node of its flagged region)
#   errbed_same     - frag's 1 2 flagged, where the donor is the same; the 3/4 bubble is past the region, so
#                     it stays frag's 3: the donor run spells exactly what it replaces, nothing to patch
#   errbed_donor    - as errbed_swap, but the donor's 4 is flagged too: frag keeps its own sequence
#   errbed_collapse - frag 10 0 1 2 5 6 8 9 with 2 5 flagged, and the donor has 2kb more (7) between them: the
#                     donor's insertion goes in (a collapse fix, which N-masking the target's unaligned
#                     sequence could not make: frag has none there)
#   errbed_edge     - as errbed_collapse, but only frag's 5 flagged: the insertion right where the region
#                     starts goes in, replacing nothing
eb_rng = random.Random(23)
eb_rnd = lambda n: ''.join(eb_rng.choice('ACGT') for _ in range(n))
eb_nodes = {0: eb_rnd(10000), 1: eb_rnd(10000), 2: eb_rnd(10000), 3: eb_rnd(1000), 4: eb_rnd(1000),
            5: eb_rnd(10000), 6: eb_rnd(10000), 7: eb_rnd(2000), 8: eb_rnd(10000), 9: eb_rnd(10000), 10: eb_rnd(10000)}
eb_swap_walks = {'frag': [10, 0, 1, 2, 3, 5, 6, 8, 9], 'donor': [10, 0, 1, 2, 4, 5, 6, 8, 9]}
eb_ins_walks = {'frag': [10, 0, 1, 2, 5, 6, 8, 9], 'donor': [10, 0, 1, 2, 7, 5, 6, 8, 9]}
eb_patched = lambda d: 'donor#1#ctg#0[1000-{}]'.format(1000 + d)
eb_cases = {
    'errbed_swap'     : (eb_swap_walks, 'frag#1#ctg\t30000\t51000\n', [10, 0, 1, 2, 4, 5, 6, 8, 9],
                         'ctg\t1\tgap-fill\tfrag#1#ctg#0\t81000\t' + eb_patched(81000) + '\t81000\t11000\t91.0\t100.0\t100.0\taccepted\t.\t30000\t41000\n'
                         '#Patched assembly on ctg for frag#1:\n'
                         'frag#1#ctg#0\t0\t20000\t+\n' + eb_patched(81000) + '\t20000\t41000\t+\n'
                         'frag#1#ctg#0\t41000\t81000\t+\n\n'),
    'errbed_same'     : (eb_swap_walks, 'frag#1#ctg\t20000\t40000\n', None,
                         'ctg\t1\tgap-fill\tfrag#1#ctg#0\t81000\t' + eb_patched(81000) + '\t81000\t10000\t100.0\t100.0\t97.6\trejected\t'
                         "donor sequence is identical to the target's: nothing to patch\t20000\t30000\n"
                         '#Unpatched (input contig) on ctg for frag#1:\n'
                         'frag#1#ctg#0\t0\t81000\t+\n\n'),
    'errbed_donor'    : (eb_swap_walks, 'frag#1#ctg\t30000\t51000\ndonor#1#ctg\t41000\t42000\n', None,
                         '#Unpatched (input contig) on ctg for frag#1:\n'
                         'frag#1#ctg#0\t0\t81000\t+\n\n'),
    # the donor's interval named with its phase block (donor#1#ctg#0), as the path is less its [1000-82000]
    'errbed_donor_phase' : (eb_swap_walks, 'frag#1#ctg\t30000\t51000\ndonor#1#ctg#0\t41000\t42000\n', None,
                         '#Unpatched (input contig) on ctg for frag#1:\n'
                         'frag#1#ctg#0\t0\t81000\t+\n\n'),
    'errbed_collapse' : (eb_ins_walks, 'frag#1#ctg\t30000\t50000\n', [10, 0, 1, 2, 7, 5, 6, 8, 9],
                         'ctg\t1\tgap-fill\tfrag#1#ctg#0\t80000\t' + eb_patched(82000) + '\t82000\t10000\t100.0\t100.0\t100.0\taccepted\t.\t30000\t40000\n'
                         '#Patched assembly on ctg for frag#1:\n'
                         'frag#1#ctg#0\t0\t20000\t+\n' + eb_patched(82000) + '\t20000\t42000\t+\n'
                         'frag#1#ctg#0\t40000\t80000\t+\n\n'),
    'errbed_edge'     : (eb_ins_walks, 'frag#1#ctg\t40000\t50000\n', [10, 0, 1, 2, 7, 5, 6, 8, 9],
                         'ctg\t1\tgap-fill\tfrag#1#ctg#0\t80000\t' + eb_patched(82000) + '\t82000\t0\t.\t100.0\t100.0\taccepted\t.\t40000\t40000\n'
                         '#Patched assembly on ctg for frag#1:\n'
                         'frag#1#ctg#0\t0\t30000\t+\n' + eb_patched(82000) + '\t30000\t42000\t+\n'
                         'frag#1#ctg#0\t40000\t80000\t+\n\n'),
}
for base, (walks, errors, patched_walk, truth) in eb_cases.items():
    with open(base + '.gfa', 'w') as f:
        f.write('H\tVN:Z:1.1\n')
        for n in sorted(eb_nodes):
            if any(n in w for w in walks.values()):
                f.write('S\t{}\t{}\n'.format(n, eb_nodes[n]))
        for a, b in sorted(set(e for w in walks.values() for e in zip(w, w[1:]))):
            f.write('L\t{}\t+\t{}\t+\t0M\n'.format(a, b))
        for sample, w in walks.items():
            st = 1000 if sample == 'donor' else 0
            f.write('W\t{}\t1\tctg\t{}\t{}\t{}\n'.format(sample, st, st + sum(len(eb_nodes[n]) for n in w),
                                                       ''.join('>{}'.format(n) for n in w)))
    with open(base + '.errors.bed', 'w') as f:
        f.write(errors)
    with open(base + '.truth', 'w') as f:
        f.write(truth)
    run('vg convert {0}.gfa > {0}.vg'.format(base))
    run('panpatch {0}.vg -r frag -s frag -s donor --error-bed {0}.errors.bed --bed {0}.int -f {0}.fa 2>/dev/null | grep -v "^chrom" | {1} > {0}.rep ; cat {0}.rep {0}.int > {0}.out'.format(base, telo_filter))
    run('diff {0}.out {0}.truth'.format(base))
    if patched_walk:
        # the patched contig spells the donor's walk across the swap, the target's elsewhere
        with open(base + '.hap1.fa') as f:
            got = ''.join(l.strip() for l in f if not l.startswith('>'))
        if got != ''.join(eb_nodes[n] for n in patched_walk):
            print('failed: {} patched sequence is not {}'.format(base, patched_walk))
            fail_count += 1
    temp_files += ['{0}.{1}'.format(base, x) for x in ['gfa', 'vg', 'errors.bed', 'truth', 'out', 'hap1.fa']]

# A donor that visits a node twice (as through a satellite array) is followed on only from the visit the
# threading reached it on.  frag 1 2 3 4 5 6 (10kb each) with 3 4 5 flagged; donor 1 4 5 6 2 3 4 N 8.  The
# threading leaves frag at 2 for the donor's second half (2 3 4), which runs into its N after 4.  Its first
# visit of 4 does go on to 5, and taking it put donor pieces [2 3 4] and [4 5] side by side, which were merged
# into one interval running backwards along the donor (a crash in intervals_to_sequence).  Now frag carries
# on across its own flagged sequence, and the donor's 3, shared with frag, changes nothing
rv_rng = random.Random(29)
rv_rnd = lambda n: ''.join(rv_rng.choice('ACGT') for _ in range(n))
rv_nodes = {1: rv_rnd(10000), 2: rv_rnd(10000), 3: rv_rnd(10000), 4: rv_rnd(10000), 5: rv_rnd(10000), 6: rv_rnd(10000),
            7: 'N' * 1000, 8: rv_rnd(10000)}
rv_walks = {'frag': [1, 2, 3, 4, 5, 6], 'donor': [1, 4, 5, 6, 2, 3, 4, 7, 8]}
with open('donor_revisit.gfa', 'w') as f:
    f.write('H\tVN:Z:1.1\n')
    for n in sorted(rv_nodes):
        f.write('S\t{}\t{}\n'.format(n, rv_nodes[n]))
    for a, b in sorted(set(e for w in rv_walks.values() for e in zip(w, w[1:]))):
        f.write('L\t{}\t+\t{}\t+\t0M\n'.format(a, b))
    for sample, w in rv_walks.items():
        f.write('W\t{}\t1\tctg\t0\t{}\t{}\n'.format(sample, sum(len(rv_nodes[n]) for n in w), ''.join('>{}'.format(n) for n in w)))
with open('donor_revisit.errors.bed', 'w') as f:
    f.write('frag#1#ctg\t20000\t50000\n')
with open('donor_revisit.truth', 'w') as f:
    f.write("ctg\t1\tgap-fill\tfrag#1#ctg#0\t60000\tdonor#1#ctg#0\t81000\t10000\t100.0\t100.0\t100.0\trejected\t"
            "donor sequence is identical to the target's: nothing to patch\t20000\t30000\n"
            '#Unpatched (input contig) on ctg for frag#1:\n'
            'frag#1#ctg#0\t0\t60000\t+\n\n')
run('vg convert donor_revisit.gfa > donor_revisit.vg')
run('panpatch donor_revisit.vg -r frag -s frag -s donor --error-bed donor_revisit.errors.bed --bed donor_revisit.int 2>/dev/null | grep -v "^chrom" | {0} > donor_revisit.rep ; cat donor_revisit.rep donor_revisit.int > donor_revisit.out'.format(telo_filter))
run('diff donor_revisit.out donor_revisit.truth')
temp_files += ['donor_revisit.{}'.format(x) for x in ['gfa', 'vg', 'errors.bed', 'truth', 'out']]

# --telomere-report: the standalone FASTA scan must agree exactly with the in-graph measurement, since
# both go through print_contig_telomere_line.  Patch a graph with only the target sample, so the output
# FASTA is that one path verbatim and the two measurements describe identical sequence -- then the
# "#Contig" line panpatch printed from the graph and the one it prints from the FASTA must be the same.
# (A caller like cactus-panpatch relies on this: it measures its post-processed output with the FASTA
# mode and compares against runs that were measured in-graph.)
run('vg convert telopatch.gfa > telorep.vg')
run('panpatch telorep.vg -r x -s frag -f telorep.fa 2>/dev/null | grep "^#Contig" > telorep.graph.txt')
run('panpatch --telomere-report telorep.hap1.fa > telorep.fasta.txt')
run('diff telorep.graph.txt telorep.fasta.txt')
# and the scan must not depend on repeat-masking: soft-masked input gives the same answer
run("awk '/^>/ {print; next} {print tolower($0)}' telorep.hap1.fa > telorep.lower.fa")
run('panpatch --telomere-report telorep.lower.fa > telorep.lower.txt')
run('diff telorep.fasta.txt telorep.lower.txt')
temp_files += ['telorep.vg', 'telorep.hap1.fa', 'telorep.lower.fa',
               'telorep.graph.txt', 'telorep.fasta.txt', 'telorep.lower.txt']

# Regressions from a review of the threading, graft and telomere code (Nov 2026).  Each graph is built from
# signed walks; the output FASTA is decoded back into signed node ids and compared with what should be
# written, and the report must carry the reason given.  Before the fixes, each of these wrote a false patch
# (deleted, duplicated, inverted or appended sequence, a lost telomere), a misleading report, or aborted.
def rg_rnd(seed, n):
    rng = random.Random(seed)
    return ''.join(rng.choice('ACGT') for _ in range(n))

def rg_rc(s):
    return s[::-1].translate(str.maketrans('ACGTN', 'TGCAN'))

def rg_case(base, nodes, walks, args, errors, expect, reasons=()):
    """walks: [(sample, hap, contig, ['+1', '-2', ...])]; expect: {output record: 'signed node ids'}"""
    global fail_count
    pw = lambda w: [(int(x[1:]), x[0] == '-') for x in w]
    with open(base + '.gfa', 'w') as f:
        f.write('H\tVN:Z:1.1\n')
        for n in sorted(set(n for _, _, _, w in walks for n, _ in pw(w))):
            f.write('S\t{}\t{}\n'.format(n, nodes[n]))
        edges = set()
        for _, _, _, w in walks:
            for (a, ar), (b, br) in zip(pw(w), pw(w)[1:]):
                edges.add((a, '-' if ar else '+', b, '-' if br else '+'))
        for e in sorted(edges):
            f.write('L\t{}\t{}\t{}\t{}\t0M\n'.format(*e))
        for sample, hap, contig, w in walks:
            f.write('W\t{}\t{}\t{}\t0\t{}\t{}\n'.format(sample, hap, contig, sum(len(nodes[n]) for n, _ in pw(w)),
                                                      ''.join(('<' if r else '>') + str(n) for n, r in pw(w))))
    if errors:
        with open(base + '.errors.bed', 'w') as f:
            f.write(errors)
        args += ' --error-bed {}.errors.bed'.format(base)
    run('vg convert {0}.gfa > {0}.vg'.format(base))
    if os.path.exists(base + '.hap1.fa'):
        os.remove(base + '.hap1.fa')
    p = subprocess.run('panpatch {0}.vg {1} -f {0}.fa 2>/dev/null'.format(base, args), shell=True, capture_output=True, text=True)
    recs, cur = {}, None
    if os.path.exists(base + '.hap1.fa'):
        for l in open(base + '.hap1.fa'):
            l = l.strip()
            if l.startswith('>'):
                cur = l[1:]; recs[cur] = ''
            elif cur is not None:
                recs[cur] += l
    by_len = sorted(nodes.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    got = {}
    for name, seq in recs.items():
        ids, i = [], 0
        while i < len(seq):
            for n, s in by_len:
                if seq.startswith(s, i):
                    ids.append('+%d' % n); i += len(s); break
                if seq.startswith(rg_rc(s), i):
                    ids.append('-%d' % n); i += len(s); break
            else:
                ids.append('?'); break
        got[name] = ' '.join(ids)
    if p.returncode != 0 or got != expect or not all(r in p.stdout for r in reasons):
        print('failed: {} (exit {}) wrote {}, expected {}'.format(base, p.returncode, got, expect))
        for r in reasons:
            if r not in p.stdout:
                print('  report lacks: {}'.format(r))
        fail_count += 1
    temp_files.extend([base + x for x in ['.gfa', '.vg', '.errors.bed', '.hap1.fa']])

rg_K = lambda seeds, n: {k: rg_rnd(s, n) for k, s in seeds.items()}
rg_fwd = lambda ids: ['+%d' % k for k in ids]

# strand: the donor carries an inversion of 3-5, and frag's 3 is flagged.  The donor was searched both ways
# and followed into its inverted copy, deleting frag's 3 4 with nothing put back; frag now keeps its own
rg_n = rg_K({1: 101, 2: 102, 3: 103, 4: 104, 5: 105, 6: 106}, 1000); rg_n[1] = rg_rnd(107, 50000); rg_n[6] = rg_rnd(108, 50000)
rg_case('rg_inversion', rg_n, [('frag', 1, 'ctg', rg_fwd(range(1, 7))), ('donor', 1, 'ctg', ['+1', '+2', '-5', '-4', '-3', '+6'])],
        '-r frag -s frag -s donor', 'frag#1#ctg\t51000\t52000\n', {'frag#1#ctg#0': '+1 +2 +3 +4 +5 +6'})
# strand: the donor is a rotation of the target (3 4 1 2): it was followed backwards across frag's gap,
# writing an inverted chimera with 50kb repeated.  frag keeps its gap
rg_n = {1: rg_rnd(201, 50000), 2: rg_rnd(202, 1000), 3: rg_rnd(203, 1000), 4: rg_rnd(204, 50000), 7: 'N' * 1000}
rg_case('rg_rotation', rg_n, [('frag', 1, 'ctg', ['+1', '+2', '+7', '+3', '+4']), ('donor', 1, 'ctg', ['+3', '+4', '+1', '+2'])],
        '-r frag -s frag -s donor', None, {'frag#1#ctg#0': '+1 +2 +7 +3 +4'})
# a donor visiting the node it is taken up on twice (a stray copy of 5 at its start): its first visit was
# taken, and the run went back over 0-4, turning a 1kb gap fill into a 51kb duplication
rg_n = rg_K({k: 300 + k for k in range(9)}, 10000); rg_n[5] = rg_rnd(399, 1000); rg_n[10] = 'N' * 1000
rg_case('rg_stray_copy', rg_n, [('frag', 1, 'ctg', rg_fwd([0, 1, 2, 3, 4, 5, 10, 6, 7, 8])), ('donor', 1, 'ctg', rg_fwd([5, 0, 1, 2, 3, 4, 5, 6, 7, 8]))],
        '-r frag -s frag -s donor', None, {'frag#1#ctg#0': '+0 +1 +2 +3 +4 +5 +10 +6 +7 +8'})
# a donor run back over anchors already patched (donor 3 1 2 4 5 across frag's gap after 3): 1 2 written twice
rg_n = rg_K({k: 400 + k for k in range(6)}, 10000); rg_n[9] = 'N' * 1000
rg_case('rg_backtrack', rg_n, [('frag', 1, 'ctg', rg_fwd([0, 1, 2, 3, 9, 4, 5])), ('donor', 1, 'ctg', rg_fwd([3, 1, 2, 4, 5]))],
        '-r frag -s frag -s donor', None, {'frag#1#ctg#0': '+0 +1 +2 +3 +9 +4 +5'})
# tip graft: frag's last 5kb flagged and the donor running 300kb past frag's end: the whole tail was appended
rg_n = rg_K({1: 501, 2: 502, 3: 503}, 10000); rg_n.update(rg_K({4: 504, 5: 505, 6: 506}, 100000)); rg_n[0] = rg_rnd(500, 200000)
rg_case('rg_tip_tail', rg_n, [('frag', 1, 'ctg', rg_fwd([0, 1, 2, 3])), ('donor', 1, 'ctg', rg_fwd(range(7)))],
        '-r frag -s frag -s donor', 'frag#1#ctg\t225000\t230000\n', {'frag#1#ctg#0': '+0 +1 +2 +3'},
        ['tip graft adds 310000bp of donor sequence for 10000bp'])
# an identical donor across a flagged region of a target stored antisense to the reference: compared strand
# against strand, it passed as a change (and the contig came out reverse-complemented)
rg_n = rg_K({k: 600 + k for k in range(1, 9)}, 10000)
rg_case('rg_same_reverse', rg_n, [('x', 0, 'chrX', rg_fwd(range(1, 9))), ('frag', 1, 'ctg', ['-%d' % k for k in range(8, 0, -1)]),
                                  ('donor', 1, 'ctg', rg_fwd(range(1, 9)))],
        '-r x -s frag -s donor', 'frag#1#ctg\t30000\t50000\n', {'frag#1#ctg#0': '-8 -7 -6 -5 -4 -3 -2 -1'},
        ["donor sequence is identical to the target's: nothing to patch"])
# a target contig (ctgB) flagged end to end, covered by the donor: it was written twice, once as the donor's
# copy in the patched record and once as itself
rg_n = rg_K({k: 700 + k for k in range(1, 9)}, 10000)
rg_case('rg_flagged_contig', rg_n, [('x', 0, 'chrX', rg_fwd([1, 2, 3, 5, 6, 7, 8])), ('frag', 1, 'ctgA', rg_fwd([1, 2, 3, 5])),
                                    ('frag', 1, 'ctgB', rg_fwd([6, 7, 8])), ('donor', 1, 'ctg', rg_fwd([1, 2, 3, 5, 6, 7, 8]))],
        '-r x -s frag -s donor', 'frag#1#ctgB\t0\t30000\n', {'frag#1#ctgA#0': '+1 +2 +3 +5', 'frag#1#ctgB#0': '+6 +7 +8'})
# frag carries its own copy of 5 6 at its start; after the donor bridges its gap, frag was taken up again at
# that first copy and 60kb was written twice.  Now only the gap is filled
rg_n = rg_K({k: 800 + k for k in range(1, 9)}, 10000); rg_n[9] = 'N' * 1000
rg_case('rg_target_early', rg_n, [('x', 0, 'chr', rg_fwd(range(1, 9))), ('frag', 1, 'ctg', rg_fwd([5, 6, 1, 2, 3, 4, 9, 5, 6, 7, 8])),
                                  ('donor', 1, 'ctg', rg_fwd(range(1, 9)))],
        '-r x -s frag -s donor', None, {'chr_hap_1': '+5 +6 +1 +2 +3 +4 +5 +6 +7 +8'})
# frag's capped 5' end (a 6kb telomere) would be trimmed off by the donor's run from its own end: the
# telomere guard only looked at trims of 20kb or more
rg_n = {1: 'CCCTAA' * 1000, 2: rg_rnd(902, 3000), 3: rg_rnd(903, 3000), 4: rg_rnd(904, 3000), 5: rg_rnd(905, 100000),
        6: rg_rnd(906, 100000), 7: rg_rnd(907, 100000), 8: 'TTAGGG' * 1001}   # (not 1000: the decoder would read it as -1)
rg_case('rg_telomere_trim', rg_n, [('x', 0, 'chrX', rg_fwd([2, 5, 6, 7, 8])), ('frag', 1, 'ctg', rg_fwd([1, 2, 3, 5, 6, 7, 8])),
                                   ('donor', 1, 'ctg', rg_fwd([2, 4, 5, 6, 7, 8]))],
        '-r x -s frag -s donor', 'frag#1#ctg\t9000\t12000\n', {'frag#1#ctg#0': '+1 +2 +3 +5 +6 +7 +8'},
        ["tip graft would take off the telomere at the target's end"])
# a scaffold join through a donor bridge that carries a telomere array (9): with -T it was accepted, the
# internal-telomere check only looking at the interior's overall density
rg_n = {20: 'CCCTAA' * 1400, 1: rg_rnd(1001, 60000), 2: rg_rnd(1002, 60000), 9: 'TTAGGG' * 1500, 5: rg_rnd(1005, 2000),
        3: rg_rnd(1003, 60000), 4: rg_rnd(1004, 60000), 21: 'TTAGGG' * 1600}
rg_case('rg_join_telomere', rg_n, [('x', 0, 'chrX', rg_fwd([20, 1, 2, 9, 5, 3, 4, 21])), ('frag', 1, 'ctgA', rg_fwd([20, 1, 2, 9])),
                                   ('frag', 1, 'ctgB', rg_fwd([3, 4, 21])), ('donor', 1, 'ctg', rg_fwd([20, 1, 2, 9, 5, 3, 4, 21]))],
        '-r x -s frag -s donor -T', None, {'frag#1#ctgA#0': '+20 +1 +2 +9', 'frag#1#ctgB#0': '+3 +4 +21'},
        ['internal telomere in patched contig (scaffold misjoin)'])
# a 25kb contig whose front and back telomere handoffs cross: the interval ran backwards
rg_n = rg_K({1: 1101, 2: 1102, 3: 1103, 4: 1104, 5: 1105}, 5000); rg_n.update({10: 'CCCTAA' * 1000, 11: rg_rnd(1111, 500), 12: rg_rnd(1112, 500), 13: 'TTAGGG' * 1000})
rg_case('rg_cross_handoff', rg_n, [('x', 0, 'chrX', rg_fwd(range(1, 6))), ('frag', 1, 'ctg', rg_fwd(range(1, 6))),
                                   ('donor', 1, 'ctg', rg_fwd([10, 11, 4])), ('donor2', 1, 'ctg', rg_fwd([2, 12, 13]))],
        '-r x -s frag -s donor -s donor2', None, {'frag#1#ctg#0': '+1 +2 +3 +4 +5'},
        ["the contig is too short for both ends' telomere handoffs (they cross)"])
# a scaffold bridge across two donors, the second sharing only the junction node with ctgB: both flanks were
# measured against the first donor, and it passed at 100%/100%
rg_n = rg_K({k: 1200 + k for k in range(1, 14)}, 10000); rg_n[10] = 'N' * 1000; rg_n[11] = 'N' * 1001
rg_case('rg_two_donor_bridge', rg_n, [('x', 0, 'chrX', rg_fwd(range(1, 10))), ('frag', 1, 'ctgA', rg_fwd([1, 2, 3, 10, 5])),
                                      ('frag', 1, 'ctgB', rg_fwd([7, 8, 9])), ('d1', 1, 'ctg', rg_fwd([1, 2, 3, 4, 5, 11, 7, 8, 9])),
                                      ('d2', 1, 'ctg', rg_fwd([5, 6, 7, 12, 13]))],
        '-r x -s frag -s d1 -s d2', None, {'frag#1#ctgA#0': '+1 +2 +3 +10 +5', 'frag#1#ctgB#0': '+7 +8 +9'},
        ['foreign bridge anchoring 100.0%/0.0% < 50.0% (wrong-locus join)'])
# --min-cover counted a 40kb fragment the patch never uses, and reverted a good telomere patch
rg_n = {1: 'CCCTAA' * 334, 2: rg_rnd(1302, 10000), 3: rg_rnd(1303, 10000), 4: rg_rnd(1304, 10000), 5: 'TTAGGG' * 335,
        6: rg_rnd(1306, 3000), 17: rg_rnd(1317, 40000)}
rg_case('rg_min_cover', rg_n, [('x', 0, 'chrX', rg_fwd([1, 2, 3, 4, 5])), ('frag', 1, 'ctg', rg_fwd([1, 2, 3, 6])),
                               ('frag', 1, 'ctg2', rg_fwd([17])), ('donor', 1, 'ctg', rg_fwd([2, 3, 4, 5]))],
        '-r x -s frag -s donor --max-subtelomere-graft 0', None, {'chrX_hap_1': '+1 +2 +3 +4 +5', 'frag#1#ctg2#0': '+17'})
# the reference's first node traversed in reverse: an assert aborted the whole run
rg_n = rg_K({1: 1401, 2: 1402, 3: 1403}, 10000); rg_n[9] = 'N' * 100
rg_case('rg_reverse_first', rg_n, [('frag', 1, 'ctg', ['-1', '+2', '+9', '+3']), ('donor', 1, 'ctg', ['-1', '+2', '+3'])],
        '-r frag -s frag -s donor', None, {'ctg_hap_1': '-1 +2 +3'})
# a graft spliced from two donor contigs (frag's flagged 3 4: ctgA has 3 then runs off, ctgB goes on from 3 through
# 9 to 5): sequence no donor has in one piece, excised by default and kept with --max-graft-pieces 0
rg_n = rg_K({1: 1501, 2: 1502, 5: 1505, 6: 1506, 7: 1507}, 10000); rg_n.update(rg_K({3: 1503, 4: 1504, 8: 1508, 9: 1509}, 1000))
rg_w = [('frag', 1, 'ctg', rg_fwd(range(1, 8))), ('donor', 1, 'ctgA', rg_fwd([1, 2, 3, 8])), ('donor', 1, 'ctgB', rg_fwd([3, 9, 5, 6, 7]))]
rg_case('rg_two_piece_graft', rg_n, rg_w, '-r frag -s frag -s donor', 'frag#1#ctg\t20000\t22000\n',
        {'frag#1#ctg#0': '+1 +2 +3 +4 +5 +6 +7'},
        ['spliced from 2 pieces of donor sequence, over --max-graft-pieces 1 (no donor has it in one piece) (unverifiable)'])
rg_case('rg_two_piece_graft_off', rg_n, rg_w, '-r frag -s frag -s donor --max-graft-pieces 0', 'frag#1#ctg\t20000\t22000\n',
        {'ctg_hap_1': '+1 +2 +3 +9 +5 +6 +7'})
# a patched contig keeps a gap (21) that no donor crosses: it gets a row of its own, not silence
rg_n = rg_K({1: 1601, 2: 1602, 5: 1605, 6: 1606, 7: 1607}, 10000); rg_n.update(rg_K({3: 1603, 4: 1604}, 30000))
rg_n[20] = 'N' * 1000; rg_n[21] = 'N' * 2000
rg_case('rg_unfilled_gap', rg_n, [('frag', 1, 'ctg', rg_fwd([1, 2, 20, 3, 4, 21, 5, 6])), ('donor', 1, 'ctg', rg_fwd([1, 2, 7, 3, 4]))],
        '-r frag -s frag -s donor', None, {'ctg_hap_1': '+1 +2 +7 +3 +4 +21 +5 +6'},
        ['ctg\t1\tgap-fill\tfrag#1#ctg#0\t103000\t.\t0\t.\t.\t.\t.\trejected\tno donor could be threaded across this gap; kept as it is\t81000\t83000'])
# --max-swap-size-change: the donor's 2kb insertion into frag's flagged region (errbed_collapse, with no N gap
# to fill) is a swap that changes the length, so with a 1kb cap frag keeps its own sequence; the 10kb default
# keeps the insertion
rg_ins = [('frag', 1, 'ctg', rg_fwd([10, 0, 1, 2, 5, 6, 8, 9])), ('donor', 1, 'ctg', rg_fwd([10, 0, 1, 2, 7, 5, 6, 8, 9]))]
rg_case('rg_swap_cap', eb_nodes, rg_ins, '-r frag -s frag -s donor --max-swap-size-change 1000', 'frag#1#ctg\t30000\t50000\n',
        {'frag#1#ctg#0': '+10 +0 +1 +2 +5 +6 +8 +9'}, ['swap changes the contig\'s length by +2000bp, over --max-swap-size-change 1000 (unverifiable)'])
rg_case('rg_swap_cap_default', eb_nodes, rg_ins, '-r frag -s frag -s donor', 'frag#1#ctg\t30000\t50000\n',
        {'ctg_hap_1': '+10 +0 +1 +2 +7 +5 +6 +8 +9'})
# -e naming a sample with two contigs here: a failed patch aborted the run (assert) ...
run('vg convert interior_splice.gfa > interior_splice.e.vg')
run('panpatch interior_splice.e.vg -r x -s frag -e frag -f interior_splice.e.fa 2>/dev/null > interior_splice.e.rep && test $(grep -c "^>" interior_splice.e.hap1.fa) -eq 2')
# ... and a telomere revert wrote only the first of them, dropping ctgB
run('vg convert scaffold.gfa > scaffold.e.vg')
run('panpatch scaffold.e.vg -r x -s frag -T -e frag -f scaffold.e.fa 2>/dev/null > scaffold.e.rep && test $(grep -c "^>" scaffold.e.hap1.fa) -eq 2')
temp_files += ['interior_splice.e.vg', 'interior_splice.e.hap1.fa', 'scaffold.e.vg', 'scaffold.e.hap1.fa']

temp_files += [f for f in os.listdir('.') if f.endswith('.int') or f.endswith('.rep')]
for f in temp_files:
    if os.path.isfile(f):
        os.remove(f)

if fail_count:
    print('{} tests failed'.format(fail_count))
else:
    print('all tests passed')
    
sys.exit(fail_count)
