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
telo_filter = 'grep -v "^#Contig"'
for base in ['telopatch', 'telopatch_rev', 'telopatch_front', 'telopatch_multi', 'telopatch_buried']:
    run('vg convert {0}.gfa > {0}.vg'.format(base))
    run('panpatch {0}.vg -r x -s frag -s donor -T --bed {0}.int -f {0}.fa 2>/dev/null | {1} > {0}.rep && cat {0}.rep {0}.int > {0}.patch.bed'.format(base, telo_filter))
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
    'ntelo_skip'    : (n_telo, '-r x -s frag -s ndonor -s donor', (
                       'chrX\t1\ttelomere\tfrag#1#ctg#0\t25004\tdonor#1#ctg#0\t32004\t3000\t0.0\t.\t.\taccepted\t.\t22004\t25004\n'
                       '#Patched assembly on chrX for frag#1:\n'
                       'frag#1#ctg#0\t0\t12004\t+\n'
                       'donor#1#ctg#0\t10000\t32004\t+\n'
                       '\n')),
    'ntelo_only'    : (n_telo, '-r x -s frag -s ndonor', (
                       "chrX\t1\ttelomere\tfrag#1#ctg#0\t25004\tndonor#1#ctg#0\t33004\t3000\t.\t.\t.\trejected\tdonor's run to its telomere contains N (an assembly gap or masked error region)\t22004\t25004\n"
                       '#Patched assembly on chrX for frag#1:\n'
                       'frag#1#ctg#0\t0\t25004\t+\n'
                       '\n')),
    'ntelo_handoff' : (dict(n_telo, frag=[1, 2, 3, 7]), '-r x -s frag -s ndonor -s donor', (
                       'chrX\t1\ttelomere\tfrag#1#ctg#0\t23004\tdonor#1#ctg#0\t32004\t1000\t0.0\t.\t.\taccepted\t.\t22004\t23004\n'
                       '#Patched assembly on chrX for frag#1:\n'
                       'frag#1#ctg#0\t0\t12004\t+\n'
                       'donor#1#ctg#0\t10000\t32004\t+\n'
                       '\n')),
    'ngap_shared'   : ({'frag' : [11, 12, 13, 14, 16], 'donor' : [11, 12, 13, 15, 16]}, '-r frag -s frag -s donor', (
                       '#Patched assembly on ctg for frag#1:\n'
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

temp_files += [f for f in os.listdir('.') if f.endswith('.int') or f.endswith('.rep')]
for f in temp_files:
    if os.path.isfile(f):
        os.remove(f)

if fail_count:
    print('{} tests failed'.format(fail_count))
else:
    print('all tests passed')
    
sys.exit(fail_count)
