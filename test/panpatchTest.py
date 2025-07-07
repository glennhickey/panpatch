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



for f in temp_files:
    if os.path.isfile(f):
        os.remove(f)

if fail_count:
    print('{} tests failed'.format(fail_count))
else:
    print('all tests passed')
    
sys.exit(fail_count)
