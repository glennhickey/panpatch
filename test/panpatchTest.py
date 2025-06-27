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
temp_files += ['tiny.x2.vg', 'tiny.vg', 'tiny.truth.fa']

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

for f in temp_files:
    if os.path.isfile(f):
        os.remove(f)

if fail_count:
    print('{} tests failed'.format(fail_count))
else:
    print('all tests passed')
    
sys.exit(fail_count)
