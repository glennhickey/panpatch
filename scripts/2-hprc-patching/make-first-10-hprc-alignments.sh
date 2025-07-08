#!/usr/bin/bash

for SAMPLE in HG00097 HG00099 HG00126 HG00128 HG00133 HG00140 HG00146 HG00232 HG00235 HG00253; do
    ./make-hprc-alignments.sh ${SAMPLE} &
done
wait
