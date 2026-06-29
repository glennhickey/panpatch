#!/usr/bin/env bash
# Self-reference (no-CHM13) patching: ONE multi-graph panpatch run per (sample, haplotype).
# reference == target == the verkko haplotype itself; donor = diploid hifiasm of the same sample.
# Writes the deliverable layout directly:
#   deliverable/self-ref/<id>.hap<h>/report      (the TSV report, stdout)
#   deliverable/self-ref/<id>.hap<h>/out.bed     (patched-assembly intervals)
#   deliverable/self-ref/<id>.hap<h>/out.hap0.fa (FASTA, from -f; self-ref haplotype field is 0)
#   deliverable/self-ref/<id>.hap<h>/stderr
# Resumable: skipped if the report already exists.
#
# usage: ./run-nochm13.sh <chroms_dir> <ref/target_sample> <donor_sample> [out_dir]
#   e.g. ./run-nochm13.sh HG01074.no-chm13.1.chroms verkko-R3-HG01074_1 hifiasm-R2-HG01074
#   env: PANPATCH (default ~/dev/panpatch/panpatch), THREADS (default 8)
set -uo pipefail

PANPATCH="${PANPATCH:-$HOME/dev/panpatch/panpatch}"
DIR="$1"; REF="$2"; DONOR="$3"
OUT="${4:-deliverable}"
THREADS="${THREADS:-8}"

# derive <id> and <h> from REF = verkko-R3-<id>_<h>
idh="${REF#verkko-R3-}"            # <id>_<h>
id="${idh%_*}"; h="${idh##*_}"
dest="$OUT/self-ref/$id.hap$h"
mkdir -p "$dest"

report="$dest/report"
if [ -f "$report" ]; then echo "[skip] $id hap$h (have $report)"; exit 0; fi
echo "[run ] $id hap$h ..."
# *.full.vg includes the chrOther graph (many unplaced contigs) -> panpatch passes it through; keep it.
"$PANPATCH" "$DIR"/*.full.vg -r "$REF" -s "$REF" -s "$DONOR" \
    -T -t "$THREADS" -p --bed "$dest/out.bed" -f "$dest/out" \
    > "$report" 2> "$dest/stderr"
rc=$?
echo "[done] $id hap$h rc=$rc"
[ "$rc" -ne 0 ] && { echo "  ERROR (see $dest/stderr)"; exit "$rc"; }
# per-run summary, straight from the report TSV (type=$3, decision=$12)
awk -F'\t' 'NR>1 && $3!="" && $3!="type" {c[$3" "$12]++}
     END{for (k in c) printf "    %-22s %d\n", k, c[k]}' "$report" | sort
