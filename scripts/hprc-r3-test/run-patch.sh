#!/usr/bin/env bash
# CHM13-referenced patching: ONE multi-graph panpatch run per sample over all its chromosome graphs.
# Target = verkko-R3-<id> (patched), donor = hifiasm-R2-<id>, reference = CHM13.
# Writes the deliverable layout directly:
#   deliverable/chm13/<id>.report   (the TSV report, stdout)
#   deliverable/chm13/<id>.bed      (patched-assembly intervals)
#   deliverable/chm13/<id>.hap1.fa, <id>.hap2.fa   (per-haplotype FASTA, from -f)
#   deliverable/chm13/<id>.stderr
# Resumable: a sample whose report already exists is skipped.
#
# usage: ./run-patch.sh [HG00126 HG00235 HG01074]      (no args = all three)
#   env: ROOT (dir holding <id>.chroms/, default .), PATCHOUT (output dir, default deliverable),
#        THREADS (default 8), PANPATCH (default ~/dev/panpatch/panpatch)
set -uo pipefail

PANPATCH="${PANPATCH:-$HOME/dev/panpatch/panpatch}"
ROOT="${ROOT:-.}"
OUT="${PATCHOUT:-deliverable}"
THREADS="${THREADS:-8}"
mkdir -p "$OUT/chm13"

if [ "$#" -gt 0 ]; then SAMPLES=("$@"); else SAMPLES=(HG00126 HG00235 HG01074); fi

for id in "${SAMPLES[@]}"; do
  dir="$ROOT/$id.chroms"
  [ -d "$dir" ] || { echo "[warn] missing $dir"; continue; }
  report="$OUT/chm13/$id.report"
  if [ -f "$report" ]; then echo "[skip] $id (have $report)"; continue; fi
  echo "[run ] $id ..."
  "$PANPATCH" "$dir"/chr*.full.vg -r CHM13 -s "verkko-R3-$id" -s "hifiasm-R2-$id" \
      -T -t "$THREADS" -p --bed "$OUT/chm13/$id.bed" -f "$OUT/chm13/$id" \
      > "$report" 2> "$OUT/chm13/$id.stderr"
  rc=$?
  echo "[done] $id rc=$rc"
  if [ "$rc" -ne 0 ]; then echo "  ERROR (see $OUT/chm13/$id.stderr)"; continue; fi
  # per-sample summary, straight from the report TSV (type=$3, decision=$12)
  awk -F'\t' 'NR>1 && $3!="" && $3!="type" {c[$3" "$12]++}
       END{for (k in c) printf "    %-22s %d\n", k, c[k]}' "$report" | sort
done
