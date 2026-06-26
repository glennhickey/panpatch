#!/usr/bin/env bash
# Run panpatch in self-reference (no-CHM13) mode on every chromosome graph in a directory.
# reference == target (the verkko haplotype itself); donor is the diploid hifiasm of the same sample.
# usage: ./run-nochm13.sh <chroms_dir> <ref/target_sample> <donor_sample> [out_dir]
#   e.g. ./run-nochm13.sh HG01074.no-chm13.1.chroms verkko-R3-HG01074_1 hifiasm-R2-HG01074
set -uo pipefail
PANPATCH="$HOME/dev/panpatch/panpatch"
DIR="$1"; REF="$2"; DONOR="$3"
OUT="${4:-${DIR%.chroms}.runs}"
THREADS=8
JOBS=3
mkdir -p "$OUT"

run_one() {
  local vg base
  vg="$1"; base=$(basename "$vg" .full.vg)
  [ "$base" = "chrOther" ] && return 0
  [ -f "$OUT/$base.status" ] && { echo "[skip] $base"; return 0; }
  echo "[run ] $base"
  "$PANPATCH" "$vg" -r "$REF" -s "$REF" -s "$DONOR" -T -t "$THREADS" -p > "$OUT/$base.bed" 2> "$OUT/$base.stderr"
  echo "$?" > "$OUT/$base.status"
  echo "[done] $base"
}
export -f run_one; export PANPATCH OUT REF DONOR THREADS

ls "$DIR"/haplotype*.full.vg | xargs -P "$JOBS" -I{} bash -c 'run_one "$@"' _ {}

# classify (one bed == one chromosome/haplotype in self-reference mode)
printf "contig\texit\tresult\n" > "$OUT/summary.tsv"
for b in "$OUT"/haplotype*.bed; do
  base=$(basename "$b" .bed); rc=$(cat "$OUT/$base.status" 2>/dev/null || echo "?")
  if   [ "$rc" != 0 ];                                          then res="ERROR"
  elif grep -q '^#Reverting patch'             "$b";           then res="GUARD_REVERTED"
  elif grep -q '^#Reverting failed patch'      "$b";           then res="LENGTH_REVERTED"
  elif grep -q '^#Telomere validation failed'  "$b";           then res="TELO_VALIDATION_FAILED"
  elif grep -q '^#Reverting to input'          "$b";           then res="REVERTED_HAS_GAPS"
  elif grep -q '^#Telomere patch'              "$b";           then res="TELOMERE_PATCHED"
  elif grep -q '^#Interior graft'              "$b";           then res="GAP_FILLED"
  elif grep -q '^#No patching is required'     "$b";           then res="NO_PATCH_NEEDED"
  else                                                              res="OTHER"
  fi
  printf "%s\t%s\t%s\n" "$base" "$rc" "$res" >> "$OUT/summary.tsv"
done
echo "=== SUMMARY ($OUT/summary.tsv) ==="; column -t "$OUT/summary.tsv"
