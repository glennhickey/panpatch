#!/usr/bin/env bash
# Run panpatch -T on every chromosome of the given sample graph(s).
# Target = verkko-R3-<id> (patched), donor = hifiasm-R2-<id>, reference = CHM13.
# Keeps the BED (stdout) and stderr per run; does NOT write FASTA.
# Resumable: a run that finished (ok or error) is recorded and skipped on re-run.
#
# usage: ./run-patch.sh [HG00126 HG00235 HG01074]   (no args = all three)
set -uo pipefail

PANPATCH="$HOME/dev/panpatch/panpatch"
ROOT="$HOME/dev/work/panpatch-jun24"
OUT="${PATCHOUT:-$ROOT/patch-runs}"
THREADS=8
mkdir -p "$OUT"

if [ "$#" -gt 0 ]; then SAMPLES=("$@"); else SAMPLES=(HG00126 HG00235 HG01074); fi

summary="$OUT/summary.${1:-all}.tsv"   # per-invocation so parallel runs don't race on it
[ -f "$summary" ] || printf "sample\tchrom\texit\tresult\n" > "$summary"

# natural chromosome order (1..22, X, Y, M)
chrom_key() { local c=${1#chr}; case "$c" in X) echo 23;; Y) echo 24;; M) echo 25;; *) echo "$c";; esac; }

for id in "${SAMPLES[@]}"; do
  dir="$ROOT/$id.chroms"
  verkko="verkko-R3-$id"
  hifiasm="hifiasm-R2-$id"
  [ -d "$dir" ] || { echo "[warn] missing $dir"; continue; }

  # iterate chromosomes in natural order
  mapfile -t vgs < <(for f in "$dir"/chr*.full.vg; do printf "%s\t%s\n" "$(chrom_key "$(basename "$f" .full.vg)")" "$f"; done | sort -n -k1,1 | cut -f2)

  for vg in "${vgs[@]}"; do
    chr=$(basename "$vg" .full.vg)
    bed="$OUT/$id.$chr.bed"
    err="$OUT/$id.$chr.stderr"
    mark="$OUT/$id.$chr.status"      # records exit code; presence => done
    if [ -f "$mark" ]; then
      echo "[skip] $id $chr (done: $(cat "$mark"))"
      continue
    fi
    echo "[run ] $id $chr ..."
    "$PANPATCH" "$vg" -r CHM13 -s "$verkko" -s "$hifiasm" -T -t "$THREADS" -p > "$bed" 2> "$err"
    rc=$?

    if   [ $rc -ne 0 ];                                      then res="ERROR(rc=$rc)"
    elif grep -q '^#Telomere patch'           "$bed";       then res="TELOMERE_PATCHED"
    elif grep -q '^#Telomere not patched'     "$bed";       then res="TELOMERE_SKIPPED_CAP"
    elif grep -q '^#Telomere validation failed' "$bed";     then res="TELO_VALIDATION_FAILED"
    elif grep -q '^#Reverting'                "$bed";       then res="REVERTED"
    elif grep -q '^#No patching is required'  "$bed";       then res="NO_PATCH_NEEDED"
    else                                                          res="PATCHED_OTHER"
    fi

    printf "%s\t%s\t%s\t%s\n" "$id" "$chr" "$rc" "$res" >> "$summary"
    echo "$rc" > "$mark"
    echo "[done] $id $chr -> rc=$rc  $res"
  done
done

echo "=== SUMMARY ($summary) ==="
column -t "$summary"
