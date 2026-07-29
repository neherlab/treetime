#!/usr/bin/env bash
# One-shot: perf stat (HW counters) at thread endpoints + perf record profiles.
# Run while the system is idle: bash perf_capture.sh
set -euo pipefail

THIS_DIR="$(dirname "$(realpath "${BASH_SOURCE[0]}")")"
REPO_ROOT="$(realpath "${THIS_DIR}/../../../..")"
BEFORE_WT="${BEFORE_WT:?Set BEFORE_WT to the baseline worktree path}"
AFTER_WT="${AFTER_WT:?Set AFTER_WT to the feature worktree path}"
BREL="${BEFORE_WT}/.out/treetime"
AREL="${AFTER_WT}/.out/treetime"
BPROF="${BEFORE_WT}/.build/docker/profiling/treetime"
APROF="${AFTER_WT}/.build/docker/profiling/treetime"
DATA="${DATA:-${REPO_ROOT}/data}"
OUT="${OUT:-${REPO_ROOT}/tmp/perf-pangraph/bench}"
mkdir -p "$OUT/perfstat" "$OUT/timev" "$OUT/profile" "$OUT/scratch/perf"
EV=task-clock,context-switches,cpu-migrations,page-faults,cycles,instructions,branches,branch-misses,cache-references,cache-misses

# Build treetime argument list for a (cmd, dataset, jobs) into global array ARGS.
build_args() {
  local cmd="$1" ds="$2" jobs="$3"
  local tree="$DATA/$ds/tree.nwk" aln="$DATA/$ds/aln.fasta.xz" meta="$DATA/$ds/metadata.tsv"
  local idcol=()
  [[ "$ds" == dengue/* ]] && idcol=(--metadata-id-columns genbank_accession)
  case "$cmd" in
    timetree)      ARGS=(--jobs="$jobs" timetree -t "$tree" -a "$aln" --dates "$meta" "${idcol[@]}" -O "$OUT/scratch/perf") ;;
    anc-marginal)  ARGS=(--jobs="$jobs" ancestral --method-anc=marginal  -t "$tree" -a "$aln" -O "$OUT/scratch/perf") ;;
    anc-parsimony) ARGS=(--jobs="$jobs" ancestral --method-anc=parsimony -t "$tree" -a "$aln" -O "$OUT/scratch/perf") ;;
  esac
}

# perf stat + /usr/bin/time -v at endpoints {1,16} for the cells with real signal.
for cell in "timetree dengue/2000" "anc-marginal dengue/2000" "anc-parsimony dengue/2000"; do
  read -r cmd ds <<<"$cell"
  dstag="${ds//\//-}"
  for j in 1 16; do
    for rev in before after; do
      bin="$BREL"; [[ "$rev" == after ]] && bin="$AREL"
      build_args "$cmd" "$ds" "$j"
      base="${cmd}__${dstag}__${rev}__j${j}"
      "$bin" "${ARGS[@]}" >/dev/null 2>&1 || true   # warm caches
      perf stat -x, -e "$EV" -o "$OUT/perfstat/${base}.csv" -- "$bin" "${ARGS[@]}" >/dev/null 2>&1 || true
      /usr/bin/time -v -o "$OUT/timev/${base}.txt" "$bin" "${ARGS[@]}" >/dev/null 2>&1 || true
      echo "perfstat $base"
    done
  done
done

# perf record call-graph profiles at representative points (profiling binaries).
profile() {
  local label="$1" bin="$2" cmd="$3" ds="$4" jobs="$5"
  build_args "$cmd" "$ds" "$jobs"
  local data="$OUT/profile/${label}.perf.data"
  perf record --call-graph dwarf -F 1999 -o "$data" -- "$bin" "${ARGS[@]}" >/dev/null 2>"$OUT/profile/${label}.record.log" || true
  perf report -i "$data" --stdio --sort overhead,symbol 2>/dev/null | grep -vE '^\s*#' | head -40 > "$OUT/profile/${label}.top-symbols.txt" || true
  echo "profiled $label"
}
profile tt-dengue2000-j8-before "$BPROF" timetree      dengue/2000 8
profile tt-dengue2000-j8-after  "$APROF" timetree      dengue/2000 8
profile marg-dengue2000-j8-after "$APROF" anc-marginal  dengue/2000 8
profile fitch-dengue2000-j16-after "$APROF" anc-parsimony dengue/2000 16
echo "perf_capture complete"
