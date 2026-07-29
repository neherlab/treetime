#!/usr/bin/env bash
# Run the hyperfine wall-clock matrix: each (command,dataset) over
# rev={before,after} x jobs={1,2,4,8,16}. Exports one JSON per (command,dataset).
set -euo pipefail

THIS_DIR="$(dirname "$(realpath "${BASH_SOURCE[0]}")")"
REPO_ROOT="$(realpath "${THIS_DIR}/../../../..")"
OUTDIR="${1:-${REPO_ROOT}/tmp/perf-pangraph/bench}"
JOBS_LIST="${JOBS_LIST:-1,2,4,8,16}"

mkdir -p "$OUTDIR/hyperfine"
SCRATCH="$OUTDIR/scratch"
mkdir -p "$SCRATCH"

# matrix rows: "<cmd> <dataset> <runs> <warmup>"
MATRIX=(
  "anc-marginal  dengue/2000   12 3"
  "anc-marginal  dengue/1000   12 3"
  "anc-parsimony dengue/2000   12 3"
  "anc-parsimony dengue/1000   12 3"
  "timetree      dengue/2000   10 2"
  "timetree      flu/h3n2/500  12 3"
  "clock         dengue/2000   15 3"
)

for row in "${MATRIX[@]}"; do
  read -r cmd ds runs warmup <<<"$row"
  tag="${cmd}__${ds//\//-}"
  echo ">>> hyperfine $cmd $ds (runs=$runs warmup=$warmup)"
  hyperfine \
    --shell=none \
    --warmup "$warmup" --runs "$runs" \
    --parameter-list jobs "$JOBS_LIST" \
    --parameter-list rev before,after \
    --command-name "{rev} j={jobs}" \
    --export-json "$OUTDIR/hyperfine/${tag}.json" \
    "$THIS_DIR/run_bench.sh {rev} $cmd $ds {jobs} $SCRATCH/{rev}-{jobs}"
done

echo "hyperfine matrix complete -> $OUTDIR/hyperfine"
