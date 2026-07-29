#!/usr/bin/env bash
# Dense-marginal wall-clock matrix. Dense per-node work uses full probability
# vectors and OpenBLAS.
set -euo pipefail

THIS_DIR="$(dirname "$(realpath "${BASH_SOURCE[0]}")")"
REPO_ROOT="$(realpath "${THIS_DIR}/../../../..")"
OUTDIR="${1:-${REPO_ROOT}/tmp/perf-pangraph/bench}"
JOBS_LIST="${JOBS_LIST:-1,2,4,8,16}"

mkdir -p "$OUTDIR/hyperfine"
SCRATCH="$OUTDIR/scratch"
mkdir -p "$SCRATCH"

# rows: "<cmd> <dataset> <runs> <warmup>"
# dengue/2000 dense dropped: genome-scale dense is atypical (that is sparse's job),
# single-run probe already put it at +/-5% (break-even, BLAS-bound), and j1 alone is
# ~16-34s/run. flu/500 is the representative dense scale; dengue/1000 is the stress point.
MATRIX=(
  "anc-marginal-dense flu/h3n2/500 10 2"
  "anc-marginal-dense dengue/1000   6 1"
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

echo "dense hyperfine matrix complete -> $OUTDIR/hyperfine"
