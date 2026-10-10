#!/usr/bin/env bash
# Generate the full validation grid: 3 coalescent schedules x 10 replicate trees,
# each with sequences for L in {1000, 10000} x mu in {0.01, 0.003, 0.001, 0.0003}.
# Usage: ./generate_all.sh [n_reps] [n_jobs]
set -euo pipefail
cd "$(dirname "$0")"
export PYTHONPATH="$PWD/../packages/legacy/treetime" OPENBLAS_NUM_THREADS=1
reps=${1:-10}; jobs=${2:-8}
for coal in const exp piecewise; do for rep in $(seq 0 $((reps - 1))); do echo "$coal $rep"; done; done |
  xargs -P "$jobs" -L 1 sh -c 'uv run -q --no-project --with biopython,numpy,scipy,pandas python simulate.py --coal $0 --rep $1 --out data'
