#!/usr/bin/env bash
# Generate the full validation grid: 3 coalescent schedules x 10 replicate trees,
# each with sequences for L in {1000, 10000} x mu in {0.01, 0.003, 0.001, 0.0003}.
# Usage: pixi run ./generate_all.sh [n_reps] [n_jobs (default: available CPUs)]
set -euo pipefail
cd "$(dirname "$0")"
cpus=$(awk '$1 != "max" {print int($1 / $2); exit} {exit 1}' /sys/fs/cgroup/cpu.max 2>/dev/null || nproc) # cgroup quota, which nproc does not see
reps=${1:-10}
jobs=${2:-$cpus}
for coal in const exp piecewise; do for rep in $(seq 0 $((reps - 1))); do echo "--coal $coal --rep $rep"; done; done |
  xargs -P "$jobs" -L 1 python simulate.py --out data
