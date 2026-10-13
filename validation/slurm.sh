#!/usr/bin/env bash
# Submit the whole validation pipeline as Slurm jobs, chained by dependencies:
#   generate (array: tips x coal x rep) -> run (array per density and suite: one task per sequence set)
#   -> evaluate (one job per density and suite) -> analyze (one job per density).
# Every step skips finished work, so resubmitting continues an interrupted pipeline.
#
# Usage: ./slurm.sh            (build ../target/release/treetime first)
# Settings (environment, defaults in brackets):
#   REPS [10]  TIPS ["10 30 100"]  SUITES ["ancestral timetree optimize"]  TIMETREE_SET [core]
#   RUN_CPUS [4] cores per run task, RUN_MEM_PER_CPU [2G], RUN_TIME [2:00:00]
#   GEN_TIME [1:00:00]  EVAL_CPUS [32]  EVAL_TIME [4:00:00]
#   OPT_FILTER ['^(?!.*fixed-topo-dense)']: dense optimize needs > 8 GB per run at 4000 tips
#   CONDA_ENV [nextstrain]: environment with biopython, numpy, scipy, pandas, matplotlib, tabulate
#   PY [python of CONDA_ENV]: resolved on the submitting (head) node; compute nodes run it directly, without network
# Partition, account, QOS: Slurm's own SBATCH_PARTITION, SBATCH_ACCOUNT, SBATCH_QOS.
# Simulations use the legacy TreeTime of this repository (PYTHONPATH), as for the existing data.
set -euo pipefail
cd "$(dirname "$0")"

REPS=${REPS:-10}
TIPS=${TIPS:-"10 30 100"}
SUITES=${SUITES:-"ancestral timetree optimize"}
TIMETREE_SET=${TIMETREE_SET:-core}
RUN_CPUS=${RUN_CPUS:-4}
RUN_MEM_PER_CPU=${RUN_MEM_PER_CPU:-2G}
RUN_TIME=${RUN_TIME:-2:00:00}
GEN_TIME=${GEN_TIME:-1:00:00}
EVAL_CPUS=${EVAL_CPUS:-32}
EVAL_TIME=${EVAL_TIME:-4:00:00}
OPT_FILTER=${OPT_FILTER:-'^(?!.*fixed-topo-dense)'}
CONDA_ENV=${CONDA_ENV:-nextstrain}
PY=${PY:-$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)')}
N_SETS=$((3 * REPS * 8))  # coal x rep x (L, mu); run.py shards over the sets that exist

test -x ../target/release/treetime || { echo "build first: cargo build --release --bin treetime" >&2; exit 1; }
mkdir -p slurm_logs
export PYTHONPATH="$PWD/../packages/legacy/treetime" OPENBLAS_NUM_THREADS=1
"$PY" -c 'import Bio, numpy, scipy, pandas, matplotlib, tabulate; from treetime.seqgen import SeqGen' ||
  { echo "missing Python dependencies in $PY" >&2; exit 1; }
echo "python: $PY"

submit() {  # submit <name> <dependency or ''> <sbatch options...> -- <command>
  local name=$1 dep=$2; shift 2
  local opts=()
  while [[ $1 != -- ]]; do opts+=("$1"); shift; done
  shift
  sbatch --parsable --job-name="tt-$name" --output="slurm_logs/%x_%A_%a.out" ${dep:+--dependency=afterany:$dep} \
    "${opts[@]}" --wrap="export PYTHONPATH=$PYTHONPATH OPENBLAS_NUM_THREADS=1; $*"  # also with SBATCH_EXPORT=NONE
}

# generate: one task per (tips, coal, rep), largest trees first
jobs_file=slurm_logs/generate_jobs.txt
for n in $(printf '%s\n' $TIPS | sort -rn); do
  out=data; [[ $n == 10 ]] || out=data_n$n
  for coal in const exp piecewise; do for rep in $(seq 0 $((REPS - 1))); do echo "$coal $rep $n $out"; done; done
done >"$jobs_file"
gen=$(submit generate '' --array=1-"$(wc -l <"$jobs_file")" --cpus-per-task=1 --mem=4G --time="$GEN_TIME" -- \
  'set -- $(sed -n "${SLURM_ARRAY_TASK_ID}p" '"$jobs_file"'); '"$PY"' simulate.py --coal $1 --rep $2 --tips-per-sample $3 --out $4')
echo "generate: job $gen"

for n in $TIPS; do
  data=data; results=results
  [[ $n == 10 ]] || { data=data_n$n; results=results_n$n; }
  evals=()
  for suite in $SUITES; do
    extra=()
    [[ $suite == timetree ]] && extra=(--set "$TIMETREE_SET") || extra=(--set core)
    [[ $suite == optimize && -n $OPT_FILTER ]] && extra+=(--filter "'$OPT_FILTER'")
    run=$(submit "run-$suite-n$n" "$gen" --array=0-$((N_SETS - 1)) --cpus-per-task="$RUN_CPUS" \
      --mem-per-cpu="$RUN_MEM_PER_CPU" --time="$RUN_TIME" -- \
      "$PY run.py --suite $suite ${extra[*]} --data $data --results $results -j $RUN_CPUS" \
      '--shard ${SLURM_ARRAY_TASK_ID}/'"$N_SETS")
    ev=$(submit "evaluate-$suite-n$n" "$run" --cpus-per-task="$EVAL_CPUS" --mem-per-cpu=1G --time="$EVAL_TIME" -- \
      "$PY evaluate.py --suite $suite --data $data --results $results -j $EVAL_CPUS")
    evals+=("$ev")
    echo "n=$n $suite: run $run, evaluate $ev"
  done
  an=$(submit "analyze-n$n" "$(IFS=:; echo "${evals[*]}")" --cpus-per-task=1 --mem=8G --time=1:00:00 -- \
    "$PY analyze.py --results $results")
  echo "n=$n analyze: job $an -> $results/report/"
done
