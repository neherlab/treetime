#!/usr/bin/env bash
# Benchmark treetime optimize and ancestral across core counts.
# Usage: ./bench_scaling.sh [dataset_dir]
#
# dataset_dir must contain tree.nwk and aln.fasta.xz
# Default: data/mpox/clade-ii/2000
#
# Example with a different dataset:
#   ./bench_scaling.sh data/rsv/a/2000
#   ./bench_scaling.sh data/flu/h3n2/200

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BINARY="${BINARY:-target/release/treetime}"
DATASET="${1:-data/mpox/clade-ii/2000}"
CORE_COUNTS="${CORE_COUNTS:-1 2 4 8}"   # space-separated list
RUNS="${RUNS:-3}"                        # repetitions per configuration (median is reported)
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------
TREE="$DATASET/tree.nwk"
ALN="$DATASET/aln.fasta.xz"

if [[ ! -f "$BINARY" ]]; then
  echo "ERROR: binary not found: $BINARY" >&2
  echo "Run: cargo build --release" >&2
  exit 1
fi
if [[ ! -f "$TREE" ]]; then
  echo "ERROR: tree not found: $TREE" >&2
  exit 1
fi
if [[ ! -f "$ALN" ]]; then
  echo "ERROR: alignment not found: $ALN" >&2
  exit 1
fi

if ! command -v /usr/bin/time &>/dev/null; then
  echo "ERROR: /usr/bin/time not found (need GNU time for -v flag)" >&2
  exit 1
fi

echo "Binary:  $BINARY"
echo "Dataset: $DATASET"
echo "Cores:   $CORE_COUNTS"
echo "Runs:    $RUNS"
echo

# ---------------------------------------------------------------------------
# Helper: run one timed measurement, echo "wall_s rss_kb cpu_pct"
# ---------------------------------------------------------------------------
measure_one() {
  local cmd=("$@")
  local timefile="$TMP_DIR/time_$$_${RANDOM}.txt"

  # Discard stdout/stderr from the binary; capture /usr/bin/time -v output
  /usr/bin/time -v "${cmd[@]}" \
    >"$TMP_DIR/out.txt" 2>"$timefile" \
    || true   # don't abort on non-zero exit from the benchmarked command

  local wall_raw rss_kb cpu_pct wall_s
  wall_raw=$(grep "Elapsed (wall clock)" "$timefile" | awk '{print $NF}')
  rss_kb=$(grep "Maximum resident set size" "$timefile" | awk '{print $NF}')
  cpu_pct=$(grep "Percent of CPU this job got" "$timefile" | tr -d '%' | awk '{print $NF}')

  # Convert m:ss.ss or h:mm:ss to seconds
  if [[ "$wall_raw" =~ ^([0-9]+):([0-9]+):([0-9.]+)$ ]]; then
    wall_s=$(awk "BEGIN{printf \"%.2f\", ${BASH_REMATCH[1]}*3600 + ${BASH_REMATCH[2]}*60 + ${BASH_REMATCH[3]}}")
  elif [[ "$wall_raw" =~ ^([0-9]+):([0-9.]+)$ ]]; then
    wall_s=$(awk "BEGIN{printf \"%.2f\", ${BASH_REMATCH[1]}*60 + ${BASH_REMATCH[2]}}")
  else
    wall_s="$wall_raw"
  fi

  echo "$wall_s $rss_kb $cpu_pct"
}

# ---------------------------------------------------------------------------
# Helper: run RUNS repetitions and report median wall, max rss, median cpu
# ---------------------------------------------------------------------------
benchmark() {
  local label="$1"; shift
  local cores="$1"; shift
  local cmd=("$@")

  local wall_vals=() rss_vals=() cpu_vals=()
  for (( i=1; i<=RUNS; i++ )); do
    read -r w r c < <(measure_one "${cmd[@]}")
    wall_vals+=("$w")
    rss_vals+=("$r")
    cpu_vals+=("$c")
  done

  # Sort and pick median index
  local sorted_wall sorted_rss sorted_cpu
  sorted_wall=($(printf '%s\n' "${wall_vals[@]}" | sort -n))
  sorted_rss=($(printf '%s\n' "${rss_vals[@]}" | sort -n))
  sorted_cpu=($(printf '%s\n' "${cpu_vals[@]}" | sort -n))

  local mid=$(( (RUNS - 1) / 2 ))
  local med_wall="${sorted_wall[$mid]}"
  local med_rss="${sorted_rss[$mid]}"
  local med_cpu="${sorted_cpu[$mid]}"
  local rss_mb
  rss_mb=$(awk "BEGIN{printf \"%.0f\", $med_rss/1024}")

  printf "%-12s  j=%-3s  wall=%8.2fs  rss=%6s MiB  cpu=%s%%\n" \
    "$label" "$cores" "$med_wall" "$rss_mb" "$med_cpu"
}

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
echo "=== treetime optimize ==="
for j in $CORE_COUNTS; do
  OUT_NWK="$TMP_DIR/annotated_tree_j${j}.nwk"
  benchmark "optimize" "$j" \
    "$BINARY" optimize \
      --tree "$TREE" \
      --aln "$ALN" \
      -j "$j" \
      --output-tree-nwk-annotated "$OUT_NWK"
done

echo
echo "=== treetime ancestral ==="
for j in $CORE_COUNTS; do
  OUT_NWK="$TMP_DIR/annotated_tree_j${j}.nwk"
  benchmark "ancestral" "$j" \
    "$BINARY" ancestral \
      --tree "$TREE" \
      --aln "$ALN" \
      -j "$j" \
      --output-tree-nwk-annotated "$OUT_NWK"
done
