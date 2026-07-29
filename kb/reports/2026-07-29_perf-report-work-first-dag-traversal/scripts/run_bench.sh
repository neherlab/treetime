#!/usr/bin/env bash
# Map a revision label to its binary and dispatch a benchmark invocation.
# Usage: run_bench.sh <rev> <command> <dataset> <jobs> <outdir> [flavor]
#   <rev>: before | after
#   <flavor>: release (default) | profiling
set -euo pipefail

THIS_DIR="$(dirname "$(realpath "${BASH_SOURCE[0]}")")"

REV="$1"; CMD="$2"; DS="$3"; JOBS="$4"; OUT="$5"; FLAVOR="${6:-release}"

BEFORE_WT="${BEFORE_WT:?Set BEFORE_WT to the baseline worktree path}"
AFTER_WT="${AFTER_WT:?Set AFTER_WT to the feature worktree path}"

case "$REV" in
  before) WT="$BEFORE_WT" ;;
  after)  WT="$AFTER_WT" ;;
  *) echo "unknown rev: $REV" >&2; exit 2 ;;
esac

case "$FLAVOR" in
  release)   BIN="$WT/.out/treetime" ;;
  profiling) BIN="$WT/.build/docker/profiling/treetime" ;;
  *) echo "unknown flavor: $FLAVOR" >&2; exit 2 ;;
esac

exec "$THIS_DIR/run_one.sh" "$BIN" "$CMD" "$DS" "$JOBS" "$OUT"
