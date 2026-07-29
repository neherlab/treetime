#!/usr/bin/env bash
# Dispatch a single treetime invocation for benchmarking.
# Usage: run_one.sh <binary> <command> <dataset> <jobs> <outdir>
#   <command>: anc-marginal | anc-parsimony | clock | timetree
#   <dataset>: path under DATA_DIR, e.g. flu/h3n2/500
set -euo pipefail

THIS_DIR="$(dirname "$(realpath "${BASH_SOURCE[0]}")")"
REPO_ROOT="$(realpath "${THIS_DIR}/../../../..")"

BIN="$1"; CMD="$2"; DS="$3"; JOBS="$4"; OUT="$5"
DATA_DIR="${DATA_DIR:-${REPO_ROOT}/data}"

TREE="$DATA_DIR/$DS/tree.nwk"
ALN="$DATA_DIR/$DS/aln.fasta.xz"
META="$DATA_DIR/$DS/metadata.tsv"

# Per-dataset metadata id column (tree tip labels). Dengue tips are genbank
# accessions, not the default `strain` column; flu tips match the default.
IDCOL=()
case "$DS" in
  dengue/*) IDCOL=(--metadata-id-columns genbank_accession) ;;
esac

mkdir -p "$OUT"

case "$CMD" in
  anc-marginal)
    exec "$BIN" --jobs="$JOBS" ancestral --method-anc=marginal \
      -t "$TREE" -a "$ALN" -O "$OUT" ;;
  anc-marginal-dense)
    exec "$BIN" --jobs="$JOBS" ancestral --method-anc=marginal --dense=true \
      -t "$TREE" -a "$ALN" -O "$OUT" ;;
  anc-parsimony)
    exec "$BIN" --jobs="$JOBS" ancestral --method-anc=parsimony \
      -t "$TREE" -a "$ALN" -O "$OUT" ;;
  clock)
    exec "$BIN" --jobs="$JOBS" clock \
      -t "$TREE" --dates "$META" "${IDCOL[@]}" -O "$OUT" ;;
  timetree)
    exec "$BIN" --jobs="$JOBS" timetree \
      -t "$TREE" -a "$ALN" --dates "$META" "${IDCOL[@]}" -O "$OUT" ;;
  *)
    echo "unknown command: $CMD" >&2; exit 2 ;;
esac
