#!/usr/bin/env bash
# Benchmark script for treetime timetree command
# Cycles through small datasets and coalescent settings, kills any run > 30s

set -euo pipefail

BIN="./target/release/treetime"
OUTDIR="auspice"
TIMEOUT_SEC=30

# Datasets: "label tree aln dates"
DATASETS=(
    "ebola_20    data/ebola/20/tree.nwk       data/ebola/20/aln.fasta.xz       data/ebola/20/metadata.tsv"
    "ebola_100   data/ebola/100/tree.nwk      data/ebola/100/aln.fasta.xz      data/ebola/100/metadata.tsv"
    "ebola_362   data/ebola/362/tree.nwk      data/ebola/362/aln.fasta.xz      data/ebola/362/metadata.tsv"
    "flu_h3n2_20  data/flu/h3n2/20/tree.nwk   data/flu/h3n2/20/aln.fasta.xz   data/flu/h3n2/20/metadata.tsv"
    "flu_h3n2_200 data/flu/h3n2/200/tree.nwk  data/flu/h3n2/200/aln.fasta.xz  data/flu/h3n2/200/metadata.tsv"
)

# Coalescent settings: "label flag_args"
# "off" = no coalescent flag; numeric/skyline passed via --coalescent
COALESCENT=(
    "off      "
    # "skyline  --coalescent skyline"
    "0.001    --coalescent 0.001"
    "0.01     --coalescent 0.01"
    "0.1      --coalescent 0.1"
    # "1.0      --coalescent 1.0"
)

mkdir -p "$OUTDIR"

# Result tracking: "label status elapsed"
declare -a RESULTS

pad() { printf "%-${2}s" "$1"; }

echo "=========================================================="
echo " treetime timetree benchmark  (timeout: ${TIMEOUT_SEC}s)"
echo "=========================================================="
printf "%-20s %-10s %-8s %s\n" "DATASET" "COALESCENT" "STATUS" "ELAPSED"
echo "----------------------------------------------------------"

for ds_entry in "${DATASETS[@]}"; do
    read -r label tree aln dates <<< "$ds_entry"

    for coal_entry in "${COALESCENT[@]}"; do
        coal_label=$(awk '{print $1}' <<< "$coal_entry")
        coal_flag=$(sed 's/^[^ ]*//' <<< "$coal_entry" | xargs)

        out_json="${OUTDIR}/${label}_coal${coal_label}.json"

        cmd=(timeout "$TIMEOUT_SEC" "$BIN" timetree -j1 \
            --tree "$tree" \
            --aln "$aln" \
            --dates "$dates" \
            --output-tree-auspice "$out_json" \
            --ladderize ascending)
        [[ -n "$coal_flag" ]] && cmd+=($coal_flag)

        t_start=$(date +%s%N)
        set +e
        "${cmd[@]}" > "${out_json%.json}.log" 2>&1
        exit_code=$?
        set -e
        t_end=$(date +%s%N)

        elapsed_ms=$(( (t_end - t_start) / 1000000 ))
        elapsed_s=$(printf "%.1fs" "$(echo "scale=1; $elapsed_ms/1000" | bc)")

        if [[ $exit_code -eq 124 ]]; then
            status="TIMEOUT"
        elif [[ $exit_code -eq 0 ]]; then
            status="ok"
        else
            status="FAIL($exit_code)"
        fi

        printf "%-20s %-10s %-8s %s\n" "$label" "$coal_label" "$status" "$elapsed_s"
        RESULTS+=("${label} ${coal_label} ${status} ${elapsed_s}")
    done
    echo ""
done

echo "=========================================================="
echo "Output JSON files in: $OUTDIR/"
echo "Logs:                 ${OUTDIR}/<name>.log"
echo "=========================================================="
