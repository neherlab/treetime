#!/usr/bin/env python3
"""
Generate coalescent contribution reference data using Python TreeTime v0.

This script generates snapshot files for testing the Rust implementation
of compute_coalescent_contributions() against the Python v0 baseline.

The script loads pre-computed node times from data files, then uses the
Coalescent class to compute contributions for each node.

Data generation steps (reproducible):
1. Run clock command to get rerooted tree:
   ./dev/docker/python treetime clock --tree data/flu/h3n2/20/tree.nwk \
       --dates data/flu/h3n2/20/metadata.tsv --sequence-length 1407 \
       --outdir tmp/coalescent_test
   Output: tmp/coalescent_test/rerooted.newick, clock_rate=0.002826

2. Run timetree command to get node times:
   ./dev/docker/python treetime --tree tmp/coalescent_test/rerooted.newick \
       --dates data/flu/h3n2/20/metadata.tsv --clock-rate 0.002826 \
       --sequence-len 1407 --keep-root --keep-polytomies --branch-length-mode input \
       --outdir tmp/coalescent_test/timetree
   Output: tmp/coalescent_test/timetree/dates.tsv (numeric date column)

3. Copy output files to test_scripts/data/:
   cp tmp/coalescent_test/rerooted.newick test_scripts/data/coalescent_flu_h3n2_20_tree.nwk
   cut -f1,3 tmp/coalescent_test/timetree/dates.tsv | sed 's/#node/name/' > test_scripts/data/coalescent_flu_h3n2_20_dates.tsv

4. Run this script to generate snapshots:
   ./dev/docker/python python3 test_scripts/test_coalescent_contributions.py

Input data:
    - test_scripts/data/coalescent_flu_h3n2_20_tree.nwk: Rerooted tree file
    - test_scripts/data/coalescent_flu_h3n2_20_dates.tsv: Node dates file

Output:
    - test_scripts/snapshots/coalescent_flu_h3n2_20_tc{tc}.json: Snapshot files

The snapshot files contain:
    - inputs: tree_path, metadata_path, Tc value, present_time
    - intermediate_data: tbp_grid, lineage_counts, integral_merger_rate, total_merger_rate
    - node_contributions: per internal node contribution data (n_children, tbp, contribution)
"""

import json
import sys
from pathlib import Path

import numpy as np
from Bio import Phylo

sys.path.insert(0, "packages/legacy/treetime")

from treetime.merger_models import Coalescent

DATA_DIR = Path("test_scripts/data")
TREE_PATH = DATA_DIR / "coalescent_flu_h3n2_20_tree.nwk"
DATES_PATH = DATA_DIR / "coalescent_flu_h3n2_20_dates.tsv"


def load_node_times(dates_path: Path) -> dict[str, float]:
    """Load node times from TSV file."""
    node_times = {}
    with open(dates_path) as f:
        header = f.readline()
        if not header.startswith("name"):
            raise ValueError(f"Expected header starting with 'name', got: {header}")
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            name = parts[0]
            date = float(parts[1])
            node_times[name] = date
    return node_times


def load_tree_with_times(tree_path: Path, node_times: dict[str, float]):
    """Load tree and set time_before_present on all nodes."""
    tree = Phylo.read(tree_path, "newick")

    present_time = max(node_times.values())

    for node in tree.find_clades():
        node_name = node.name
        if node_name and node_name in node_times:
            node_time = node_times[node_name]
            node.time_before_present = present_time - node_time
        else:
            node.time_before_present = None
        node.bad_branch = node.time_before_present is None

    return tree, present_time


def generate_snapshot(
    tree_path: Path,
    node_times: dict[str, float],
    tc: float,
    output_path: Path,
):
    """Generate a snapshot file for the given Tc value."""
    tree, present_time = load_tree_with_times(tree_path, node_times)
    coal = Coalescent(tree, Tc=tc)

    max_tbp = max(
        n.time_before_present
        for n in tree.find_clades()
        if n.time_before_present is not None
    )

    padding = 0.1
    tbp_grid_start = -max_tbp * padding
    tbp_grid_end = max_tbp * (1 + padding)
    tbp_grid_n_points = 101
    tbp_grid = np.linspace(tbp_grid_start, tbp_grid_end, tbp_grid_n_points)
    lineage_counts = [float(coal.nbranches(t)) for t in tbp_grid]
    integral_merger_rate = [float(coal.integral_merger_rate(t)) for t in tbp_grid]
    total_merger_rate = [float(coal.total_merger_rate(t)) for t in tbp_grid]

    node_contributions = {}
    for node in tree.find_clades():
        if node.time_before_present is None or node.bad_branch:
            continue
        name = node.name
        if node.is_terminal():
            # Leaf nodes get survival probability contribution: exp(I(t))
            # Python stores -I(t) in neg-log space.
            i_t = coal.integral_merger_rate(tbp_grid)
            # Store NegLog values: -ln(P) with P=exp(I(t)) => -I(t)
            y_neglog = (-i_t).astype(float).tolist()
        else:
            # Internal nodes: Legacy 'node_contribution' returns a distribution in NegLog space.
            # Store NegLog values directly.
            contrib = coal.node_contribution(node, tbp_grid)
            y_neglog = contrib.y.astype(float).tolist()

        node_contributions[name] = y_neglog

    result = {
        "description": f"Coalescent contributions for flu/h3n2/20 with Tc={tc}",
        "virus_path": "flu/h3n2/20",
        "inputs": {
            "tree_path": str(tree_path),
            "metadata_path": str(DATES_PATH),
            "tc": tc,
            "present_time": present_time,
        },
        "tbp_grid": {
            "start": tbp_grid_start,
            "end": tbp_grid_end,
            "padding": padding,
            "n_points": tbp_grid_n_points,
        },
        "lineage_counts": lineage_counts,
        "integral_merger_rate": integral_merger_rate,
        "total_merger_rate": total_merger_rate,
        "node_contributions": node_contributions,
    }

    with open(output_path, "w") as f:
        json.dump(result, f, indent=2)

    print(f"Generated: {output_path}")
    print(f"  Nodes with contributions: {len(node_contributions)}")


def main():
    if not TREE_PATH.exists():
        print(f"Error: Tree file not found: {TREE_PATH}", file=sys.stderr)
        print("Run the TreeTime commands in the docstring first.", file=sys.stderr)
        sys.exit(1)

    if not DATES_PATH.exists():
        print(f"Error: Dates file not found: {DATES_PATH}", file=sys.stderr)
        print("Run the TreeTime commands in the docstring first.", file=sys.stderr)
        sys.exit(1)

    node_times = load_node_times(DATES_PATH)
    print(f"Loaded {len(node_times)} node times from {DATES_PATH}")

    output_dir = Path("test_scripts/snapshots")
    output_dir.mkdir(parents=True, exist_ok=True)

    for tc in [0.01, 0.1, 1.0, 10.0]:
        output_name = f"coalescent_flu_h3n2_20_tc{tc}.json"
        generate_snapshot(TREE_PATH, node_times, tc, output_dir / output_name)

    print("\nAll snapshots generated.")


if __name__ == "__main__":
    main()
