#!/usr/bin/env python3
"""
Generate coalescent contribution reference data using Python TreeTime v0.

This script generates snapshot files for testing the Rust implementation
of compute_coalescent_contributions() against the Python v0 baseline.

The script uses pre-computed node times from TreeTime's timetree command,
then uses the Coalescent class to compute contributions for each node.

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

3. Run this script to generate snapshots:
   ./dev/docker/python python3 test_scripts/test_coalescent_contributions.py

Input data (from the flu/h3n2/20 dataset):
    - data/flu/h3n2/20/tree.nwk: Input phylogenetic tree
    - data/flu/h3n2/20/metadata.tsv: Sample dates

Output:
    - test_scripts/snapshots/coalescent_flu_h3n2_20_tc{tc}.json: Snapshot files
    - test_scripts/data/coalescent_flu_h3n2_20_tree.nwk: Rerooted tree file
    - test_scripts/data/coalescent_flu_h3n2_20_dates.tsv: Node dates file

The snapshot files contain:
    - inputs: tree_path, metadata_path, Tc value, present_time
    - intermediate_data: tbp_grid, lineage_counts, integral_merger_rate, total_merger_rate
    - node_contributions: per-node contribution data (is_leaf, n_children, tbp, contribution)
"""

import json
import sys
from io import StringIO
from pathlib import Path

import numpy as np
from Bio import Phylo

sys.path.insert(0, "packages/legacy/treetime")

from treetime.merger_models import Coalescent


# Clock rate estimated by TreeTime clock command
CLOCK_RATE = 0.002826

# Sequence length for the flu/h3n2/20 dataset
SEQUENCE_LENGTH = 1407

# Rerooted tree from TreeTime clock command
REROOTED_TREE = r"""((A/Scotland/76/2003|CY088128|11/03/2003|United_Kingdom|03_04|H3N2/1-1409:0.00143,A/Denmark/107/2003|EU103941|2003|Denmark||H3N2/1-1409:0.00070)NODE_00000161.00:0.01778,(A/Canterbury/58/2000|CY009150|09/05/2000|New_Zealand||H3N2/8-1416:0.00647,(A/New_York/182/2000|CY001279|02/18/2000|USA|99_00|H3N2/1-1409:0.00210,(A/Mexico/InDRE940/2003|CY100628|2003|Mexico||H3N2/15-1423:0.00578,(A/Managua/25/2007|CY032439|06/27/2007|Nicaragua||H3N2/1-1409:0.00350,(A/DaNang/DN434/2008|CY104616|11/11/2008|Viet_Nam||H3N2/4-1412:0.00428,(A/Boston/57/2008|CY044710|02/24/2008|USA|07_08|H3N2/1-1409:0.00214,((A/Oregon/15/2009|GQ895004|06/25/2009|USA|08_09|H3N2/1-1409:0.00214,A/Hong_Kong/H090_695_V10/2009|CY115546|07/10/2009|Hong_Kong||H3N2/8-1416:0.00357)NODE_00000090.98:0.00426,(A/Minab/797/2011|KC865620|12/24/2011|Iran||H3N2/20-1428:0.00761,(A/Peru/PER247/2011|CY162234|08/26/2011|Peru||H3N2/8-1416:0.00517,((A/Nebraska/15/2011|KC892583|12/15/2011|USA|11_12|H3N2/1-1409:0.00142,A/Maryland/21/2011|KC892695|12/26/2011|USA|11_12|H3N2/1-1409:0.00214)NODE_00000050.93:0.00214,(A/Indiana/03/2012|KC892731|04/03/2012|USA|11_12|H3N2/1-1409:0.00285,((A/Maryland/03/2013|KF789621|02/10/2013|USA|12_13|H3N2/1-1409:0.00500,A/New_Hampshire/12/2012|KF790252|11/08/2012|USA|12_13|H3N2/1-1409:0.00071)NODE_00000010.47:0.00055,(A/Hawaii/02/2013|KF789866|05/28/2013|USA|12_13|H3N2/1-1409:0.00357,A/Boston/DOA2_107/2012|CY148382|11/01/2012|USA|12_13|H3N2/1-1409:0.00142)NODE_00000020.45:0.00055)NODE_0000000:0.00071)NODE_00000030.81:0.00286)NODE_00000040.96:0.00486)NODE_00000060.98:0.00184)NODE_00000070.79:0.00199)NODE_00000080.81:0.00074)NODE_00000100.81:0.00071)NODE_00000110.80:0.00055)NODE_00000120.87:0.01374)NODE_00000131.00:0.01010)NODE_00000141.00:0.00066)NODE_00000150.73:0.00408)NODE_0000017:0.00100;"""

# Node times from TreeTime timetree command (numeric_date column from dates.tsv)
NODE_TIMES = {
    "NODE_0000017": 1996.940903,
    "NODE_00000161.00": 2002.844210,
    "A/Scotland/76/2003|CY088128|11/03/2003|United_Kingdom|03_04|H3N2/1-1409": 2003.840520,
    "A/Denmark/107/2003|EU103941|2003|Denmark||H3N2/1-1409": 2003.002738,
    "NODE_00000150.73": 1998.507899,
    "A/Canterbury/58/2000|CY009150|09/05/2000|New_Zealand||H3N2/8-1416": 2000.681725,
    "NODE_00000141.00": 1998.775604,
    "A/New_York/182/2000|CY001279|02/18/2000|USA|99_00|H3N2/1-1409": 2000.134155,
    "NODE_00000131.00": 2001.478893,
    "A/Mexico/InDRE940/2003|CY100628|2003|Mexico||H3N2/15-1423": 2003.002738,
    "NODE_00000120.87": 2006.569842,
    "A/Managua/25/2007|CY032439|06/27/2007|Nicaragua||H3N2/1-1409": 2007.487337,
    "NODE_00000110.80": 2006.906449,
    "A/DaNang/DN434/2008|CY104616|11/11/2008|Viet_Nam||H3N2/4-1412": 2008.865161,
    "NODE_00000100.81": 2007.219376,
    "A/Boston/57/2008|CY044710|02/24/2008|USA|07_08|H3N2/1-1409": 2008.150582,
    "NODE_00000080.81": 2007.480821,
    "NODE_00000090.98": 2008.633701,
    "A/Oregon/15/2009|GQ895004|06/25/2009|USA|08_09|H3N2/1-1409": 2009.481862,
    "A/Hong_Kong/H090_695_V10/2009|CY115546|07/10/2009|Hong_Kong||H3N2/8-1416": 2009.522929,
    "NODE_00000070.79": 2008.485779,
    "A/Minab/797/2011|KC865620|12/24/2011|Iran||H3N2/20-1428": 2011.980151,
    "NODE_00000060.98": 2009.180734,
    "A/Peru/PER247/2011|CY162234|08/26/2011|Peru||H3N2/8-1416": 2011.651608,
    "NODE_00000040.96": 2010.631447,
    "NODE_00000050.93": 2011.356438,
    "A/Nebraska/15/2011|KC892583|12/15/2011|USA|11_12|H3N2/1-1409": 2011.955510,
    "A/Maryland/21/2011|KC892695|12/26/2011|USA|11_12|H3N2/1-1409": 2011.985626,
    "NODE_00000030.81": 2011.533989,
    "A/Indiana/03/2012|KC892731|04/03/2012|USA|11_12|H3N2/1-1409": 2012.257358,
    "NODE_0000000": 2011.877473,
    "NODE_00000010.47": 2012.068368,
    "A/Maryland/03/2013|KF789621|02/10/2013|USA|12_13|H3N2/1-1409": 2013.112252,
    "A/New_Hampshire/12/2012|KF790252|11/08/2012|USA|12_13|H3N2/1-1409": 2012.856947,
    "NODE_00000020.45": 2012.146432,
    "A/Hawaii/02/2013|KF789866|05/28/2013|USA|12_13|H3N2/1-1409": 2013.405202,
    "A/Boston/DOA2_107/2012|CY148382|11/01/2012|USA|12_13|H3N2/1-1409": 2012.837782,
}


def load_tree_with_times():
    """Load the rerooted tree and set time_before_present on all nodes."""
    tree = Phylo.read(StringIO(REROOTED_TREE), "newick")

    # Present time is the most recent sample date
    present_time = max(NODE_TIMES.values())

    # Set time_before_present on all nodes
    for node in tree.find_clades():
        node_name = node.name
        if node_name and node_name in NODE_TIMES:
            node_time = NODE_TIMES[node_name]
            node.time_before_present = present_time - node_time
        else:
            node.time_before_present = None
        node.bad_branch = node.time_before_present is None

    return tree, present_time


def compute_leaf_contribution(coal: Coalescent, tbp: float) -> dict:
    """Compute contribution for a leaf node.

    For a leaf at time t, the contribution to the coalescent likelihood is:
    exp(I(t)) where I(t) is the integral of the merger rate from present to t.

    This accounts for the probability of no coalescence happening before this
    lineage appears.
    """
    i_t = float(coal.integral_merger_rate(tbp))
    neg_log = -i_t
    prob = np.exp(-neg_log)  # = exp(i_t)
    # Handle overflow
    if np.isinf(prob):
        prob = 1e308
    return {
        "integral_merger_rate": i_t,
        "neg_log": neg_log,
        "prob": prob,
    }


def compute_internal_contribution(coal: Coalescent, node, tbp: float) -> dict:
    """Compute contribution for an internal node.

    For an internal node with m = n_children - 1 mergers at time t:
    neg_log_prob = m * (I(t) - log(λ(t)))

    where:
    - I(t) is the integral of the merger rate
    - λ(t) is the total merger rate at time t
    - m is the number of mergers (n_children - 1)
    """
    i_t = float(coal.integral_merger_rate(tbp))
    lambda_t = float(coal.total_merger_rate(tbp))
    m = len(node.clades) - 1
    neg_log_prob = m * (i_t - np.log(lambda_t))
    prob = np.exp(-neg_log_prob)
    return {
        "integral_merger_rate": i_t,
        "total_merger_rate": lambda_t,
        "multiplicity": m,
        "neg_log": neg_log_prob,
        "prob": prob,
    }


def generate_snapshot(tc: float, output_path: Path):
    """Generate a snapshot file for the given Tc value."""
    tree, present_time = load_tree_with_times()
    coal = Coalescent(tree, Tc=tc)

    node_contributions = {}
    for node in tree.find_clades():
        if node.time_before_present is None or node.bad_branch:
            continue

        name = node.name
        is_leaf = node.is_terminal()
        tbp = node.time_before_present
        n_children = len(node.clades)

        if is_leaf:
            contrib = compute_leaf_contribution(coal, tbp)
        else:
            contrib = compute_internal_contribution(coal, node, tbp)

        # Add tbp to contribution dict for Rust test compatibility
        contrib["tbp"] = tbp

        node_contributions[name] = {
            "is_leaf": is_leaf,
            "n_children": n_children,
            "time_before_present": tbp,
            "contribution": contrib,
        }

    max_tbp = max(n.time_before_present for n in tree.find_clades() if n.time_before_present is not None)

    tbp_grid = np.linspace(0.0, max_tbp * 1.1, 101).tolist()
    lineage_counts = [float(coal.nbranches(t)) for t in tbp_grid]
    integral_merger_rate = [float(coal.integral_merger_rate(t)) for t in tbp_grid]
    total_merger_rate = [float(coal.total_merger_rate(t)) for t in tbp_grid]

    result = {
        "description": f"Coalescent contributions for flu/h3n2/20 with Tc={tc}",
        "inputs": {
            "tree_path": "test_scripts/data/coalescent_flu_h3n2_20_tree.nwk",
            "metadata_path": "test_scripts/data/coalescent_flu_h3n2_20_dates.tsv",
            "tc": tc,
            "present_time": present_time,
        },
        "intermediate_data": {
            "tbp_grid": tbp_grid,
            "lineage_counts": lineage_counts,
            "integral_merger_rate": integral_merger_rate,
            "total_merger_rate": total_merger_rate,
        },
        "node_contributions": node_contributions,
    }

    with open(output_path, "w") as f:
        json.dump(result, f, indent=2)

    n_leaves = sum(1 for n in node_contributions.values() if n["is_leaf"])
    n_internal = sum(1 for n in node_contributions.values() if not n["is_leaf"])
    print(f"Generated: {output_path}")
    print(f"  Nodes: {len(node_contributions)} (leaves: {n_leaves}, internal: {n_internal})")


def main():
    output_dir = Path("test_scripts/snapshots")
    output_dir.mkdir(parents=True, exist_ok=True)

    data_dir = Path("test_scripts/data")
    data_dir.mkdir(parents=True, exist_ok=True)

    # Save tree file
    tree_path = data_dir / "coalescent_flu_h3n2_20_tree.nwk"
    with open(tree_path, "w") as f:
        f.write(REROOTED_TREE)
    print(f"Saved tree: {tree_path}")

    # Save dates file (TSV format: name, date)
    dates_path = data_dir / "coalescent_flu_h3n2_20_dates.tsv"
    with open(dates_path, "w") as f:
        f.write("name\tdate\n")
        for name, time in NODE_TIMES.items():
            f.write(f"{name}\t{time:.6f}\n")
    print(f"Saved dates: {dates_path}")

    for tc in [0.01, 0.1, 1.0, 10.0]:
        output_name = f"coalescent_flu_h3n2_20_tc{tc}.json"
        generate_snapshot(tc, output_dir / output_name)

    print("\nAll snapshots generated.")


if __name__ == "__main__":
    main()
