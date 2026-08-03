#!/usr/bin/env python3
"""Analyze the decay off/extrapolated/real comparison and regenerate its figures.

    ./dev/docker/python python kb/reports/2026-08-03_branch-grid-tail-tradeoff/data/decay_analysis.py

Reads the per-case node-date dumps under data/decay/<regime>/<label>.json and the v0 golden
master, prints the root-error and pairwise-shift tables, and writes two figures to plots/.
"""

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
PLOTS = HERE.parent / "plots"
PLOTS.mkdir(exist_ok=True)

ROOT = "NODE_0000017"
LABELS = ["decay_off", "decay_extrapolated", "decay_real"]
PRETTY = {"decay_off": "off (flat)", "decay_extrapolated": "extrapolated", "decay_real": "real"}
BLUE, ORANGE, GREEN = "#1f77b4", "#d62728", "#2ca02c"

ORACLE = json.loads(
    (REPO / "packages/treetime/src/timetree/inference/__tests__/__fixtures__/gm_runner_outputs.json").read_text()
)["flu_h3n2_20"]["marginal_dense"]


def load(regime: str) -> dict[str, dict[str, float]]:
    return {label: json.loads((HERE / "decay" / regime / f"{label}.json").read_text()) for label in LABELS}


def max_shift(a: dict[str, float], b: dict[str, float]) -> tuple[float, str, float]:
    shared = set(a) & set(b)
    diffs = [(abs(a[k] - b[k]), k) for k in shared]
    worst, node = max(diffs)
    rms = math.sqrt(sum(d * d for d, _ in diffs) / len(diffs))
    return worst, node, rms


def print_tables(reroot: dict, keeproot: dict) -> None:
    print("Root error vs v0 (reroot path):")
    for label in LABELS:
        v1 = reroot[label][ROOT]
        print(f"  {PRETTY[label]:<14} root {v1:.4f}  err {abs(v1 - ORACLE[ROOT]):.4f} y")
    for regime, data in (("reroot", reroot), ("keep-root@v0", keeproot)):
        print(f"\nNode-date shift from decay_off ({regime}):")
        for label in ("decay_extrapolated", "decay_real"):
            worst, node, rms = max_shift(data["decay_off"], data[label])
            print(f"  off -> {PRETTY[label]:<14} max {worst:.4f} y ({node}), rms {rms:.4f} y")


def plot_root_error(reroot: dict) -> None:
    errs = [abs(reroot[label][ROOT] - ORACLE[ROOT]) for label in LABELS]
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    bars = ax.bar([PRETTY[label] for label in LABELS], errs, color=[BLUE, ORANGE, GREEN])
    ax.bar_label(bars, fmt="%.3f", padding=3)
    ax.set_ylabel("root date error vs v0 (years)")
    ax.set_title("Far-tail treatment vs v0 root error (reroot path)")
    ax.set_ylim(0, max(errs) * 1.2)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / "decay_root_error_reroot.png", dpi=130)
    plt.close(fig)


def plot_shifts(reroot: dict, keeproot: dict) -> None:
    regimes = [("reroot", reroot), ("keep-root@v0", keeproot)]
    extrap = [max_shift(d["decay_off"], d["decay_extrapolated"])[0] for _, d in regimes]
    real = [max_shift(d["decay_off"], d["decay_real"])[0] for _, d in regimes]
    x = range(len(regimes))
    w = 0.35
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    b1 = ax.bar([i - w / 2 for i in x], extrap, w, label="off -> extrapolated", color=ORANGE)
    b2 = ax.bar([i + w / 2 for i in x], real, w, label="off -> real", color=GREEN)
    ax.bar_label(b1, fmt="%.3f", padding=3)
    ax.bar_label(b2, fmt="%.3f", padding=3)
    ax.set_xticks(list(x))
    ax.set_xticklabels([name for name, _ in regimes])
    ax.set_ylabel("max node-date shift from flat tail (years)")
    ax.set_title("Turning the far-tail decay on barely moves any date")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(PLOTS / "decay_shift_from_flat.png", dpi=130)
    plt.close(fig)


def main() -> None:
    reroot = load("reroot")
    keeproot = load("keeproot")
    print_tables(reroot, keeproot)
    plot_root_error(reroot)
    plot_shifts(reroot, keeproot)
    print(f"\nwrote 2 figures to {PLOTS}")


if __name__ == "__main__":
    main()
