#!/usr/bin/env python3
"""Regenerate the branch grid/tail tradeoff figures from results.csv.

    ./dev/docker/python python kb/reports/2026-08-03_branch-grid-tail-tradeoff/data/plot.py
"""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
PLOTS = HERE.parent / "plots"
PLOTS.mkdir(exist_ok=True)

BLUE, ORANGE, GREEN = "#1f77b4", "#d62728", "#2ca02c"


def load_rows():
    with open(HERE / "results.csv") as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        for k in ("grid_size", "tail_max_grid_growth", "nodes_over_tol"):
            r[k] = int(r[k])
        for k in ("tail_rel_floor", "max_abs_diff", "rms_abs_diff", "time_median_s"):
            r[k] = float(r[k])
    return rows


def by_label(rows):
    return {r["label"]: r for r in rows}


def plot_accuracy_vs_grid_size(rows):
    labels = [f"n{n}_E10_F1e-10" for n in (150, 300, 600, 1000, 2000, 3000)]
    d = by_label(rows)
    sub = [d[x] for x in labels]
    n = [r["grid_size"] for r in sub]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(n, [r["max_abs_diff"] for r in sub], "o-", color=BLUE, label="max abs diff")
    ax.plot(n, [r["rms_abs_diff"] for r in sub], "s-", color=ORANGE, label="RMS abs diff")
    ax.set_xscale("log")
    ax.set_xlabel("branch grid size n (points)")
    ax.set_ylabel("node-date error vs v0 (years)")
    ax.set_title("Accuracy vs grid size (E=10, F=1e-10)")
    ax.set_ylim(bottom=0)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(PLOTS / "accuracy_vs_grid_size.png", dpi=130)
    plt.close(fig)


def plot_runtime_vs_grid_size(rows):
    labels = [f"n{n}_E10_F1e-10" for n in (150, 300, 600, 1000, 2000, 3000)]
    d = by_label(rows)
    sub = [d[x] for x in labels]
    n = [r["grid_size"] for r in sub]
    t = [r["time_median_s"] for r in sub]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(n, t, "o-", color=BLUE, label="measured (median)")
    ref = [t[0] * (x / n[0]) ** 2 for x in n]
    ax.plot(n, ref, "--", color="gray", label="O(n^2) reference")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("branch grid size n (points)")
    ax.set_ylabel("run_timetree wall-clock (s)")
    ax.set_title("Runtime vs grid size (E=10, F=1e-10, 1 thread)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(PLOTS / "runtime_vs_grid_size.png", dpi=130)
    plt.close(fig)


def plot_accuracy_vs_tail_extent(rows):
    labels = [f"n300_E{e}_F1e-10" for e in (1, 2, 5, 10, 20)]
    d = by_label(rows)
    sub = [d[x] for x in labels]
    e = [r["tail_max_grid_growth"] for r in sub]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(e, [r["max_abs_diff"] for r in sub], "o-", color=BLUE, label="max abs diff")
    ax.plot(e, [r["rms_abs_diff"] for r in sub], "s-", color=ORANGE, label="RMS abs diff")
    ax.set_xlabel("tail extent budget E (multiple of n)")
    ax.set_ylabel("node-date error vs v0 (years)")
    ax.set_title("Accuracy vs tail extent (n=300, F=1e-10)")
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(PLOTS / "accuracy_vs_tail_extent.png", dpi=130)
    plt.close(fig)


def plot_pareto(rows):
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    for r in rows:
        ax.scatter(r["time_median_s"], r["max_abs_diff"], color=GREEN, zorder=3)
        ax.annotate(
            r["label"].replace("_F1e-10", ""),
            (r["time_median_s"], r["max_abs_diff"]),
            textcoords="offset points",
            xytext=(5, 3),
            fontsize=7,
        )
    ax.set_xscale("log")
    ax.set_xlabel("run_timetree wall-clock (s, median, 1 thread)")
    ax.set_ylabel("max node-date error vs v0 (years)")
    ax.set_title("Accuracy/cost tradeoff (flu_h3n2_20, marginal_dense)")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / "pareto_accuracy_vs_time.png", dpi=130)
    plt.close(fig)


def main():
    rows = load_rows()
    plot_accuracy_vs_grid_size(rows)
    plot_runtime_vs_grid_size(rows)
    plot_accuracy_vs_tail_extent(rows)
    plot_pareto(rows)
    print(f"wrote 4 figures to {PLOTS}")


if __name__ == "__main__":
    main()
