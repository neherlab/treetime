#!/usr/bin/env python3
"""Render report charts from the hyperfine JSON into assets/ as SVG.

Charts:
  1. timetree parallel speedup (T1/Tj), before vs after, with ideal reference.
  2. timetree dengue/2000 wall clock, before vs after, grouped by thread count.
  3. speedup heatmap (command x threads), diverging.

Colours from the validated data-viz palette (categorical after=blue, before=orange;
diverging blue<->red for the heatmap). Light surface for portability on GitHub.
"""
import json
import statistics
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data" / "hyperfine"
ASSETS = HERE.parent / "assets"
ASSETS.mkdir(exist_ok=True)

# palette (light)
AFTER = "#2a78d6"     # blue, slot 1  (the new build, highlighted)
BEFORE = "#eb6834"    # orange, slot 8
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
SURFACE = "#fcfcfb"
RED = "#d03b3b"
MIDGRAY = "#f0efec"

JOBS = [1, 2, 4, 8, 16]
PHYS_CORES = 10

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "svg.fonttype": "none",
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "axes.edgecolor": "#c3c2b7",
    "text.color": INK,
    "axes.labelcolor": INK2,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
})


def load(name):
    """Return {(rev, jobs): mean_seconds} for a hyperfine JSON file."""
    d = json.loads((DATA / name).read_text())
    cell = {}
    for r in d["results"]:
        p = r["parameters"]
        rev = "before" if ("before" in p.get("rev", "") or "chore-perf-baseline" in p.get("bin", "")) else "after"
        cell[(rev, int(p["jobs"]))] = statistics.fmean(r["times"])
    return cell


def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def chart_scaling():
    c = load("timetree__dengue-2000.json")
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    x = list(range(len(JOBS)))
    # ideal linear speedup: speedup == thread count
    ax.plot(x, JOBS, color=MUTED, lw=1.2, ls=(0, (4, 3)), zorder=1)
    ax.text(x[-1], JOBS[-1], "  ideal linear", color=MUTED, va="center", fontsize=9)
    for rev, col in (("before", BEFORE), ("after", AFTER)):
        t1 = c[(rev, 1)]
        y = [t1 / c[(rev, j)] for j in JOBS]
        ax.plot(x, y, color=col, lw=2.2, marker="o", ms=7, zorder=3,
                markeredgecolor=SURFACE, markeredgewidth=1.5)
        ax.text(x[-1] + 0.12, y[-1], f"{rev}  {y[-1]:.1f}x", color=col,
                va="center", fontsize=11, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(JOBS)
    ax.set_xlim(-0.2, len(JOBS) - 1 + 1.1)
    ax.set_ylim(0, JOBS[-1] + 0.5)
    ax.set_xlabel("worker threads (--jobs)")
    ax.set_ylabel("parallel speedup  T1 / Tj")
    ax.set_title("timetree (dengue/2000): parallel speedup vs thread count",
                 color=INK, fontsize=12.5, fontweight="bold", pad=12, loc="left")
    style_axes(ax)
    fig.tight_layout()
    [fig.savefig(ASSETS / f"01-timetree-scaling.{ext}", dpi=150) for ext in ("svg", "png")]
    plt.close(fig)


def chart_runtime():
    c = load("timetree__dengue-2000.json")
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    x = list(range(len(JOBS)))
    w = 0.38
    b = [c[("before", j)] for j in JOBS]
    a = [c[("after", j)] for j in JOBS]
    bars_b = ax.bar([i - w / 2 for i in x], b, w, color=BEFORE, label="before", zorder=3)
    bars_a = ax.bar([i + w / 2 for i in x], a, w, color=AFTER, label="after", zorder=3)
    for bars, vals in ((bars_b, b), (bars_a, a)):
        for rect, v in zip(bars, vals):
            ax.text(rect.get_x() + rect.get_width() / 2, v + 0.3, f"{v:.1f}",
                    ha="center", va="bottom", fontsize=8.5, color=INK2)
    ax.set_xticks(x)
    ax.set_xticklabels(JOBS)
    ax.set_xlabel("worker threads (--jobs)")
    ax.set_ylabel("wall clock (seconds, lower is better)")
    ax.set_title("timetree (dengue/2000): wall clock roughly halves at 8-16 threads",
                 color=INK, fontsize=12, fontweight="bold", pad=12, loc="left")
    ax.legend(frameon=False, loc="upper right")
    style_axes(ax)
    fig.tight_layout()
    [fig.savefig(ASSETS / f"02-timetree-runtime.{ext}", dpi=150) for ext in ("svg", "png")]
    plt.close(fig)


def chart_heatmap():
    rows = [
        ("timetree dengue/2000", "timetree__dengue-2000.json"),
        ("timetree flu/500", "timetree__flu-h3n2-500.json"),
        ("marginal (sparse) dengue/2000", "anc-marginal__dengue-2000.json"),
        ("marginal (sparse) dengue/1000", "anc-marginal__dengue-1000.json"),
        ("marginal (dense) flu/500", "dense__flu-h3n2-500.json"),
        ("marginal (dense) dengue/1000", "dense__dengue-1000.json"),
        ("parsimony dengue/2000", "anc-parsimony__dengue-2000.json"),
        ("parsimony dengue/1000", "anc-parsimony__dengue-1000.json"),
    ]
    grid = []
    for _, fn in rows:
        c = load(fn)
        grid.append([(c[("before", j)] - c[("after", j)]) / c[("before", j)] * 100 for j in JOBS])

    cmap = LinearSegmentedColormap.from_list("div", [RED, MIDGRAY, AFTER])
    vmax = max(abs(v) for r in grid for v in r)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    ax.imshow(grid, cmap=cmap, norm=norm, aspect="auto")
    ax.set_xticks(range(len(JOBS)))
    ax.set_xticklabels(JOBS)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows], fontsize=9.5)
    ax.set_xlabel("worker threads (--jobs)")
    for i, r in enumerate(grid):
        for j, v in enumerate(r):
            ax.text(j, i, f"{v:+.0f}%", ha="center", va="center", fontsize=9,
                    color=INK if abs(v) < vmax * 0.55 else "#ffffff", fontweight="bold")
    ax.set_title("Speedup (%): blue = faster, red = slower",
                 color=INK, fontsize=12, fontweight="bold", pad=10, loc="left")
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    fig.tight_layout()
    [fig.savefig(ASSETS / f"03-speedup-heatmap.{ext}", dpi=150) for ext in ("svg", "png")]
    plt.close(fig)


if __name__ == "__main__":
    chart_scaling()
    chart_runtime()
    chart_heatmap()
    print("wrote charts to", ASSETS)
