#!/usr/bin/env python3
"""Aggregate hyperfine JSON + perf-stat CSV into comparison statistics.

Reads:
  <bench_dir>/hyperfine/*.json   one file per (command,dataset), each holding
                                 before/after x jobs parameter combos.
  <bench_dir>/perfstat/*.csv     optional; perf stat -x, output per single run.

Writes:
  <bench_dir>/derived/summary.json
  <bench_dir>/derived/runtime_table.md
  <bench_dir>/derived/scaling_table.md
  <bench_dir>/derived/hwcounters_table.md   (if perfstat present)
"""
import glob
import json
import math
import os
import statistics
import sys
from collections import defaultdict

RNG_SEED = 20260729


def bootstrap_ratio_ci(before, after, n_boot=10000, alpha=0.05, seed=RNG_SEED):
    """Bootstrap CI for ratio before_mean/after_mean (speedup > 1 = after faster)."""
    import random
    rng = random.Random(seed)
    nb, na = len(before), len(after)
    ratios = []
    for _ in range(n_boot):
        b = sum(before[rng.randrange(nb)] for _ in range(nb)) / nb
        a = sum(after[rng.randrange(na)] for _ in range(na)) / na
        ratios.append(b / a)
    ratios.sort()
    lo = ratios[int((alpha / 2) * n_boot)]
    hi = ratios[int((1 - alpha / 2) * n_boot)]
    return lo, hi


def stats_of(times):
    n = len(times)
    mean = statistics.fmean(times)
    sd = statistics.stdev(times) if n > 1 else 0.0
    med = statistics.median(times)
    mad = statistics.median([abs(t - med) for t in times]) if n > 1 else 0.0
    return {
        "n": n,
        "mean": mean,
        "stddev": sd,
        "cv": (sd / mean) if mean else 0.0,
        "median": med,
        "mad": mad,
        "min": min(times),
        "max": max(times),
    }


def load_hyperfine(bench_dir):
    # cell[(cmd,ds,rev,jobs)] = {"times":[...], ...}
    cells = {}
    for path in sorted(glob.glob(os.path.join(bench_dir, "hyperfine", "*.json"))):
        base = os.path.basename(path)[:-5]  # strip .json
        # filename convention: <cmd>__<ds-with-dashes>.json
        cmd, ds = base.split("__", 1)
        ds = ds.replace("-", "/")
        with open(path) as f:
            data = json.load(f)
        for r in data["results"]:
            p = r.get("parameters", {})
            rev = p.get("rev")
            jobs = int(p.get("jobs"))
            cells[(cmd, ds, rev, jobs)] = r["times"]
    return cells


def load_perfstat(bench_dir):
    # returns counters[(cmd,ds,rev,jobs)] = {event: value}
    counters = defaultdict(dict)
    pdir = os.path.join(bench_dir, "perfstat")
    if not os.path.isdir(pdir):
        return counters
    for path in sorted(glob.glob(os.path.join(pdir, "*.csv"))):
        base = os.path.basename(path)[:-4]
        # filename: <cmd>__<ds-dashes>__<rev>__j<jobs>.csv
        parts = base.split("__")
        cmd, ds, rev, jtok = parts[0], parts[1].replace("-", "/"), parts[2], parts[3]
        jobs = int(jtok[1:])
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                fields = line.split(",")
                # perf stat -x, format: value,unit,event,run_time,pct,...
                if len(fields) < 3:
                    continue
                val_s, _unit, event = fields[0], fields[1], fields[2]
                try:
                    val = float(val_s)
                except ValueError:
                    continue
                counters[(cmd, ds, rev, jobs)][event] = val
    return counters


def fmt(x, nd=4):
    return f"{x:.{nd}f}"


def main():
    bench_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    cells = load_hyperfine(bench_dir)
    counters = load_perfstat(bench_dir)
    derived = os.path.join(bench_dir, "derived")
    os.makedirs(derived, exist_ok=True)

    keys = sorted({(c, d) for (c, d, _, _) in cells})
    jobset = sorted({j for (_, _, _, j) in cells})
    summary = {}

    # Runtime + ratio table
    rt_lines = [
        "| command | dataset | jobs | before mean (s) | after mean (s) | before CV | after CV | speedup (before/after) | 95% CI | % change |",
        "|---|---|---|---:|---:|---:|---:|---:|---|---:|",
    ]
    for (cmd, ds) in keys:
        for j in jobset:
            b = cells.get((cmd, ds, "before", j))
            a = cells.get((cmd, ds, "after", j))
            if not b or not a:
                continue
            sb, sa = stats_of(b), stats_of(a)
            speedup = sb["mean"] / sa["mean"]
            lo, hi = bootstrap_ratio_ci(b, a)
            pct = (sa["mean"] - sb["mean"]) / sb["mean"] * 100.0
            summary[f"{cmd}|{ds}|{j}"] = {
                "before": sb, "after": sa,
                "speedup_before_over_after": speedup,
                "speedup_ci95": [lo, hi],
                "pct_change_after_vs_before": pct,
            }
            rt_lines.append(
                f"| {cmd} | {ds} | {j} | {fmt(sb['mean'])} | {fmt(sa['mean'])} | "
                f"{fmt(sb['cv'],3)} | {fmt(sa['cv'],3)} | {fmt(speedup,3)} | "
                f"[{fmt(lo,3)}, {fmt(hi,3)}] | {fmt(pct,1)}% |"
            )

    # Scaling table (self-relative speedup T1/Tj and efficiency)
    sc_lines = [
        "| command | dataset | rev | jobs | mean (s) | speedup T1/Tj | efficiency |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for (cmd, ds) in keys:
        for rev in ("before", "after"):
            base1 = cells.get((cmd, ds, rev, 1))
            t1 = statistics.fmean(base1) if base1 else None
            for j in jobset:
                cj = cells.get((cmd, ds, rev, j))
                if not cj:
                    continue
                tj = statistics.fmean(cj)
                sp = (t1 / tj) if t1 else float("nan")
                eff = sp / j if t1 else float("nan")
                sc_lines.append(
                    f"| {cmd} | {ds} | {rev} | {j} | {fmt(tj)} | {fmt(sp,3)} | {fmt(eff,3)} |"
                )

    with open(os.path.join(derived, "runtime_table.md"), "w") as f:
        f.write("\n".join(rt_lines) + "\n")
    with open(os.path.join(derived, "scaling_table.md"), "w") as f:
        f.write("\n".join(sc_lines) + "\n")

    # Performance counters table
    if counters:
        hw_lines = [
            "| command | dataset | rev | jobs | IPC | cache-miss % | branch-miss % | task-clock (s) |",
            "|---|---|---|---:|---:|---:|---:|---:|",
        ]
        for (cmd, ds, rev, j) in sorted(counters):
            c = counters[(cmd, ds, rev, j)]
            ins = c.get("instructions")
            cyc = c.get("cycles")
            cref = c.get("cache-references")
            cmis = c.get("cache-misses")
            bins = c.get("branches") or c.get("branch-instructions")
            bmis = c.get("branch-misses")
            taskclock_ms = c.get("task-clock") or c.get("task-clock:u")
            ipc = (ins / cyc) if (ins and cyc) else float("nan")
            cmr = (cmis / cref * 100) if (cmis and cref) else float("nan")
            bmr = (bmis / bins * 100) if (bmis and bins) else float("nan")
            tc_s = (taskclock_ms / 1000.0) if taskclock_ms else float("nan")
            hw_lines.append(
                f"| {cmd} | {ds} | {rev} | {j} | {fmt(ipc,3)} | {fmt(cmr,3)} | "
                f"{fmt(bmr,3)} | {fmt(tc_s,3)} |"
            )
        with open(os.path.join(derived, "hwcounters_table.md"), "w") as f:
            f.write("\n".join(hw_lines) + "\n")

    with open(os.path.join(derived, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print(f"wrote derived tables to {derived}")
    print(f"cells: {len(cells)}  perfstat cells: {len(counters)}")


if __name__ == "__main__":
    main()
