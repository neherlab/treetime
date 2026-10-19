# Timetree reports one weak-signal failure with different errors

> [!IMPORTANT]
> **Investigation required.** Only one of the three errors is recorded. From the failing runs of the timetree suite in [validation/README.md](../../validation/README.md) (the `error` column of `results/timetree.tsv`), collect the three messages, the inputs that produce them, and the code path of each.

## Problem

On simulated data with weak temporal signal, the same case fails with three different errors, depending on the code path. One of them is internal: `Cannot divide by empty distribution` [packages/treetime-distribution/src/distribution_ops/divide.rs#L30](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L30). None of the messages tells the user that the date signal is weak.

## Related errors

> [!WARNING]
> **Needs review.** The recorded reproductions of the first two errors (`data/dengue/500`, `data/dengue/1000`, `data/lassa/L/20`, and `keep-skyline-fixed/const/rep01` at $L = 10^4$ of the timetree suite) complete at commit `dd2eb7428`. Collect the current errors of the suite before tracing the paths.

- `Cannot divide by empty distribution`: division by an empty distribution
- a `NaN` time distribution, raised with and without round and node context
- an undated node, reported as a stale coalescent model by the point-time lineage counts, which phase 2 of [timetree-coalescent-prior-and-clock-rate-streamlining.md](../proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md) removed; the expected lineage counts count an undated node as surviving at all times

## Fix direction

Unknown until the paths are traced. A node whose date the rest of the tree gives no probability gets an empty posterior; the forward pass keeps the date the node was given and warns ([packages/treetime/src/timetree/inference/forward_pass.rs#L43-L53](../../packages/treetime/src/timetree/inference/forward_pass.rs#L43-L53)).

## Validation

Each weak-signal input of the validation suite ends with the same result on every code path: one error that names the cause, or a completed run.
