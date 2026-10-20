# Timetree reports one weak-signal failure with different errors

> [!IMPORTANT]
> **Investigation required.** Only one of the errors is recorded with its code path. The validation results that list the failing runs are not in the checkout. Produce them with the timetree suite of [validation/README.md](../../validation/README.md): build the binary with `./dev/docker/run just build-release treetime`, then in `validation/` run `pixi run python run.py --suite timetree --set core` and `pixi run python evaluate.py --suite timetree`, which writes `results/timetree.tsv` with the TreeTime error message of each failed run in its `error` column. From that table, collect the messages, the inputs that produce them, and the code path of each. The recorded reproductions of the `NaN` and undated-node errors (`data/dengue/500`, `data/dengue/1000`, `data/lassa/L/20`, and `keep-skyline-fixed/const/rep01` at $L = 10^4$ of the timetree suite) complete at commit `dd2eb7428`, so these two errors need new reproductions.

## Problem

On simulated data with weak temporal signal, the same case fails with different errors, depending on the code path. One of them is internal: `Cannot divide by empty distribution` [packages/treetime-distribution/src/distribution_ops/divide.rs#L30](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L30), raised by the division that forms the cavity of a node in the forward pass [packages/treetime/src/timetree/inference/forward_pass.rs#L153](../../packages/treetime/src/timetree/inference/forward_pass.rs#L153). The same division raises a second internal error for two point distributions [packages/treetime-distribution/src/distribution_ops/divide.rs#L31-L34](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L31-L34). None of the messages tells the user that the date signal is weak.

## Related errors

- `Cannot divide by empty distribution`: division by an empty distribution
- `Cannot divide point by point: operation not well-defined`: division of two point distributions. Every run with `--branch-length-mode=input` fails with it, see [H-timetree-input-branch-lengths-abort-on-point-division.md](H-timetree-input-branch-lengths-abort-on-point-division.md)
- a `NaN` time distribution, raised with and without round and node context
- an undated node, reported as a stale coalescent model by the point-time lineage counts, which phase 2 of [timetree-coalescent-prior-and-clock-rate-streamlining.md](../proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md) removed; the expected lineage counts count an undated node as surviving at all times

## Fix direction

Unknown until the paths are traced. A node whose date the rest of the tree gives no probability gets an empty posterior; the forward pass keeps the date the node was given and warns ([packages/treetime/src/timetree/inference/forward_pass.rs#L43-L53](../../packages/treetime/src/timetree/inference/forward_pass.rs#L43-L53)).

## Validation

Each weak-signal input of the validation suite ends with the same result on every code path: one error that names the cause, or a completed run.
