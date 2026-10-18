# Timetree reports one weak-signal failure with different errors

> [!IMPORTANT]
> **Investigation required.** Only one of the three errors is recorded. From the failing runs of the timetree suite in [validation/README.md](../../validation/README.md) (the `error` column of `results/timetree.tsv`), collect the three messages, the inputs that produce them, and the code path of each.

## Problem

On simulated data with weak temporal signal, the same case fails with three different errors, depending on the code path. One of them is internal: `Cannot divide by empty distribution` [packages/treetime-distribution/src/distribution_ops/divide.rs#L30](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L30). None of the messages tells the user that the date signal is weak.

## Related errors

- [H-timetree-empty-distribution-division.md](H-timetree-empty-distribution-division.md): division by an empty distribution on real datasets
- [H-timetree-backward-pass-nan-time-distribution.md](H-timetree-backward-pass-nan-time-distribution.md): a `NaN` time distribution, raised with and without round and node context
- an undated node, reported as a stale coalescent model by the point-time lineage counts, which phase 2 of [timetree-coalescent-prior-and-clock-rate-streamlining.md](../proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md) removed; the expected lineage counts count an undated node as surviving at all times

## Fix direction

Unknown until the three paths are traced. What an empty posterior should resolve to is described in the Mechanism section of [H-timetree-backward-pass-nan-time-distribution.md](H-timetree-backward-pass-nan-time-distribution.md).

## Validation

Each weak-signal input of the validation suite ends with the same result on every code path: one error that names the cause, or a completed run.
