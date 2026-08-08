# Timetree convergence metric deficiencies

## Summary

Make the reported total log-likelihood comparable across iterations, then decide whether and how it
enters the convergence decision. Fix sequence-diff undercounting across topology changes.

The convergence criterion itself has already been changed to node-time movement and is out of scope
here; see
[timetree-convergence-on-node-times.md](../decisions/timetree-convergence-on-node-times.md).

## Scope

### 1. Make `log_lh_total` a fixed composition

[`metrics.rs`](../../packages/treetime/src/timetree/convergence/metrics.rs) sums whichever of
`log_lh_seq`, `log_lh_pos`, `log_lh_coal` are `Some`. Decide the contract: either the total is
`None` unless every component is present, or the component set is fixed per run and a missing
component is an error. A component silently dropping out mid-run must not look like a likelihood
improvement.

Note that `log_lh_coal` goes missing for two different reasons — the run has no coalescent, or
`collect_coalescent_edges` failed on an inverted edge — and only the second is a defect. Coordinate
with
[M-coalescent-edge-collection-nan-bypass-and-unreachable-fallback.md](../issues/M-coalescent-edge-collection-nan-bypass-and-unreachable-fallback.md).

### 2. Decide whether likelihood enters `has_converged`

Requires a decision first, because the coalescent term is evaluated against live lineage counts
while the times were inferred under frozen ones, so the reported total is not the maximized
objective. Options: test only the components that are part of the objective; test the total as a
plateau guard rather than a criterion; or leave the criterion on node times and treat likelihood as
reporting only. Do not implement a stop rule until this is settled.

### 3. Count sequence changes across topology changes

[`sequence_changes.rs`](../../packages/treetime/src/timetree/convergence/sequence_changes.rs)
should attribute diffs for nodes present in only one snapshot rather than logging and discarding
them.

## Validation

`data/ebola/20` and `data/mpox/clade-ii/1000`, with and without `--coalescent-opt` and
`--resolve-polytomies`, comparing the tracelog CSV before and after. Round counts must not regress.

## Related issues

- Source: [kb/issues/M-timetree-convergence-metric-deficiencies.md](../issues/M-timetree-convergence-metric-deficiencies.md) -- delete after full resolution
