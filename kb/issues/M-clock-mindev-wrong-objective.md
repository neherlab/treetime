# MinDev reroot uses wrong objective (EstimatedRate instead of FixedRate(0))

v1's MinDev reroot uses `RootObjective::EstimatedRate` (freely estimated clock rate) instead of `RootObjective::FixedRate(0.0)` (zero rate). This minimizes the full WLS regression residual rather than root-to-tip distance variance.

## v0 behavior

`treetime.py:575`: `slope = 0.0 if ... root.startswith('min_dev') ...`

v0 passes `slope=0` to `optimal_reroot` -> `find_best_root(slope=0)`. With slope=0, `chisq_fixed_rate(0)` reduces to:

$(S_{dd} \cdot n - S_d^2) / (2n^2)$

This is proportional to the variance of root-to-tip distances. The v0 docstring confirms: `"'min_dev' minimizes variance of root-to-tip distances"` (`argument_parser.py:83`).

## v1 behavior

`clock/reroot.rs:174-177`:

```rust
RerootSpec::Method(RerootMethod::MinDev) => {
    find_best_root(graph, options, params, false, RootObjective::EstimatedRate)
}
```

`EstimatedRate` calls `ClockSet::chisq()`, which minimizes the freely estimated-rate regression residual -- a different objective that factors out the best-fit rate at each candidate root.

## Impact

MinDev finds a different root than v0 on the same tree. The difference is small for well-clocked trees (where the freely estimated rate is close to zero at the optimal root) but can be significant for trees with strong temporal signal or asymmetric branch-length distributions.

## Fix

Change `RootObjective::EstimatedRate` to `RootObjective::FixedRate(0.0)` in the MinDev branch of `select_root` (`clock/reroot.rs:174-177`).

## Locations

- `packages/treetime/src/clock/reroot.rs:174-177`
- v0 reference: `packages/legacy/treetime/treetime/treetime.py:575`
- v0 docstring: `packages/legacy/treetime/treetime/argument_parser.py:83`
