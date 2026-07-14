# Migrate clock rerooting to generic search

Route the existing clock reroot objective through the generic `reroot` infrastructure without changing numerical behavior, defaults, or name-resolution policy.

## Required changes

- Implement `RootStats` directly for `ClockSet` using `leaf_contribution_to_parent`, `propagate_averages`, and `chisq`.
- Replace clock-only search and edge-cost calls with the corresponding generic `reroot` search and orchestration APIs.
- Extract existing per-edge `ClockSet` statistics into the map consumed by generic search.
- Preserve every caller-supplied and default `BranchPointOptimizationParams` value exactly.
- Preserve the current MinDev objective and tip-name resolution behavior; their defects are tracked separately.
- Keep the clock-specific edge fixup that swaps `to_parent` and `to_child` statistics on inverted edges.
- Update imports and delete `packages/treetime/src/clock/find_best_root/` after no callers remain.

## Validation

- For explicit grid, Brent, and golden-section parameters, compare selected edge, split fraction, score, clock model, and output before and after migration.
- Cover least-squares, min-dev, and tip reroot modes on deterministic fixtures.
- Assert all existing clock and timetree reroot tests pass without expected-value changes.
- Unit-test `RootStats` methods against the corresponding `ClockSet` methods over the same inputs.
- Full lint and test suite.

## Related issues

- Source: [kb/issues/N-reroot-clock-search-duplicates-generic-module.md](../issues/N-reroot-clock-search-duplicates-generic-module.md)
- Separate objective defect: [kb/issues/M-clock-mindev-wrong-objective.md](../issues/M-clock-mindev-wrong-objective.md)
- Separate name-policy decision: [kb/issues/N-reroot-duplicated-tip-name-resolution.md](../issues/N-reroot-duplicated-tip-name-resolution.md)
- Separate default-method decision: [kb/issues/N-reroot-split-optimizer-default-diverges-from-v0.md](../issues/N-reroot-split-optimizer-default-diverges-from-v0.md)
