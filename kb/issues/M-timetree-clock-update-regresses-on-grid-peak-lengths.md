# In-loop clock update regresses on grid-peak lengths scaled by gamma

Inside the refinement loop, v1 re-estimates the clock model from edge divergences computed as `time_length * rate * gamma` (`fn edge_divergence` in [packages/treetime/src/clock/clock_regression.rs](../../packages/treetime/src/clock/clock_regression.rs), fed by `Refinement::update_clock_model` in [packages/treetime/src/timetree/refinement.rs](../../packages/treetime/src/timetree/refinement.rs)). `time_length` is the peak of the per-edge branch-length likelihood on a 300-point time grid, so this value is the ML branch length plus grid error, multiplied by the relaxed-clock gamma.

v0 regresses on `mutation_length`, the ML branch length itself ([packages/legacy/treetime/treetime/clock_tree.py#L276](../../packages/legacy/treetime/treetime/clock_tree.py#L276)), with no gamma factor.

## Impact

- Grid quantization adds noise to the root-to-tip divergences that the clock rate is fitted from (magnitude not measured)
- With `--relax`, the regression input is scaled by the fitted rate multipliers, which v0 does not do
