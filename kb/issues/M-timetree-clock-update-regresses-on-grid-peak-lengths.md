# In-loop clock update regresses on grid-peak branch lengths

Inside the refinement loop, `fn update_clock_model()` ([packages/treetime/src/timetree/round.rs#L248-L300](../../packages/treetime/src/timetree/round.rs#L248-L300)) re-estimates the clock model. The regression takes the length of each edge from `fn edge_length()` ([packages/treetime/src/clock/clock_regression.rs#L326-L337](../../packages/treetime/src/clock/clock_regression.rs#L326-L337)): `time_length * rate * gamma` when the edge has a time length, and the branch length of the round otherwise.

`time_length` is the peak of the per-edge branch-length likelihood on its time grid of `GRID_POINTS` (300) points ([packages/treetime/src/timetree/inference/runner.rs#L29](../../packages/treetime/src/timetree/inference/runner.rs#L29), [packages/treetime/src/timetree/inference/runner.rs#L285](../../packages/treetime/src/timetree/inference/runner.rs#L285)). The grid converts branch length to time with the rate `clock_rate * gamma` ([packages/treetime/src/timetree/inference/runner.rs#L274-L282](../../packages/treetime/src/timetree/inference/runner.rs#L274-L282)), so the product `time_length * rate * gamma` gives back the branch length at the grid peak. The gamma factor cancels, and `--relax` does not scale the regression input. What remains is the grid: the regression input is the most likely branch length plus the quantization error of the grid, not the branch length itself.

v0 regresses on `mutation_length`, the ML branch length itself ([packages/legacy/treetime/treetime/clock_tree.py#L276](../../packages/legacy/treetime/treetime/clock_tree.py#L276)). The v1 counterpart is the branch length of the round in `RoundState::branch_lengths`, which `fn edge_length()` uses only for edges without a time length.

> [!IMPORTANT]
> **Investigation required.** The size of the effect is not measured. Run `timetree` on several datasets (for example `data/ebola/{20,100,362}` and `data/flu/h3n2/{20,200}`) twice: once as is, and once with `fn edge_length()` returning the branch length of the round. Compare the clock rate after each round, the final clock rate, and the node dates. The issue is negligible if the rate differs by much less than its reported standard deviation.

## Impact

- Grid quantization adds noise to the root-to-tip divergences that the clock rate is fitted from
- The regression input differs from the v0 input by the grid error on every edge with a time length
