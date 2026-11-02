# Negated tail laws carry a negative zero

`DistributionFunction::negate_arg_inplace()` negates the tail laws of both ends with `fn scale_tail_law(..., -1.0)` ([packages/treetime-distribution/src/distribution_core/function.rs#L333-L342](../../packages/treetime-distribution/src/distribution_core/function.rs#L333-L342), [#L416-L423](../../packages/treetime-distribution/src/distribution_core/function.rs#L447-L454)). The laws multiply their coordinate by the factor: `HardApproachLaw::scale_arg()` computes `t_hard * factor` ([packages/treetime-grid/src/hard_approach_law.rs#L81-L86](../../packages/treetime-grid/src/hard_approach_law.rs#L81-L86)), and `EdgeCellLaw::scale_arg()` computes `t_edge * factor` ([packages/treetime-grid/src/edge_cell_law.rs#L14-L20](../../packages/treetime-grid/src/edge_cell_law.rs#L14-L20)). A hard end at 0, such as the end of every branch-time distribution at branch length 0, becomes `-0.0` after negation: `0.0 * -1.0` is `-0.0` in IEEE 754 arithmetic.

`-0.0` compares equal to `0.0`, so the sign does not change a comparison or a sum with a nonzero value. It shows in debug output (`t_hard: -0.0`) and in every operation that reads the sign (`is_sign_negative()`, division by the value, formatting). The negation of the grid of the same function avoids the sign by subtracting from zero (`Grid::negated()`, [packages/treetime-grid/src/grid.rs#L188-L193](../../packages/treetime-grid/src/grid.rs#L188-L193)), so a negated branch distribution ends its grid at `0.0` and its tail law at `-0.0`.

## Fix direction

Negate the coordinates of the tail laws by subtraction from zero (`0.0 - t_hard`), or give the laws a `negate_arg()` that does so, so the grid and its tail laws agree on the sign of a zero end.
