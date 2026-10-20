# Grid search upper bound caps at 0.5 subs/site

> [!IMPORTANT]
> **Decision required.** The proposed actions are alternatives: accept the cap as a known limitation (actions 1 and 2), or extend the grid on the non-unimodal path (action 3). Extending the grid changes optimizer results on saturated edges and needs a v0 comparison before it can be chosen.

## Problem

`fn grid_search_branch_lengths()` [packages/treetime/src/optimize/zero_boundary.rs#L110-L118](../../packages/treetime/src/optimize/zero_boundary.rs#L110-L118) takes the grid upper bound from `fn search_upper_bound()` [zero_boundary.rs#L120-L122](../../packages/treetime/src/optimize/zero_boundary.rs#L120-L122), which returns `max(1.5 * branch_length + one_mutation, MIN_SEARCH_UPPER_BOUND)` with `MIN_SEARCH_UPPER_BOUND = 0.5` [packages/treetime/src/optimize/zero_boundary.rs#L11](../../packages/treetime/src/optimize/zero_boundary.rs#L11). The floor ensures the grid covers the biologically plausible range even when the current branch length estimate is zero or very small. The upper bound grows proportionally with the input, so for inputs above $0.5 / 1.5 \approx 0.33$ the proportional term dominates.

For a non-unimodal surface whose global maximum lies far beyond $1.5 \cdot \text{input}$ (e.g. a saturated substitution regime with true $t \gg 0.5$), the grid scan cannot reach it. `fn reconcile_zero_boundary()` [packages/treetime/src/optimize/zero_boundary.rs#L57-L77](../../packages/treetime/src/optimize/zero_boundary.rs#L57-L77) would then either return zero (if zero beats every visible grid point) or return the best visible positive mode, both of which are wrong in the saturated regime.

The solver intervals call the same function: `fn brent_bracket()` [packages/treetime/src/optimize/method_brent.rs#L81-L86](../../packages/treetime/src/optimize/method_brent.rs#L81-L86) and `fn solve_newton()` [packages/treetime/src/optimize/method_newton.rs#L19-L43](../../packages/treetime/src/optimize/method_newton.rs#L19-L43) end their search at `search_upper_bound()`, so neither method can move a branch beyond this bound in one step. Both solvers report a stop at the upper end where the log-likelihood still increases (`Convergence::SearchBound`) in the `branch-lengths-not-converged` run warning. The grid search has no such report.

## Evidence

The K80 surface of [Dinh and Matsen 2017](https://doi.org/10.1214/16-AAP1240) used in `test_dispatch_zero_boundary.rs` (`fn make_dinh_matsen_k80_contribution()`) has a local maximum at $t \approx 0.167$ ($\ell = -2.77333$) and its global maximum at $t \approx 6.37$ ($\ell = -2.77254$). Both values follow from the fixture coefficients (stationary points of $P(x) = p_1(x)\, p_2(x)$ with $x = e^{-t/2}$, $p_1 = 0.25 + 0.04522316\, x - 0.04545042\, x^2$, $p_2 = 0.25 - 0.04456328\, x + 0.04456328\, x^2$).

- With `one_mutation = 0.01` and a start $t_0 \leq 2$, the upper end $\max(1.5 t_0 + \text{one\_mutation}, 0.5)$ is at most $3.01$, so neither the grid nor a solver interval contains $6.37$. `test_dispatch_zero_boundary_solve_newton_never_worse_than_start_and_converges_on_dinh_matsen_k80` pins this: from $t_0 = 2$ the Newton method stops at $3.01$ with the status `SearchBound`, and only from $t_0 = 5$ (upper end $7.51$) does it reach $6.37$
- The maximum also lies outside the v0 search range $t \leq 4$ (bracket $s \in [-2, 2]$) and the fixture omits the v0 profile penalty, so finding it is not a v0 parity requirement
- Tests of the bounded recovery assert only that a positive point beats zero inside the grid; they are not tests of global optimization

## Impact

Negligible for real phylogenetic data. Branch lengths above $0.5$ subs/site are rare (correspond to near-saturation, where the likelihood surface is nearly flat and the optimizer's choice has little downstream impact). Branches above $\approx 1.0$ are typically an artifact of bad data or a misspecified model rather than a true biological signal.

On synthetic fixtures designed to stress the optimizer, and on multi-modal surfaces where the global max is at equilibrium (large $t$), the cap can produce incorrect results. A solver that stops at the upper end raises a run warning; a grid search result below the cap raises none.

## Proposed action

1. Record the cap as an intentional algorithmic limitation with explicit scope: "grid misses modes beyond $t = \max(1.5 \cdot \text{input} + \text{one\_mutation}, 0.5)$". The project bans explanatory comments in Rust source, so the carrier is the test of action 2 (its name and assertion message) and the commit message.
2. Add a test with a fixture where the true mode lies beyond the grid range, asserting that the test documents the known boundary of applicability (not that the optimizer finds the mode).
3. Alternative: extend the grid for the non-unimodal path specifically, e.g. cover $[\epsilon, \text{equilibrium\_t}]$ where equilibrium is the $\ln(2) \cdot 10$ cutoff used elsewhere in the code.

## Cross-references

- `fn grid_search_branch_lengths()`, `fn search_upper_bound()`, `const MIN_SEARCH_UPPER_BOUND` [packages/treetime/src/optimize/zero_boundary.rs](../../packages/treetime/src/optimize/zero_boundary.rs)
- `fn brent_bracket()` [packages/treetime/src/optimize/method_brent.rs](../../packages/treetime/src/optimize/method_brent.rs)
