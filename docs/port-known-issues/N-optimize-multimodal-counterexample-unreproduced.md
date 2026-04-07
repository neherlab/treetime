# Multi-modal surface counterexample with $\ell'(0) < 0$ is not constructed

## Problem

`reconcile_zero_boundary` handles two failure modes of the per-edge inner solver on non-unimodal models: a positive candidate worse than zero, and an exact-zero candidate from Newton step clamping. The theoretical failure mode the helper is designed to catch is a multi-modal surface where:

1. $\ell'(0) < 0$ (zero is a local maximum), AND
2. The inner solver returns exactly $0$ (observed for `newton_sqrt_inner` from $t_0 = 0.6$ on the Dinh and Matsen 2017 K80 $\kappa = 3$ counterexample), AND
3. A better positive mode exists elsewhere on the admissible interval.

The Dinh-Matsen K80 counterexample at [packages/treetime/src/commands/optimize/**tests**/test_is_zero_branch_optimal.rs#L452-L500](../../packages/treetime/src/commands/optimize/__tests__/test_is_zero_branch_optimal.rs#L452-L500) satisfies (2) and (3) but not (1): on that surface $\ell'(0) > 0$, so the function is INCREASING at the boundary rather than exhibiting zero as a local maximum. The observed sqrt-space clamp at $t_0 = 0.6$ is caused by chain rule sign flip, not by zero being a local max.

No surface in the v1 codebase satisfies all three conditions simultaneously. Constructing one from a real GTR model (K80, HKY85, TN93, general GTR) on a real alignment has not been exhibited. A synthetic contribution with hand-crafted coefficients could in principle be constructed, but the construction itself is open research.

## Impact

Negligible under the current codebase. The exact-zero gate in `reconcile_zero_boundary` is defense in depth: it correctly routes exact-zero candidates through grid verification for any non-unimodal partition set, whether or not the surface actually exhibits the (1)+(2)+(3) combination.

Without a real counterexample, the gate's empirical justification rests on the observed sqrt-space clamp and the theoretical argument that grid verification cannot make the wrong decision (it either confirms zero or finds a better positive mode).

## Proposed action

1. Attempt to construct a real-data counterexample where all three conditions hold simultaneously. Candidates: K80 with high $\kappa$, HKY85 with skewed frequencies, TN93 with distinct transition/transversion rates, on alignments with specific mutation patterns.
2. If a real counterexample is identified, add a regression test to [`test_dispatch_zero_boundary.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_dispatch_zero_boundary.rs) that exercises the full `run_optimize_mixed` flow on it and asserts `reconcile_zero_boundary` returns the correct positive mode.
3. If no real counterexample is found after exhaustive search, document the result and upgrade the gate's justification from empirical to theoretical.

## v0 comparison

v0 uses `scipy.optimize.minimize_scalar(method='brent')` in $\sqrt{t}$ space and does not have a structural step-clamping pattern that can reach exact zero. The (1)+(2)+(3) scenario is specific to v1's Newton variants.

## Cross-references

- Helper: [packages/treetime/src/commands/optimize/optimize_unified.rs](../../packages/treetime/src/commands/optimize/optimize_unified.rs) (`reconcile_zero_boundary`)
- Counterexample fixture: [packages/treetime/src/commands/optimize/**tests**/test_dispatch_zero_boundary.rs](../../packages/treetime/src/commands/optimize/__tests__/test_dispatch_zero_boundary.rs) (`make_dinh_matsen_k80_contribution`)
- Scientific basis: Dinh and Matsen 2017, "The Shape of the One-Dimensional Phylogenetic Likelihood Function", <https://doi.org/10.1214/16-aap1240>
