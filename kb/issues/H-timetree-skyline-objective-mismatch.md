# Skyline optimizer uses a different objective than the reported coalescent cost

`skyline::compute_total_neg_log_lh` at [packages/treetime/src/coalescent/skyline.rs#L261-L295](../../packages/treetime/src/coalescent/skyline.rs#L261-L295) iterates lineage-count breakpoints and adds one `-ln(lambda)` per merger event. When a polytomy at a single breakpoint merges `m` lineages simultaneously (`k(t1) < k(t0)` with `delta > 1`), the cost function adds only one `-ln(lambda)` term. The Kingman coalescent requires `(m-1)` merger factors for an `m`-merger event.

The authoritative per-edge cost function `sum_coalescent_cost` at [packages/treetime/src/coalescent/edge_data.rs#L112-L146](../../packages/treetime/src/coalescent/edge_data.rs#L112-L146) correctly uses per-edge contributions and is used by `optimize_tc` and `compute_coalescent_total_lh`. The skyline optimizer maximizes a different quantity than the one reported and cross-compared in convergence metrics.

## Details

`compute_total_neg_log_lh` detects mergers by checking `k(t1) < k(t0)` and computes `k_clamped = max(0.5, k(t1))`, recovering the pre-event pair count only for binary mergers. For polytomies, a single `delta < 0` breakpoint absorbs multiple simultaneous mergers but the cost adds only one `-ln(lambda)` term.

Skyline log-likelihood reconstruction at [packages/treetime/src/coalescent/skyline.rs#L134-L144](../../packages/treetime/src/coalescent/skyline.rs#L134-L144) mixes clamped and unclamped parameters: the likelihood and smoothness terms use `log_tc_clamped` while the boundary penalty uses unclamped `log_tc`. The subtraction to recover `neg_log_lh` from `final_cost` is correct only when `best_param` stays inside the clamping range `[-200, 100]`.

## Impact

- `optimize_skyline()` finds parameters that minimize an incorrect objective
- The reported `SkylineResult.log_likelihood` does not equal the coalescent likelihood that `sum_coalescent_cost` would compute for the same Tc values
- Trees with polytomies are disproportionately affected: the skyline underpays the cost of simultaneous mergers, biasing Tc toward values that incorrectly explain polytomies as rapid coalescence rather than unresolved topology

## Affected code

- Skyline cost: [packages/treetime/src/coalescent/skyline.rs#L261-L295](../../packages/treetime/src/coalescent/skyline.rs#L261-L295)
- Clamped/unclamped mixing: [packages/treetime/src/coalescent/skyline.rs#L134-L144](../../packages/treetime/src/coalescent/skyline.rs#L134-L144), [skyline.rs#L173-L197](../../packages/treetime/src/coalescent/skyline.rs#L173-L197)
- Authoritative per-edge cost: [packages/treetime/src/coalescent/edge_data.rs#L112-L146](../../packages/treetime/src/coalescent/edge_data.rs#L112-L146)

## Fix

Replace the merger-detection loop in `compute_total_neg_log_lh` with a formulation that counts `(m-1)` factors per merger event, matching the per-edge Kingman algebra used by `sum_coalescent_cost`. Unify the clamped/unclamped parameter handling so that the optimizer, the reported likelihood, and the per-edge cost all use the same convention.

## Related

- [M-timetree-skyline-nelder-mead-optimizer.md](M-timetree-skyline-nelder-mead-optimizer.md): skyline uses Nelder-Mead instead of SLSQP
- [M-timetree-coalescent-missing-leaf-and-root-contributions.md](M-timetree-coalescent-missing-leaf-and-root-contributions.md): missing contributions in backward pass
