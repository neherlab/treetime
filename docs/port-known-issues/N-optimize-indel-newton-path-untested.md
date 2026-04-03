# Indel-aware Newton optimization path has no direct tests

The Newton optimization loop in `run_optimize_mixed()` includes indel-aware evaluation (`evaluate_with_indels()`), a `min_branch_length` clamp for indel-bearing edges, and an `indel_count == 0` guard before the zero-branch check. None of these are tested directly. Existing integration tests assert only `bl > 0.0` and `bl.is_finite()`.

## Location

- [packages/treetime/src/commands/optimize/optimize_unified.rs#L375-L440](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L375-L440): Newton loop with indel contributions
- [packages/treetime/src/commands/optimize/optimize_unified.rs#L398](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L398): `min_branch_length` clamping for indel edges
- [packages/treetime/src/commands/optimize/optimize_unified.rs#L388](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L388): `indel_count == 0` guard

## Impact

Medium. The Poisson indel model changes Newton step behavior (additional log-likelihood term, derivative divergence at zero, minimum branch length clamping). Without direct tests, regressions in the indel-aware path are detected only through end-to-end integration tests, which lack precision.

## Tests needed

- Unit test: Newton with `indel_count > 0` converges to the Poisson MLE `k/mu` for indel-only edges (zero substitution signal)
- Unit test: `min_branch_length` clamping prevents Newton from reaching zero on indel-bearing edges
- Unit test: `evaluate_with_indels()` returns correct log-likelihood for known Poisson parameters
- Property test: `evaluate_with_indels(contributions, k, mu, t) >= evaluate_mixed(contributions, t)` for all `k >= 0` (indel term is non-positive for k=0, so this holds when indels are present)
