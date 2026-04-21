# Optimize loop likelihood excludes indel term

The branch-length optimizer maximizes `substitution_lh + poisson_indel_lh` per edge, but the outer optimize loop at [packages/treetime/src/commands/optimize/run.rs#L280-L283](../../packages/treetime/src/commands/optimize/run.rs#L280-L283) totals only the substitution likelihood from `update_marginal`. The indel contribution is computed and used inside per-edge optimization but is not accumulated into `total_lh`.

## Impact

Every stopping condition (`Converged`, `Oscillating`, `Worsened`), the `best_lh` tracking, and the worsened-state revert at [packages/treetime/src/commands/optimize/run.rs#L307-L336](../../packages/treetime/src/commands/optimize/run.rs#L307-L336) test a scalar that does not reflect the full objective. The `Worsened` branch reverts to a previous state based on substitution-only likelihood, which discards improvements in the joint objective.

On datasets with indels (sc2, ebola, mpox), the reported likelihood drifts non-monotonically because the indel term is invisible to the convergence monitor.

## Affected code

- Likelihood accumulation: [packages/treetime/src/commands/optimize/run.rs#L280-L283](../../packages/treetime/src/commands/optimize/run.rs#L280-L283)
- Convergence and revert logic: [packages/treetime/src/commands/optimize/run.rs#L307-L336](../../packages/treetime/src/commands/optimize/run.rs#L307-L336)
- Per-edge optimization includes indels: [packages/treetime/src/commands/optimize/optimize_unified.rs#L259](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L259) (`evaluate_with_indels_log_lh_only`)

## Fix

After each `run_optimize_mixed` call, accumulate `poisson_indel_log_lh` for all edges and add it to `total_lh`. The indel rate per edge is already computed during optimization and can be cached or recomputed.

## v0 comparison

v0 (`TreeTime.py`) has a similar structure where the reported likelihood does not include the indel term, so this is a shared defect rather than a v0/v1 divergence. v1 should fix it regardless.
