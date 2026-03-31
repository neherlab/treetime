# Grid search zero-comparison ignores indel likelihood

The grid search fallback in branch length optimization compares $t=0$ against the best positive grid point using substitution-only likelihood, while the grid itself is evaluated with the full indel-aware likelihood. When indels are present on an edge, this can declare $t=0$ optimal even though the Poisson indel log-likelihood is $-\infty$ at zero.

## Location

[packages/treetime/src/commands/optimize/optimize_unified.rs#L370-L374](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L370-L374)

```rust
let zero_is_better = contributions.iter().all(|c| c.all_sites_valid_at_zero()) && {
  let log_lh_zero = evaluate_mixed_log_lh_only(&contributions, 0.0);
  let log_lh_best = evaluate_mixed_log_lh_only(&contributions, best_positive);
  log_lh_zero > log_lh_best
};
```

The grid search at [L358-L364](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L358-L364) uses `evaluate_with_indels_log_lh_only`, which includes the Poisson indel contribution. The zero-vs-best comparison at L371-L372 uses `evaluate_mixed_log_lh_only`, which evaluates substitution likelihood only.

## Mathematical background

The Poisson indel model adds a log-likelihood contribution per edge:

$$\ell_{\text{indel}}(t) = k \ln(\mu t) - \mu t - \ln(k!)$$

where $k$ is the number of indel events and $\mu$ is the indel rate (see [docs/algorithms/optimize.md#L33-L36](../../docs/algorithms/optimize.md#L33-L36)).

At $t = 0$ with $k > 0$: $\ln(\mu \cdot 0) = -\infty$, so $\ell_{\text{indel}}(0) = -\infty$. The derivative $d\ell/dt = k/t - \mu$ diverges positively as $t \to 0^+$. Zero is never the optimum when indels are present.

The Newton path handles this correctly: `run_optimize_mixed()` at [L308-L314](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L308-L314) skips the zero-branch check when `indel_count > 0` and bumps the starting value away from zero. Only the grid search fallback has the inconsistency.

## Impact

Low in practice. Requires both conditions: (1) Newton fails on the edge (non-concave second derivative, falls through to grid) AND (2) indels are present on that edge. Both are individually uncommon in typical viral datasets. When triggered, assigns zero length to a branch that has indel evidence for finite length.

## v0 comparison

v0 uses Brent's method with fixed bracket [0, 4.0] for all branch length optimization. There is no grid search fallback. Indels are not modeled in v0 branch length optimization at all. This bug is v1-only.

## Cross-links

- [H-optimize-sparse-hessian-multiplicity.md](H-optimize-sparse-hessian-multiplicity.md) -- grid fallback triggers when Newton fails; if the Hessian bug is fixed, fewer edges fall through to grid
- [M-optimize-grid-search-narrow-at-zero.md](M-optimize-grid-search-narrow-at-zero.md) -- related grid search issue affecting range coverage
- [Indel contribution intentional change](../port-intentional-changes/optimize-indel-contribution-to-likelihood.md) -- Poisson model documentation
- [Indel models algorithm inventory](../port-algo-inventory/indel-models.md) -- full catalog of indel modeling approaches
- [Indel model alternatives proposal](../port-proposals/optimize-indel-model-alternatives.md) -- future extensions

## Fix

Use `evaluate_with_indels_log_lh_only` for the zero-vs-best comparison, consistent with the grid evaluation. When `indel_count > 0`, skip the zero candidate entirely (matching the Newton path logic).
