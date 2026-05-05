# normalize_inplace NEG_INFINITY masks all contributions

## Summary

`normalize_inplace` returns `f64::NEG_INFINITY` for degenerate rows and substitutes uniform distributions. A single degenerate row makes the total log-likelihood negative infinity, masking contributions from all well-determined rows.

## Details

`packages/treetime/src/representation/partition/marginal_dense.rs:390:`

The function iterates rows, normalizing each to sum to 1.0. For rows where the sum is zero or non-finite, it fills the row with a uniform distribution (1/n_cols) and adds `NEG_INFINITY` to the running `log_lh` total.

Since `log_lh` accumulates via addition, a single degenerate row produces:

```
log_lh = ... + finite_contributions + NEG_INFINITY = NEG_INFINITY
```

Convergence checks downstream see `NEG_INFINITY` even when the rest of the tree is well-determined. The degenerate row's contribution dominates and makes it impossible to detect likelihood improvements.

## Impact

- Convergence criteria based on likelihood changes become unreliable
- A single degenerate position (e.g., all-gap column) masks the entire tree's likelihood signal
- Optimizer cannot distinguish "improving" from "stuck" when total is negative infinity

## Fix

Options:

1. Return a finite penalty for degenerate rows (e.g., `ln(1/n_cols) * n_positions`) instead of `NEG_INFINITY`
2. Track degenerate row count separately and exclude from convergence arithmetic
3. Use a masked accumulator that skips degenerate contributions
