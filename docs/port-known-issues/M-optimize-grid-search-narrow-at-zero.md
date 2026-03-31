# Grid search fallback covers narrow range at zero branch length

The grid search fallback at [packages/treetime/src/commands/optimize/optimize_unified.rs#L275](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L275) uses a range proportional to the current branch length:

```rust
let branch_lengths = ndarray::Array1::linspace(0.1 * one_mutation, 1.5 * branch_length + one_mutation, 100);
```

When `branch_length = 0.0`, this becomes `linspace(0.1/L, 1.0/L, 100)` where L is total alignment length. For a 1000-site alignment, the grid covers only `[1e-4, 1e-3]`, missing branches longer than `1/L` substitutions per site.

## Impact

Medium. Branch lengths in phylogenetic trees span 3-4 orders of magnitude ($10^{-5}$ to $10^{-1}$ subs/site). A grid covering only one decade near `1/L` misses the long-branch regime. The outer coordinate-ascent loop partially compensates by re-optimizing each branch in subsequent iterations, but convergence is slower.

## Proposed fix

Use a minimum upper bound: `linspace(0.1 * one_mutation, max(1.5 * branch_length + one_mutation, 0.5), 100)`, or use a log-spaced grid to cover the full biologically plausible range.
