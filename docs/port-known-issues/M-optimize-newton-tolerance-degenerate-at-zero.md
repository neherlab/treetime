# Newton convergence tolerance degenerates at zero branch length

The Newton inner-loop convergence check at [packages/treetime/src/commands/optimize/optimize_unified.rs#L258](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L258) uses a purely relative tolerance:

```rust
while (new_branch_length - branch_length).abs() > 0.001 * branch_length && n_iter < max_iter {
```

When `branch_length = 0.0`, the tolerance becomes `0.0`, so the condition `|delta| > 0.0` is always true for any nonzero Newton step. The inner loop exhausts all 10 iterations.

## Impact

Medium. After the first inner iteration, `branch_length` is updated to a nonzero value, so the tolerance becomes meaningful for remaining iterations. Damage limited to one wasted iteration per zero-starting-branch edge. No incorrect results since the iteration cap prevents infinite loops.

## Proposed fix

Add an absolute tolerance floor: `0.001 * branch_length + 1e-10` or `max(0.001 * branch_length, one_mutation * 1e-3)`. Standard approach in phylogenetic optimizers (RAxML, IQ-TREE).
