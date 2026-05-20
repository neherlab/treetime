# initial_guess_mixed allocates Vec<Sub> per edge for count only

`initial_guess_mixed()` ([packages/treetime/src/optimize/dispatch.rs#L309](../../packages/treetime/src/optimize/dispatch.rs#L309)) calls `edge_subs()` and takes `.len()`, discarding the substitution vector. For dense partitions, `edge_subs()` ([packages/treetime/src/partition/marginal_dense.rs#L87-L118](../../packages/treetime/src/partition/marginal_dense.rs#L87-L118)) builds a `Vec<Sub>` with one allocation per differing site. The initializer runs across every edge before optimization.

A count-only method (incrementing a counter instead of pushing `Sub` values) would avoid the allocation. For small alignments the cost is negligible. For large alignments (thousands of positions, hundreds of edges) the cumulative allocation is avoidable.

## Related issues

- Source: [N-optimize-initial-guess-alloc.md](../issues/N-optimize-initial-guess-alloc.md) -- delete after full resolution
