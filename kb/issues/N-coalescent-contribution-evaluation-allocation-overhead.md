# Coalescent contribution evaluation allocates for scalar merger rates

Each internal-node coalescent contribution is represented by a separately allocated closure that clones shared interpolators through `Arc`. Every evaluation then allocates two one-element arrays to call an array-oriented merger-rate function for scalar inputs.

## Current overhead

[`compute_node_contributions()`](../../packages/treetime/src/coalescent/contributions.rs) creates one `DistributionFormula` closure per eligible node and stores the results in an `IndexMap<GraphNodeKey, Arc<DistributionNegLog>>`. Each internal closure owns cloned `Arc` handles for the same integral-merger-rate function, lineage-count function, and Tc distribution; only the node multiplicity differs.

Inside the closure, `k_t` and `tc_t` are scalars, but the code constructs `Array1::from_vec(vec![k_t])` and `Array1::from_vec(vec![tc_t])`, calls [`compute_merger_rates()`](../../packages/treetime/src/coalescent/integration.rs), and indexes the one-element total-rate result. This repeats heap allocation for every formula evaluation. Skyline likelihood evaluation in [`skyline.rs`](../../packages/treetime/src/coalescent/skyline.rs) separately duplicates the same scalar clamp and total-rate arithmetic.

## Ticket-ready scope

Extract the existing scalar merger-rate formula and make the array calculation delegate to it. Contributions and skyline evaluation can then use the scalar result directly, eliminating the one-element arrays and duplicated arithmetic without changing scientific behavior.

The broader single shared calendar-time coalescent model proposed in commit `542ac860c7cfa4bab6764aee1d1b3810a09eb54f` is not ready. It depends on unresolved time-coordinate sidedness at event breakpoints, skyline transformation semantics, multiplication ordering, and the missing leaf/root contribution design.

## Related ticket

- [kb/tickets/coalescent-extract-scalar-merger-rate-evaluation.md](../tickets/coalescent-extract-scalar-merger-rate-evaluation.md)

## Related issues

- [M-timetree-coalescent-missing-leaf-and-root-contributions.md](M-timetree-coalescent-missing-leaf-and-root-contributions.md)
- [M-timetree-coalescent-multiplication-ordering-diverges-from-v0.md](M-timetree-coalescent-multiplication-ordering-diverges-from-v0.md)
