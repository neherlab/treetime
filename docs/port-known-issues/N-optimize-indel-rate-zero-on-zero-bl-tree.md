# Indel rate estimation returns zero on all-zero-BL input trees

When an input tree has all branch lengths set to zero, `estimate_indel_rate()` returns 0 because the denominator (total branch length) is zero. With `InitialGuessMode::Auto` (default), `initial_guess_mixed()` preserves these zero values, so the denominator stays zero through the first optimization pass. The Poisson indel term is inactive for that pass because `poisson_indel_log_lh(k, 0.0, t)` returns neutral metrics.

## Location

- [packages/treetime/src/commands/optimize/optimize_indel.rs#L20-L23](../../packages/treetime/src/commands/optimize/optimize_indel.rs#L20-L23): early return when `mu == 0.0`
- [packages/treetime/src/commands/optimize/optimize_unified.rs#L289](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L289): `estimate_indel_rate` called once per `run_optimize_mixed` invocation
- [packages/treetime/src/commands/optimize/run.rs#L132-L133](../../packages/treetime/src/commands/optimize/run.rs#L132-L133): `Auto` mode passes `overwrite_valid=false`

## Self-correction

The multi-iteration loop in `run.rs` compensates:

1. First pass: `estimate_indel_rate()` returns 0, but `run_optimize_mixed` bumps indel-bearing edges from zero to `one_mutation` (the `branch_length == 0.0 && indel_count > 0` guard).
2. After first pass: edges have positive BLs, so `estimate_indel_rate()` returns a positive value.
3. Second pass: the Poisson term is active, branches converge toward the indel optimum.

The cost is one extra iteration before the indel objective is active, not incorrect final results.

## Impact

Negligible. Requires all-zero input BLs with identical sequences and indel-only evidence. In practice, input trees from phylogenetic inference tools (RAxML, IQ-TREE, FastTree) have positive branch lengths. The convergence note in the [indel contribution intentional change](../port-intentional-changes/optimize-indel-contribution-to-likelihood.md) documents this limitation.

## v0 comparison

v0 does not model indels in branch length optimization. Not applicable.
