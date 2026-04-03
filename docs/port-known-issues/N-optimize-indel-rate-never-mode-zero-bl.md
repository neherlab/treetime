# Indel rate estimation bypassed in Never mode on all-zero-BL trees

When `--branch-length-initial-guess=never` is used with a tree where all branch lengths are exactly zero and indels are present, `estimate_indel_rate()` returns 0.0 because total branch length is zero. The Poisson indel log-likelihood contributes nothing throughout the optimization.

## Location

- [packages/treetime/src/commands/optimize/run.rs#L138-L146](../../packages/treetime/src/commands/optimize/run.rs#L138-L146): `Never` mode only rejects `None`/`NaN`, not zero with indels
- [packages/treetime/src/commands/optimize/optimize_unified.rs#L358](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L358): `estimate_indel_rate()` called once before per-edge loop, returns 0 when total BL is zero
- [packages/treetime/src/commands/optimize/optimize_indel.rs#L20-L23](../../packages/treetime/src/commands/optimize/optimize_indel.rs#L20-L23): `poisson_indel_log_lh()` early return when mu == 0.0

## Mechanism

The `Auto` path treats zero-BL indel-bearing edges as invalid in `initial_guess_mixed()`, seeding positive BLs before `run_optimize_mixed()`. The `Never` path skips `initial_guess_mixed()` entirely, so `estimate_indel_rate()` sees total BL = 0 and returns zero. The per-edge bootstrap at `optimize_unified.rs:L377-L382` sets a local `branch_length = one_mutation`, but `indel_rate` was already computed as zero.

## Impact

Negligible. Requires all five conditions simultaneously:

1. All input branch lengths exactly zero
2. Sequences identical (no substitution signal)
3. Indels present on at least one edge
4. `InitialGuessMode::Never` explicitly chosen
5. The convergence check fires before `estimate_indel_rate` gets positive total BL on second iteration

Users choosing `--branch-length-initial-guess=never` provide trees with pre-computed positive branch lengths in practice.

## Fix options

- Recompute `indel_rate` after the per-edge zero-BL bootstrap in `run_optimize_mixed()`, or move the bootstrap before `estimate_indel_rate()`
- Validate in `Never` mode: reject zero-BL edges that carry indels (matching the `Auto` path logic)

## v0 comparison

v0 does not model indels in branch length optimization. Not applicable.
