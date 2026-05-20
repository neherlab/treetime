# Gate initial_guess_mixed() with --no-indels flag

`initial_guess_mixed()` at `packages/treetime/src/optimize/dispatch.rs:748` queries `edge_indel_count()` and uses the indel rate to compute initial branch length guesses in Auto and Always modes. When `--no-indels` is set, the optimize loop and per-edge optimizer correctly zero indel contributions, but the pre-loop initialization step does not receive the `no_indels` flag.

This means the initial branch length guess can still be influenced by indel counts even when `--no-indels` is set, causing the optimizer to start from a different point than expected.

## Fix

Pass `no_indels` into `initial_guess_mixed()`. When true, zero `indel_count` and `indel_rate` in the initial guess computation. Add test for Auto/Always modes with `no_indels=true` verifying indel-free initialization.

## Related issues

- Source: [M-optimize-no-indels-initial-guess-not-gated.md](../issues/M-optimize-no-indels-initial-guess-not-gated.md) -- delete after full resolution

Two additional cleanup items in the same area:

- `run_optimize_mixed_with_indel_rate()` at `packages/treetime/src/optimize/dispatch.rs:552` is unused in production after the refactor (only called by tests). Gate with `#[cfg(test)]` to eliminate the dead-code warning.
- `test_no_indels.rs` has private copies of `TREE_ZERO_BL`, `setup_dense_with_marginal`, `inject_indel_on_first_edge` that exist as `pub` in `test_initial_guess_mode.rs` helpers module. Deduplicate by importing from the shared module.
