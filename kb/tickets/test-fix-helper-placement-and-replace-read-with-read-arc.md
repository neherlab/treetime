# Fix helper placement violations and replace .read() with .read_arc()

## Helper placement violations (9 files)

Helper functions before tests and/or not wrapped in `mod helpers`:

1. `seq/__tests__/test_mutation.rs`
2. `seq/__tests__/test_find_char_ranges.rs`
3. `alphabet/__tests__/test_alphabet_config.rs`
4. `commands/timetree/coalescent/__tests__/test_skyline.rs`
5. `commands/timetree/optimization/__tests__/test_relaxed_clock.rs`
6. `commands/timetree/optimization/__tests__/test_polytomy.rs`
7. `commands/timetree/convergence/__tests__/test_metrics.rs`
8. `commands/timetree/inference/__tests__/test_branch_length_likelihood.rs`
9. `gtr/__tests__/test_gm_gtr_site_specific.rs:170-207:`

## .read() instead of .read_arc() (8 instances)

`packages/treetime/src/commands/ancestral/__tests__/test_python_parity.rs:192,240,283,327,390,391,463,473:`

Test code calls `.read()` on `parking_lot::RwLock` values wrapped in `Arc`. Per project convention, `.read_arc()` should be used to return an `ArcRwLockReadGuard` that keeps the `Arc` alive.

## Related issues

- Source: [N-test-quality-deficiencies.md](../issues/N-test-quality-deficiencies.md) -- delete after full resolution
