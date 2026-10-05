# Marginal dense golden master node key mismatch on ebola_20

> [!WARNING]
> **Needs review.** The ebola_20 marginal dense runner case is gone, not commented out. Commit `3aca7de4` (refactor(treetime): remove comments) deleted the commented-out `#[case::ebola_20("ebola_20")]` lines from the marginal dense and marginal sparse test files; [test_gm_runner_marginal_dense.rs#L32](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L32) now holds the `flu_h3n2_20` case. The ebola_20 expected values remain in the fixture. Coverage for ebola_20 in both modes must be re-added, and the 11-versus-19 key count must be re-measured when it is.

The marginal dense runner test golden master for ebola_20 has a different set of internal node keys than v1 produces. The v0-captured golden data contains 11 internal nodes while v1's rerooting produces 19. The `pretty_assert_map_abs_diff_eq!` macro requires exact key match, so the test fails on key comparison before values are checked.

The dense marginal pipeline runs correctly on ebola_20. The mismatch is in the golden master data captured from v0, not in the inference results.

## Location

Test file: [test_gm_runner_marginal_dense.rs](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs) (the ebola_20 case is absent).

Golden master data in [gm_runner_outputs.json](../../packages/treetime/src/timetree/inference/__tests__/__fixtures__/gm_runner_outputs.json).

## Resolution options

- Recapture golden master data from v0 using the same rerooted tree that v1 produces
- Compare only the intersection of node keys (leaf nodes match, internal node naming diverges)
