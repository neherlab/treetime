# Marginal dense golden master node key mismatch on ebola_20

The marginal dense runner test golden master for ebola_20 has a different set of internal node keys than v1 produces. The v0-captured golden data contains 11 internal nodes while v1's rerooting produces 19. The `pretty_assert_map_abs_diff_eq!` macro requires exact key match, so the test fails on key comparison before values are checked.

The dense marginal pipeline runs correctly on ebola_20. The mismatch is in the golden master data captured from v0, not in the inference results.

## Location

Commented-out test case in [test_gm_runner_marginal_dense.rs#L32](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L32).

Golden master data in [gm_runner_outputs.json](../../packages/treetime/src/timetree/inference/__tests__/__fixtures__/gm_runner_outputs.json).

## Resolution options

- Recapture golden master data from v0 using the same rerooted tree that v1 produces
- Compare only the intersection of node keys (leaf nodes match, internal node naming diverges)
