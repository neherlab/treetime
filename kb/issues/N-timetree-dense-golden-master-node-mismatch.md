# Marginal dense golden master node key mismatch on ebola_20

The ebola_20 `marginal_dense` golden master has a different set of internal node keys than v1 produces. It holds 11 internal nodes (`NODE_0000000`, `1`, `2`, `5`, `8`, `9`, `11`, `13`, `16`, `17`, `18`), while the rerooted tree of the same fixture holds all 19 (the `poisson` golden master that also held them was removed in 88b7caede) (`NODE_0000000` to `NODE_0000018`). The `pretty_assert_map_abs_diff_eq!` macro requires an exact key match, so the test fails on key comparison before values are checked. The mismatch is in the golden master data, not in the inference results.

## Cause

The capture script [packages/treetime/src/timetree/inference/__tests__/__fixtures__/gm_runner_capture](../../packages/treetime/src/timetree/inference/__tests__/__fixtures__/gm_runner_capture) runs v0 on a copy of the same rerooted tree. The `marginal_dense` run passes the alignment with `branch_length_mode="marginal"`, so v0 `TreeTime.run()` calls `optimize_tree()` with the default `prune_short=True` [packages/legacy/treetime/treetime/treetime.py#L215](../../packages/legacy/treetime/treetime/treetime.py#L215), [#L243](../../packages/legacy/treetime/treetime/treetime.py#L243). `TreeAnc.prune_short_branches()` then removes every internal node whose branch is shorter than 0.1 mutations and whose sequence matches its parent's with probability above 0.1 [packages/legacy/treetime/treetime/treeanc.py#L1475-L1495](../../packages/legacy/treetime/treetime/treeanc.py#L1475-L1495). The 8 missing internal nodes are the nodes this step removed.

## Location

- Test files: the ebola_20 case is commented out under a `// TODO:` line in [packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L35-L38](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L35-L38) and in [packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs#L35-L38](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs#L35-L38). The dense test is also `#[ignore]` for a separate dense-versus-v0 discrepancy on `flu_h3n2_20` ([M-timetree-branch-grid-uniform-resolution.md](M-timetree-branch-grid-uniform-resolution.md))
- Golden master data: [packages/treetime/src/timetree/inference/__tests__/__fixtures__/gm_runner_outputs.json](../../packages/treetime/src/timetree/inference/__tests__/__fixtures__/gm_runner_outputs.json)

## Resolution options

> [!IMPORTANT]
> **Decision required.** How to make the ebola_20 marginal golden master comparable with v1. Options: (a) recapture the `marginal_dense` data from v0 with `prune_short=False`, so v0 keeps the 19 internal nodes of the rerooted tree; (b) compare only the intersection of node keys, which leaves the 8 pruned nodes untested. Evidence: the rerooted input tree of the fixture holds all 19 internal nodes.
