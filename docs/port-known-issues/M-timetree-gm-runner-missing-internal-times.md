# Golden master runner tests missing internal node times for 5 datasets

Golden master runner tests for dengue_20, lassa_L_20, mpox_clade_ii_20, rsv_a_20, and tb_20 fail across all three test modes (poisson, marginal dense, marginal sparse). The actual output contains only leaf times (matching raw date constraints), while the expected output from v0 contains both leaf and internal node times refined by inference.

## Affected tests

- `test_gm_runner_poisson` for dengue_20, lassa_L_20, mpox_clade_ii_20, rsv_a_20, tb_20
- `test_gm_runner_marginal_sparse` for same datasets
- `test_gm_runner_marginal_dense` for same datasets

Test files:

- [`packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_poisson.rs`](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_poisson.rs)
- [`packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs`](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs)
- [`packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs`](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs)

## Symptoms

For the Poisson dengue_20 case, expected output has 39 entries (20 leaves + 19 internal nodes) while actual output has only 20 entries (leaves only). Leaf times in the actual output are the raw date constraint values (mid-year .5 for integer year dates) rather than inference-refined values.

## Possible causes

- `extract_node_times()` returns `None` for `time()` on internal nodes, indicating backward/forward pass does not set time distributions on them for these datasets
- The golden master rerooted tree structure may differ from what v1's rerooting produces, leading to different internal node naming or topology
- These datasets have characteristics (short branches, specific branch length distributions) that cause the inference to not converge or produce empty distributions

## History

These tests were previously disabled with TODO comments citing the grid crash (`H-timetree-crash-grid-zero-branch`). The grid crash was fixed by commit `5e4ae89a` (avoid grid x-array roundtrip in distribution operations). The tests now run without crashing but fail with value mismatches.

## Working datasets

flu_h3n2_20 and ebola_20 pass in all test modes, suggesting the issue is data-dependent.
