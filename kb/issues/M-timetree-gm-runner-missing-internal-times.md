# Golden master runner tests missing internal node times for 5 datasets

Golden master runner tests for dengue_20, lassa_L_20, mpox_clade_ii_20, rsv_a_20, and tb_20 fail across all three test modes (poisson, marginal dense, marginal sparse). The actual output contains only leaf times (matching raw date constraints), while the expected output from v0 contains both leaf and internal node times refined by inference.

## Affected tests

- `test_gm_runner_poisson` for dengue_20, lassa_L_20, mpox_clade_ii_20, rsv_a_20, tb_20
- `test_gm_runner_marginal_sparse` for same datasets
- `test_gm_runner_marginal_dense` for same datasets

Test files:

- [`packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_poisson.rs`](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_poisson.rs)
- [`packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs`](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs)
- [`packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs`](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs)

## Symptoms

For the Poisson dengue_20 case, expected output has 39 entries (20 leaves + 19 internal nodes) while actual output has only 20 entries (leaves only). Leaf times in the actual output are the raw date constraint values (mid-year .5 for integer year dates) rather than inference-refined values.

The failure is in the inference output itself. The tests reach the value-comparison stage and expose a semantic mismatch: internal node times are absent rather than merely named differently or hidden behind an earlier crash.

## Possible causes

- `extract_node_times()` returns `None` for `time()` on internal nodes, indicating backward/forward pass does not set time distributions on them for these datasets
- The golden master rerooted tree structure may differ from what v1's rerooting produces, leading to different internal node naming or topology
- These datasets have characteristics (short branches, specific branch length distributions) that cause the inference to not converge or produce empty distributions

## Working datasets

flu_h3n2_20 and ebola_20 pass in all test modes, suggesting the issue is data-dependent.

## Related tickets

- [kb/tickets/timetree-gm-runner-missing-internal-node-times.md](../tickets/timetree-gm-runner-missing-internal-node-times.md)
