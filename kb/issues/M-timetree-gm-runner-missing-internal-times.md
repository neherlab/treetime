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

## Root cause

The backward pass multiplies child parent-time messages at each internal node. Two issues caused the product to collapse to `Empty`, leaving internal time distributions unset:

1. `distribution_multiplication` ignored operand `Constant` tails, producing `Empty` when child messages had disjoint finite grids. Fixed: multiplication now extends evaluable domain on sides with `Constant` tails.
2. `normalize()` between successive child multiplications reset tails to `Error`, so the accumulated result lost its `Constant` left tail before the next multiplication. Fixed: backward pass re-applies `Constant` left / `Zero` right after each normalize.

Both fixes are in the `fix/multiply-honor-tails` branch. End-to-end verification on flu/h3n2/20 (36/36 nodes dated) and mpox/hmpxv1 (1451/1451 nodes dated) confirms internal times are now produced under `--keep-root`.

## Remaining verification

The 15 commented-out GM runner test cases (5 datasets * 3 modes) need to be uncommented and run. They may still fail for independent reasons (grid-width tolerance, rerooting topology differences) even though the missing-times root cause is addressed.

## Working datasets

flu_h3n2_20 and ebola_20 pass in all test modes, suggesting the issue is data-dependent.

## Related tickets

- [kb/tickets/timetree-gm-runner-missing-internal-node-times.md](../tickets/timetree-gm-runner-missing-internal-node-times.md)
