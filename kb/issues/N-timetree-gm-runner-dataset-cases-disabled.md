# Golden-master runner tests cover only `flu_h3n2_20`

The timetree runner tests run the marginal time inference on the inputs of the `gm_runner` fixtures (the v0-rerooted tree and clock rate of `gm_runner_outputs.json`). Both run only the case `flu_h3n2_20`:

- `test_gm_runner_marginal_dense` compares v1 dense node times with v0 [packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L9-L33](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L9-L33). It passes at `epsilon = 1e-1` (largest gap 0.070 years, [M-timetree-branch-grid-uniform-resolution.md](M-timetree-branch-grid-uniform-resolution.md)); the cases of the other seven datasets are commented out under a `// TODO:` line that cites this issue
- `test_gm_runner_marginal_sparse_dense_consistency` compares v1 sparse with v1 dense [packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs#L9-L25](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs#L9-L25). It passes at `epsilon = 1e-3` (largest gap $7.3 \cdot 10^{-4}$ years, [H-timetree-sparse-dense-node-times-differ.md](H-timetree-sparse-dense-node-times-differ.md)); the other datasets differ by 0.07 to 1.1 years and have no case

The fixtures still hold the expected values of all datasets.

## Measured state per dataset (dense against v0, include and exclude indels)

| Dataset | Outcome |
| --- | --- |
| `dengue_20` | 14.8 and 16.3 years at `NODE_0000018` |
| `lassa_L_20` | 4.9 and 6.7 years at `NODE_0000004` |
| `rsv_a_20` | 1.8 and 1.7 years at `NODE_0000017` |
| `ebola_20` | v1 dates 8 internal nodes that v0 leaves out ([N-timetree-dense-golden-master-node-mismatch.md](N-timetree-dense-golden-master-node-mismatch.md)) |
| `tb_20` | v1 dates 6 internal nodes that v0 leaves out; 33.9 years at `NODE_0000016` |
| `mpox_clade_ii_20` | exceeds the 300 s test timeout |
| `zika_20` | v1 dates 3 internal nodes that v0 leaves out ([N-timetree-dense-golden-master-node-mismatch.md](N-timetree-dense-golden-master-node-mismatch.md)); 0.077 years at `NODE_0000012` |

The gaps on dengue, lassa and rsv are far above the 0.070 years of flu. No measurement separates a v1 defect from a fixture effect.

> [!IMPORTANT]
> **Investigation required.** The v0 capture script `packages/treetime/src/timetree/inference/__tests__/__fixtures__/gm_runner_capture` does not reproduce the committed `marginal_dense` values on the development host: a rerun differs from the fixture by up to 0.30 years on `lassa_L_20`. Find out whether the fixture or the capture environment changed before reading the gaps above as v1 errors.

> [!IMPORTANT]
> **Investigation required.** Why the dense `mpox_clade_ii_20` case is slow. The dense timetree on `data/mpox/clade-ii/20` takes 66.8 s of wall time (410 s of CPU on all threads) against 0.62 s for `data/flu/h3n2/20`. Profile the dense run and find the operation that dominates.

## Fix

- Find the cause of each dataset's gap, fix v1 or the fixture, and enable the case at the tightest passing tolerance, citing the issue that holds the remaining gap
- File a separate issue for each new cause

## Validation

The runner tests hold a case for every dataset of the fixtures, each at a tolerance of at most `1e-6` or citing the issue of its measured gap.
