# Sparse and dense marginal timetree runs give different node times

The marginal time inference gives different node times in dense and in sparse mode for the same tree, alignment, clock rate, and indel policy. With JC69 and no per-site rates, the two modes compute the same branch-length likelihood and should give the same node times.

The measurements below come from a version of `test_gm_runner_marginal_sparse_dense_consistency` that runs v1 dense and v1 sparse on the same inputs and compares their node times; on `flu_h3n2_20` it passed only at `epsilon = 1e-3` (years).

> [!WARNING]
> **Needs review.** That version of the test and the measurements exist only on a stack of 2026-10-18 on top of `rust` of that day, with an exact mass model, an exact convolution of ranges, and the v0 grid rule of the FFT convolution (commit `fda9ba5e6`, local branch `worktree/fix/distribution-review-findings`). On `rust` of 2026-10-19, `test_gm_runner_marginal_sparse_dense_consistency` [packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs#L45-L58](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs#L45-L58) still compares v1 sparse with the v0 dense fixture, which the project rules forbid as an oracle for sparse mode, and is ignored with the reason "max 0.92 years at root-adjacent node". The table has not been measured on `rust`.

On `rust` at `3825359a2`, one CLI comparison confirms the difference on the default path: `timetree` on `data/ebola/20` with `--clock-rate=0.001 --keep-root --model=jc69 --seed=1` gives node dates that differ by up to 0.0495 years (at `NODE_0000013`) between `--dense=true` and `--dense=false`. That run includes the ML branch-length pre-step, so it does not separate the pre-step from the time inference.

## Severity

Sparse is the default representation, because `fn infer_dense()` always returns `false` [packages/treetime/src/partition/algo/infer_dense.rs#L1-L7](../../packages/treetime/src/partition/algo/infer_dense.rs#L1-L7). Every default `timetree` run with an alignment therefore reports dates that differ from the dense result, which the project rules name as an independent oracle for sparse mode.

## Measurements

Largest node-time difference between the modes, in years, measured with the inputs of the `test_gm_runner` fixtures (the v0-rerooted tree and clock rate of `gm_runner_outputs.json`):

| Dataset | Include indels | Exclude indels | Node |
| --- | --- | --- | --- |
| `flu_h3n2_20` | 4.9e-4 | 4.9e-4 | `NODE_0000013` |
| `ebola_20` | 0.065 | 0.098 | `NODE_0000014` |
| `tb_20` | 0.229 | 0.229 | `NODE_0000007` |
| `dengue_20` | 0.438 | 0.462 | `NODE_0000001` |
| `rsv_a_20` | 0.697 | 0.722 | `NODE_0000009` |
| `lassa_L_20` | 1.321 | 1.325 | `NODE_0000005` |

- `mpox_clade_ii_20` exceeds the 300 s test timeout for the dense and the sparse run together
- `zika_20` fails while reading the dates, before inference (see [M-timetree-date-header-hash.md](M-timetree-date-header-hash.md))

## Evidence

On `dengue_20`, all 38 branch-length distributions of the two modes have the same grid (300 points, same start and spacing), but different values: the peak moves by one grid cell on some edges (edge 0: 21.07 years dense, 20.68 years sparse), and the peak-normalized negative log-likelihood at the grid middle differs by about 3% (edge 0: 104.2 dense, 107.2 sparse). The difference therefore comes from the branch-length likelihood contributions of the two partitions, before the time passes.

The two datasets with the smallest and a large difference also differ in their alignments: `flu_h3n2_20` has 4 ambiguous characters and no gaps, `dengue_20` has 4 108 gaps and 295 ambiguous characters.

> [!IMPORTANT]
> **Investigation required.** Suspected cause: the dense and sparse partitions treat gaps or ambiguous characters differently in the branch-length likelihood. Compare the optimization contributions of one `dengue_20` edge in both modes as a function of branch length, and find the first site class where they differ.
