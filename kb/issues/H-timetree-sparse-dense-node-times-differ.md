# Sparse and dense marginal timetree runs give different node times

The marginal time inference gives different node times in dense and in sparse mode for the same tree, alignment, clock rate, and indel policy. With JC69 and no per-site rates, the two modes compute the same branch-length likelihood and should give the same node times.

`test_gm_runner_marginal_sparse_dense_consistency` [packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs#L9-L25](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs#L9-L25) runs v1 dense and v1 sparse on the same inputs and compares their node times. On `flu_h3n2_20` it passes at `epsilon = 1e-3` (largest difference $7.3 \cdot 10^{-4}$ years at `NODE_0000015`); on the other datasets of the fixtures the modes differ by 0.07 to 1.1 years, so the test holds only the flu case.

On `rust` at `3825359a2`, one CLI comparison confirms the difference on the default path: `timetree` on `data/ebola/20` with `--clock-rate=0.001 --keep-root --model=jc69 --seed=1` gives node dates that differ by up to 0.0495 years (at `NODE_0000013`) between `--dense=true` and `--dense=false`. That run includes the ML branch-length pre-step, so it does not separate the pre-step from the time inference.

## Severity

Sparse is the default representation, because `fn infer_dense()` always returns `false` [packages/treetime/src/partition/algo/infer_dense.rs#L1-L7](../../packages/treetime/src/partition/algo/infer_dense.rs#L1-L7). Every default `timetree` run with an alignment therefore reports dates that differ from the dense result, which the project rules name as an independent oracle for sparse mode.

## Measurements

Largest node-time difference between the modes, in years, measured with the inputs of the `test_gm_runner` fixtures (the v0-rerooted tree and clock rate of `gm_runner_outputs.json`) on an unmerged stack of 2026-10-18 with an exact mass model and the v0 grid rule of the FFT convolution (commit `fda9ba5e6`). On the code of 2026-10-20 `flu_h3n2_20` differs by $7.3 \cdot 10^{-4}$ years and the other datasets by 0.07 to 1.1 years; the per-dataset values below are not remeasured:

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

Per-edge log-likelihoods of the two partitions, JC69, one `marginal_update` at the tree's branch lengths, on the 4-sequence, 40-site alignment `POLYTOMY_ALIGNMENT` of `packages/treetime/src/partition/__tests__/test_set.rs`: the whole-tree log-likelihoods agree within $10^{-12}$, but the per-edge log-likelihoods (`evaluate_with_indels_log_lh_only` at the branch lengths 0.001, 0.01 and 0.1) do not. Sparse minus dense:

| Tree | Edge | 0.001 | 0.01 | 0.1 |
| --- | --- | --- | --- | --- |
| `(((A:0.01,B:0.02)N:0.005,C:0.015)P:0.01,D:0.03)R;` | N to A | 0.03127 | 0.03117 | 0.03013 |
| same | N to B | 0.01092 | 0.01148 | 0.01118 |
| same | P to N | $1.30 \cdot 10^{-4}$ | $1.29 \cdot 10^{-4}$ | $1.23 \cdot 10^{-4}$ |
| same | P to C | 0.00681 | 0.00679 | 0.00658 |
| `((A:0.01,B:0.02,C:0.015)P:0.01,D:0.03)R;` | P to A | 0.03081 | 0.03071 | 0.02975 |
| same | P to B | 0.00963 | 0.01021 | 0.00995 |
| same | P to C | 0.00759 | 0.00838 | 0.00820 |

The edges from the root (R to P, R to D) agree within $10^{-14}$. The difference changes with the branch length by up to $10^{-3}$ nats, so it changes the shape of the branch-length likelihood, not only its level, and it appears only on edges whose parent is not the root. Suspect: the sparse message from the parent (the cavity), which the root edges do not use.

On `dengue_20`, all 38 branch-length distributions of the two modes have the same grid (300 points, same start and spacing), but different values: the peak moves by one grid cell on some edges (edge 0: 21.07 years dense, 20.68 years sparse), and the peak-normalized negative log-likelihood at the grid middle differs by about 3% (edge 0: 104.2 dense, 107.2 sparse). The difference therefore comes from the branch-length likelihood contributions of the two partitions, before the time passes.

The two datasets with the smallest and a large difference also differ in their alignments: `flu_h3n2_20` has 4 ambiguous characters and no gaps, `dengue_20` has 4 108 gaps and 295 ambiguous characters.

> [!IMPORTANT]
> **Investigation required.** Suspected cause: the dense and sparse partitions treat gaps or ambiguous characters differently in the branch-length likelihood. Compare the optimization contributions of one `dengue_20` edge in both modes as a function of branch length, and find the first site class where they differ.

## Related issues

The sparse backend is the default of every run ([N-representation-infer-dense-stub.md](N-representation-infer-dense-stub.md)), and the dense backend is its independent oracle. These issues show the two backends disagreeing on default settings and may share one cause. The smallest reproducing case is the tree with one unknown leaf of [H-sparse-marginal-likelihood-depends-on-branch-of-unknown-leaf.md](H-sparse-marginal-likelihood-depends-on-branch-of-unknown-leaf.md), so an investigation of the others starts there.

- [H-ancestral-sparse-dense-internal-residues-diverge.md](H-ancestral-sparse-dense-internal-residues-diverge.md)
- [H-sparse-marginal-likelihood-depends-on-branch-of-unknown-leaf.md](H-sparse-marginal-likelihood-depends-on-branch-of-unknown-leaf.md)
- [M-ancestral-sparse-root-invariance.md](M-ancestral-sparse-root-invariance.md)
