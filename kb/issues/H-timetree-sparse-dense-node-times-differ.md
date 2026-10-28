# Sparse and dense marginal timetree runs give different node times

The marginal time inference gives different node times in dense and in sparse mode for the same tree, alignment, clock rate, and indel policy. With JC69 and no per-site rates, the two modes compute the same branch-length likelihood and should give the same node times.

> [!WARNING]
> **Needs review.** The sparse backend now leaves sites where a child subtree has no data out of its messages and its branch-length objective, and keeps the data below ancestors that Fitch resolves as gaps (commits `e380322ea`, `a084a38c8`, `ef5cc6a21`). Before this fix, sparse branch lengths were shorter than dense ones on alignments with gaps and `N`. The test `test_missing_data_sparse_matches_dense` in `packages/treetime/src/partition/marginal/sparse/__tests__/test_missing_data.rs` compares the root log-likelihood and the edge objectives of both backends on such alignments. Measure the node times of both modes again on the datasets above, and rewrite or delete this issue

`test_gm_runner_marginal_sparse_dense_consistency` [packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs#L9-L25](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs#L9-L25) runs v1 dense and v1 sparse on the same inputs and compares their node times. On `flu_h3n2_20` it passes at `epsilon = 1e-2` (largest difference $4.1 \cdot 10^{-3}$ years at `NODE_0000011`); on the other datasets of the fixtures the modes differ by 0.07 to 1.2 years, so the test holds only the flu case. The node time is the grid point with the largest posterior density, and the two modes build different posterior grids, so on `flu_h3n2_20` the test measures mostly the grid cell of about $6.4 \cdot 10^{-3}$ years (see Evidence).

On `rust` at `3825359a2`, one CLI comparison confirms the difference on the default path: `timetree` on `data/ebola/20` with `--clock-rate=0.001 --keep-root --model=jc69 --seed=1` gives node dates that differ by up to 0.0495 years (at `NODE_0000013`) between `--representation dense` and `--representation sparse`. That run includes the ML branch-length pre-step, so it does not separate the pre-step from the time inference.

## Severity

Sparse is the default representation, because `fn infer_dense()` always returns `false` [packages/treetime/src/partition/algo/infer_dense.rs#L1-L7](../../packages/treetime/src/partition/algo/infer_dense.rs#L1-L7). Every default `timetree` run with an alignment therefore reports dates that differ from the dense result, which the project rules name as an independent oracle for sparse mode.

## Measurements

Largest node-time difference between the modes, in years, measured with the inputs of the `test_gm_runner` fixtures (the v0-rerooted tree and clock rate of `gm_runner_outputs.json`), JC69, no coalescent prior. "Node time" is the time the test compares (the grid argmax of the posterior, not earlier than the parent); "median" is the largest difference of the posterior medians, which does not depend on the grid. Code: `df43a205d` (trapezoid rule in the convolutions); in parentheses the same measurement on `rust` at `f4bd506d0` (rectangle sums), where it differs:

| Dataset | Node time, include indels | Node time, exclude indels | Node | Median, include indels | Median, exclude indels |
| --- | --- | --- | --- | --- | --- |
| `flu_h3n2_20` | 4.1e-3 (7.3e-4) | 4.1e-3 (7.3e-4) | `NODE_0000011` (`NODE_0000015`) | 7.4e-4 (4.2e-4) | 7.4e-4 (4.2e-4) |
| `ebola_20` | 0.070 (0.071) | 0.099 (0.096) | `NODE_0000014` | 0.069 | 0.097 (0.098) |
| `tb_20` | 0.237 (0.238) | 0.237 (0.238) | `NODE_0000007` | 0.235 (0.262) | 0.235 (0.262) |
| `dengue_20` | 0.605 (0.576) | 0.626 (0.627) | `NODE_0000001` | 0.575 | 0.582 |
| `rsv_a_20` | 0.681 (0.692) | 0.721 (0.713) | `NODE_0000009` | 0.688 | 0.698 |
| `lassa_L_20` | 0.975 (0.931) | 1.167 (1.106) | `NODE_0000005` | 1.046 (0.971) | 0.995 (0.964) |

- `mpox_clade_ii_20` exceeds the 300 s test timeout for the dense and the sparse run together
- `zika_20` differs by $1.8 \cdot 10^{-3}$ years at `NODE_0000014` with indels included and excluded (code of 2026-10-23)

## Evidence

The time passes read the two modes only through their branch-length distributions: on every dataset of the table, both indel policies, and both code states above, the sparse setup run with the dense branch-length distributions gives the dense node times exactly, and the dense setup run with the sparse distributions gives the sparse node times exactly (difference 0). The whole difference therefore comes from the branch-length distributions of the two partitions.

The node time is the grid point with the largest posterior density (`fn likely_time()` [packages/treetime-distribution/src/distribution_core/function.rs#L297-L314](../../packages/treetime-distribution/src/distribution_core/function.rs#L297-L314)). The two modes have slightly different branch-length optima, so their posterior grids differ in start and spacing. On `flu_h3n2_20` at `df43a205d`, the posterior of `NODE_0000011` has 631 points with spacing $6.350 \cdot 10^{-3}$ in dense mode and 629 points with spacing $6.366 \cdot 10^{-3}$ in sparse mode. Its peaks, from a parabola through the three largest grid values, lie at 2006.70080 and 2006.70092, $1.2 \cdot 10^{-4}$ apart, and its medians are $4 \cdot 10^{-5}$ apart, but the grid argmaxes are 2006.70233 and 2006.69825, $4.1 \cdot 10^{-3}$ apart. On flu the difference of the test is therefore mostly the position of the peak between grid points, and it can reach one grid cell plus the continuous difference.

The trapezoid rule in the convolutions changes how strongly the per-edge difference shows. On `flu_h3n2_20` the largest median difference grows from $4.2 \cdot 10^{-4}$ (`f4bd506d0`) to $7.4 \cdot 10^{-4}$ years (`df43a205d`); on the other five datasets the median differences change by at most 11%, in both directions. Replacing one sparse branch-length distribution at a time by the dense one locates the edges that move the node of the largest difference. At `df43a205d`, the edge to the tip `A/New_Hampshire/12/2012` alone moves the grid argmax of `NODE_0000011` by $4.2 \cdot 10^{-3}$ years but its parabolic peak by $7.6 \cdot 10^{-5}$ and its median by $9 \cdot 10^{-6}$. On the four edges whose replacement moves the grid argmax of `NODE_0000011` most, the two modes differ most at branch length 0, the hard end of the grid (peak-normalized negative log-likelihood differences of $2.8 \cdot 10^{-4}$ to 0.12 at length 0, below $3 \cdot 10^{-6}$ at the peak). The trapezoid rule gives that end sample half the weight of the rectangle sum, so the end weights do not amplify the difference at the hard end.

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

The sparse backend is the default of every run ([N-representation-infer-dense-stub.md](N-representation-infer-dense-stub.md)), and the dense backend is its independent oracle. These issues show the two backends disagreeing on default settings and may share one cause.

- [H-ancestral-sparse-dense-internal-residues-diverge.md](H-ancestral-sparse-dense-internal-residues-diverge.md)
- [M-ancestral-sparse-root-invariance.md](M-ancestral-sparse-root-invariance.md)
