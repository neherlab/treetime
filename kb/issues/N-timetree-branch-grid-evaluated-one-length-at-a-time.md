# Timetree branch-length grid is evaluated one branch length at a time

For every edge in every timetree iteration, `fn evaluate_log_lh_on_grid()` evaluates the branch-length log-likelihood at 300 grid points, one call per point [packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L86-L96](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L86-L96). Each call walks all site coefficients of the edge again, so every coefficient row is read 300 times per edge.

## Cost

On the validation suite's `const/rep00/L10000_mu0.01` (400 tips, $L=10^4$, fixed clock, `--keep-root`), building the branch distributions from the alignment took about 4 s per timetree round, about 35% of the wall time, while the forward pass took about 55% ([M-timetree-forward-grids-keep-untrimmed-hard-support.md](M-timetree-forward-grids-keep-untrimmed-hard-support.md)).

Measured 2026-10-19 with the forward pass trimmed (prototype commit `6d06d211a`; 400 tips, $L=10^4$, $\mu=10^{-2}$, fixed clock): 93% of the branch-distribution step is `fn evaluate_log_lh_on_grid()`, 5-6% `fn create_edge_contribution()`, and about 1% the mass re-windowing. The step takes 0.81 s per round at $\mu=10^{-2}$ against 0.11 s at $\mu=10^{-3}$, half of the run time once the forward pass is trimmed ([kb/proposals/timetree-forward-pass-grid-support.md](../proposals/timetree-forward-pass-grid-support.md#branch-distributions-issue-4)). The grid evaluation is therefore the largest cost after the forward-pass fix.

## Locations

- `fn evaluate_log_lh_on_grid()` maps each grid point through `fn evaluate_with_indels_log_lh_only()` [packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L94](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L94); its caller is [packages/treetime/src/timetree/inference/runner.rs#L274](../../packages/treetime/src/timetree/inference/runner.rs#L274)
- `fn evaluate_site_contributions()` computes the exponentials of one branch length and loops over all sites [packages/treetime/src/optimize/eval.rs](../../packages/treetime/src/optimize/eval.rs)
- `fn create_edge_contribution()` builds the per-edge contributions of the dense and sparse partitions [packages/treetime/src/partition/marginal/dense/partition.rs#L147](../../packages/treetime/src/partition/marginal/dense/partition.rs#L147) [packages/treetime/src/partition/marginal/sparse/partition.rs#L88](../../packages/treetime/src/partition/marginal/sparse/partition.rs#L88)

## Proposed solution

Evaluate the grid as one batch per edge: build the matrix $E_{j,i} = \exp(\mu \lambda_i t_j)$ once (300 rows, one column per alphabet state), then loop over sites on the outside and grid points on the inside, with `E.row(j)` as the exponentials of grid point $j$. Each coefficient row is then read once per edge.

The results stay bit-identical when, for each grid point, the sites are still added in the same order and the dot products stay on the non-BLAS path of ndarray (fewer than 32 states). A matrix product through BLAS would be faster but changes the summation order, so outputs would change in the last bits.

## Validation

- Outputs of dense and sparse timetree runs are bit-identical before and after the change
- Timing of dense and sparse timetree on large datasets (for example `sc2/4500`, `dengue/2000`) against the one-point-at-a-time loop, with the forward pass trimmed and untrimmed
