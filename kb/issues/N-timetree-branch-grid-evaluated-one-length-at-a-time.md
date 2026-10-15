# Timetree branch-length grid is evaluated one branch length at a time

For every edge in every timetree iteration, `fn evaluate_log_lh_on_grid()` evaluates the branch-length log-likelihood at 300 grid points, one call per point [packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L86-L96](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L86-L96). Each call walks all site coefficients of the edge again, so every coefficient row is read 300 times per edge.

## Locations

- `fn evaluate_log_lh_on_grid()` maps each grid point through `fn evaluate_with_indels_log_lh_only()` [packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L94](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L94)
- `fn evaluate_site_contributions()` computes the exponentials of one branch length and loops over all sites [packages/treetime/src/optimize/eval.rs](../../packages/treetime/src/optimize/eval.rs)

## Proposed solution

Evaluate the grid as one batch per edge: build the matrix $E_{j,i} = \exp(\mu \lambda_i t_j)$ once (300 rows, one column per alphabet state), then loop over sites on the outside and grid points on the inside, with `E.row(j)` as the exponentials of grid point $j$. Each coefficient row is then read once per edge.

The results stay bit-identical when, for each grid point, the sites are still added in the same order and the dot products stay on the non-BLAS path of ndarray (fewer than 32 states). A matrix product through BLAS would be faster but changes the summation order, so outputs would change in the last bits.

> [!IMPORTANT]
> **Investigation required.** The speed-up is unmeasured. Measure dense and sparse timetree on large datasets (for example `sc2/4500`, `dengue/2000`) against the current one-point-at-a-time loop.
