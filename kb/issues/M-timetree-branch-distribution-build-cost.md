# Building branch distributions from the alignment dominates once the forward pass is fixed

On the validation suite's `const/rep00/L10000_mu0.01` (400 tips, fixed clock, `--keep-root`, base `94f921f8e7`), computing branch distributions from the alignment takes about 4 s per timetree round, about 35% of the wall time. The forward pass takes about 55% (see [M-timetree-forward-grids-keep-untrimmed-hard-support.md](M-timetree-forward-grids-keep-untrimmed-hard-support.md)); once that is fixed, branch distributions are the largest cost. An earlier profile found this step about 10 times more expensive at mu = 0.01 than at mu = 0.001 for L = 10⁴.

> [!IMPORTANT]
> **Investigation required.** Profile `fn compute_branch_length_distribution()` [packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L17-L65](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L17-L65) and its caller [packages/treetime/src/timetree/inference/runner.rs#L274](../../packages/treetime/src/timetree/inference/runner.rs#L274): split the time between evaluating the per-site likelihood contributions on the 300-point branch-length grid, the sparse/dense message computation that feeds them, and the mass re-windowing; and measure how each scales with L * mu.

Measured 2026-10-19 with the forward pass trimmed (prototype commit `6d06d211a`; 400 tips, $L=10^4$, $\mu=10^{-2}$, fixed clock): 93% of the step is `fn evaluate_log_lh_on_grid()`, 5-6% `fn create_edge_contribution()`, about 1% re-windowing; 0.81 s per round at $\mu=10^{-2}$ against 0.11 s at $\mu=10^{-3}$, half of the run time once the forward pass is fixed ([kb/proposals/timetree-forward-pass-grid-support.md](../proposals/timetree-forward-pass-grid-support.md#branch-distributions-issue-4)).

## Locations

- `fn evaluate_log_lh_on_grid()` [packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L82-L92](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L82-L92): evaluates the per-site log-likelihood once per grid branch length
- `fn create_edge_contribution()` [packages/treetime/src/partition/marginal/dense/partition.rs#L147](../../packages/treetime/src/partition/marginal/dense/partition.rs#L147) and [packages/treetime/src/partition/marginal/sparse/partition.rs#L88](../../packages/treetime/src/partition/marginal/sparse/partition.rs#L88): the per-edge contributions of the dense and sparse partitions
- `fn rewindow_to_mass()` [packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L23-L37](../../packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L23-L37): the mass re-windowing
