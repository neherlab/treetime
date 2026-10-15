# Building branch distributions from the alignment dominates once the forward pass is fixed

On the validation suite's `const/rep00/L10000_mu0.01` (400 tips, fixed clock, `--keep-root`, base `94f921f8e7`), computing branch distributions from the alignment takes about 4 s per timetree round, about 35% of the wall time. The forward pass takes about 55% (see [M-timetree-forward-grids-keep-untrimmed-hard-support.md](M-timetree-forward-grids-keep-untrimmed-hard-support.md)); once that is fixed, branch distributions are the largest cost. An earlier profile found this step about 10 times more expensive at mu = 0.01 than at mu = 0.001 for L = 10⁴.

> [!IMPORTANT]
> **Investigation required.** Profile `fn compute_branch_length_distribution()` and its callers in [packages/treetime/src/timetree/inference/branch_length_likelihood.rs](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs): split the time between evaluating the per-site likelihood contributions on the 300-point branch-length grid, the sparse/dense message computation that feeds them, and the mass re-windowing; and measure how each scales with L × mu.
