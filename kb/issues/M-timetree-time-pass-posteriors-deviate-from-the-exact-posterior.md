# Time-pass posteriors deviate from the exact posterior by up to 1.4e-2 in CDF

Without a coalescent prior, the marginal time inference is exact message passing on a tree, so its node posteriors should match the exact posterior up to the discretization of the passes. `test_marginal_exact_time_passes_match_the_integrated_posterior_without_a_prior` [packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs#L36-L70](../../packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs#L36-L70) compares them on `((A,B)X,C)R` with tips at 2012, 2010 and 2011, JC69 branches of 1000 sites with known substitution counts, and an oracle that integrates the exact joint density on a 0.005-year grid. The largest CDF difference is $1.4 \cdot 10^{-2}$, at node R with zero substitutions on the branches near the tips; each case asserts its own limit between $10^{-4}$ and $10^{-1}$ per node.

The time passes run with the production discretization: branch distributions of `fn compute_branch_length_distribution()`, backward and forward messages of 300 points (`GRID_POINTS`) with mass windows at `EPS` $= 5 \cdot 10^{-4}$, and the FFT convolution.

> [!IMPORTANT]
> **Investigation required.** Which part of the discretization causes the difference. Candidates, each with its own issue or decision:
>
> - the first-order error of the FFT convolution at the hard edge of a zero-count branch distribution, which jumps from zero to its peak at length 0 ([M-distribution-fft-convolution-hard-edge-first-order-error.md](M-distribution-fft-convolution-hard-edge-first-order-error.md)); the largest difference occurs with zero counts
> - the 300-point messages and their mass windows ([H-timetree-marginal-likelihood-of-node-times-biased-at-300-grid-points.md](H-timetree-marginal-likelihood-of-node-times-biased-at-300-grid-points.md) shows a first-order bias of the same messages)
> - the uniform branch grid ([M-timetree-branch-grid-uniform-resolution.md](M-timetree-branch-grid-uniform-resolution.md))
> - the oracle's own quadrature on the 0.005-year grid
>
> Refine one at a time and record how the difference changes.

The mean-field test of the same file ([packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs#L77-L112](../../packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs#L77-L112)) measures up to $1.62 \cdot 10^{-2}$ against the exact Kingman posterior, close to the error without a prior. The discretization may therefore dominate there too, but no measurement separates it from the mean-field approximation ([M-timetree-marginal-node-posteriors-miscalibrated-at-low-signal.md](M-timetree-marginal-node-posteriors-miscalibrated-at-low-signal.md)).

## Validation

The test without a prior passes at `1e-4` or tighter on every case, and its `TODO` line goes away.
