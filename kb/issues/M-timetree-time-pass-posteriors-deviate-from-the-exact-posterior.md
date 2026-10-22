# Time-pass posteriors deviate from the exact posterior by up to 5.3e-3 in CDF

Without a coalescent prior, the marginal time inference is exact message passing on a tree, so its node posteriors should match the exact posterior up to the discretization of the passes. `test_marginal_exact_time_passes_match_the_integrated_posterior_without_a_prior` [packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs#L36-L70](../../packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs#L36-L70) compares them on `((A,B)X,C)R` with tips at 2012, 2010 and 2011, JC69 branches of 1000 sites with known substitution counts, and an oracle that integrates the exact joint density on a 0.005-year grid. Each case asserts its own limit per node.

The time passes run with the production discretization: branch distributions of `fn compute_branch_length_distribution()`, backward and forward messages of 300 points (`GRID_POINTS`) with mass windows at `EPS` $= 5 \cdot 10^{-4}$, and the convolution of two functions by the trapezoid rule over each overlap.

Largest CDF difference per case and node, the change of the oracle CDF when its grid step is halved to 0.0025 years, and the difference to that finer oracle:

| Case (counts A, B, X, C) | Node | Difference | Oracle change at half step | Difference to the half-step oracle | Limit |
| --- | --- | --- | --- | --- | --- |
| mixed counts (3, 0, 2, 6) | X | $2.3 \cdot 10^{-5}$ | $1.2 \cdot 10^{-6}$ | $2.4 \cdot 10^{-5}$ | `1e-4` |
| mixed counts (3, 0, 2, 6) | R | $8.6 \cdot 10^{-4}$ | $4.0 \cdot 10^{-4}$ | $4.7 \cdot 10^{-4}$ | `1e-2` |
| zero counts near the tips (1, 0, 0, 2) | X | $2.4 \cdot 10^{-4}$ | $7.9 \cdot 10^{-5}$ | $1.7 \cdot 10^{-4}$ | `1e-3` |
| zero counts near the tips (1, 0, 0, 2) | R | $2.0 \cdot 10^{-4}$ | $2.2 \cdot 10^{-5}$ | $2.0 \cdot 10^{-4}$ | `1e-3` |
| long branches (8, 4, 5, 12) | X | $2.2 \cdot 10^{-4}$ | $8.9 \cdot 10^{-7}$ | $2.2 \cdot 10^{-4}$ | `1e-3` |
| long branches (8, 4, 5, 12) | R | $2.9 \cdot 10^{-4}$ | $3.6 \cdot 10^{-7}$ | $2.9 \cdot 10^{-4}$ | `1e-3` |
| both tips at zero (0, 0, 1, 1) | X | $6.2 \cdot 10^{-5}$ | $1.3 \cdot 10^{-7}$ | $6.2 \cdot 10^{-5}$ | `1e-4` |
| both tips at zero (0, 0, 1, 1) | R | $5.3 \cdot 10^{-3}$ | $1.3 \cdot 10^{-3}$ | $4.1 \cdot 10^{-3}$ | `1e-2` |

The trapezoid rule of the convolution removed the largest difference of the rectangle sum: zero counts near the tips at node R fell from $1.4 \cdot 10^{-2}$ to $2.0 \cdot 10^{-4}$, and mixed counts at node R from $1.1 \cdot 10^{-3}$ to $8.6 \cdot 10^{-4}$. With both tips at zero counts, node R rose from $1.3 \cdot 10^{-3}$ to $5.3 \cdot 10^{-3}$ (from $2.5 \cdot 10^{-3}$ to $4.1 \cdot 10^{-3}$ against the half-step oracle).

> [!IMPORTANT]
> **Investigation required.** Which part of the discretization causes the remaining differences at node R. Candidates:
>
> - the oracle's own quadrature: at mixed counts and at both tips at zero, halving its step changes its CDF by $4.0 \cdot 10^{-4}$ and $1.3 \cdot 10^{-3}$, the same order as the differences. The oracle leaves out the zero-length branch above X (`j < i` in `fn exact_marginals()`) and sums on a finite grid; a converged oracle (Richardson extrapolation over two step halvings, or an exact term for the zero-length boundary) must come first
> - the end cells of the trapezoid convolution at the hard edges of zero-count branch distributions, which jump from zero to their peak at length 0: the case with both tips at zero counts is the only one that got worse with the trapezoid rule
> - the 300-point messages and their mass windows ([H-timetree-marginal-likelihood-of-node-times-biased-at-300-grid-points.md](H-timetree-marginal-likelihood-of-node-times-biased-at-300-grid-points.md) shows a first-order bias of the same messages)
> - the uniform branch grid ([M-timetree-branch-grid-uniform-resolution.md](M-timetree-branch-grid-uniform-resolution.md))
>
> Refine one at a time and record how the difference changes.

The mean-field test of the same file ([packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs#L77-L112](../../packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs#L77-L112)) measures up to $1.62 \cdot 10^{-2}$ against the exact Kingman posterior (node X, overlapping clades, strong prior). Halving the step of its 0.02-year oracle changes the oracle CDFs by $1.1 \cdot 10^{-5}$ to $7.1 \cdot 10^{-4}$, and the differences change by less than $5 \cdot 10^{-4}$ against the finer oracle, so the oracle does not explain them. They are the error of the mean-field approximation together with the discretization ([M-timetree-marginal-node-posteriors-miscalibrated-at-low-signal.md](M-timetree-marginal-node-posteriors-miscalibrated-at-low-signal.md)).

## Validation

The test without a prior passes at `1e-4` or tighter on every case, and its `TODO` line goes away.
