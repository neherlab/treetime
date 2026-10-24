# Time-pass posteriors deviate from the exact posterior by up to 2.9e-3 in CDF

Without a coalescent prior, the marginal time inference is exact message passing on a tree, so its node posteriors should match the exact posterior up to the discretization of the passes. `test_marginal_exact_time_passes_match_the_integrated_posterior_without_a_prior` [packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs#L36-L70](../../packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs#L36-L70) compares them on `((A,B)X,C)R` with tips at 2012, 2010 and 2011, JC69 branches of 1000 sites with known substitution counts, and an oracle that integrates the exact joint density on a 0.005-year grid. Each case asserts its own limit per node.

The time passes run with the production discretization: branch distributions of `fn compute_branch_length_distribution()`, backward and forward messages of 300 points (`GRID_POINTS`) with mass windows at `EPS` $= 5 \cdot 10^{-4}$, and the convolution of two functions by the trapezoid rule over each overlap.

The oracle integrates the exact joint density by the trapezoid rule on its grid, with every hard end (the dates of the tips and the zero-length branch above X) on a grid point. Its CDFs converge at second order in its step (observed order 2.000 for every case and node between the steps 0.005, 0.0025 and 0.00125 years), and the error of the 0.005-year oracle of the test, estimated from the Richardson extrapolation of the three steps, is at most $1.6 \cdot 10^{-5}$ (node R, zero counts near the tips), below one tenth of every limit.

Largest CDF difference per case and node to the oracle of the test, with the trapezoid rule in the convolutions and, in parentheses, with rectangle sums in the convolutions:

| Case (counts A, B, X, C) | Node | Difference | Oracle error | Limit |
| --- | --- | --- | --- | --- |
| mixed counts (3, 0, 2, 6) | X | $2.3 \cdot 10^{-5}$ ($2.3 \cdot 10^{-5}$) | $1.6 \cdot 10^{-6}$ | `1e-4` |
| mixed counts (3, 0, 2, 6) | R | $7.2 \cdot 10^{-5}$ ($1.9 \cdot 10^{-3}$) | $3.8 \cdot 10^{-6}$ | `1e-4` |
| zero counts near the tips (1, 0, 0, 2) | X | $1.3 \cdot 10^{-4}$ ($4.8 \cdot 10^{-4}$) | $2.5 \cdot 10^{-6}$ | `1e-3` |
| zero counts near the tips (1, 0, 0, 2) | R | $2.0 \cdot 10^{-4}$ ($1.4 \cdot 10^{-2}$) | $1.6 \cdot 10^{-5}$ | `1e-3` |
| long branches (8, 4, 5, 12) | X | $2.2 \cdot 10^{-4}$ ($2.2 \cdot 10^{-4}$) | $1.1 \cdot 10^{-6}$ | `1e-3` |
| long branches (8, 4, 5, 12) | R | $2.9 \cdot 10^{-4}$ ($2.9 \cdot 10^{-4}$) | $4.7 \cdot 10^{-7}$ | `1e-3` |
| both tips at zero (0, 0, 1, 1) | X | $6.2 \cdot 10^{-5}$ ($6.2 \cdot 10^{-5}$) | $1.8 \cdot 10^{-7}$ | `1e-4` |
| both tips at zero (0, 0, 1, 1) | R | $2.9 \cdot 10^{-3}$ ($3.7 \cdot 10^{-3}$) | $1.2 \cdot 10^{-5}$ | `1e-2` |

The trapezoid rule makes no case worse: it removes most of the error at node R with mixed counts and with zero counts near the tips, and lowers node R with both tips at zero from $3.7 \cdot 10^{-3}$ to $2.9 \cdot 10^{-3}$. The limits of node R with mixed counts and with zero counts near the tips fail with rectangle sums.

> [!IMPORTANT]
> **Investigation required.** Which part of the discretization of the time passes causes the remaining differences, the largest at node R with both tips at zero ($2.9 \cdot 10^{-3}$), where both tip branches have zero counts and the subtree density of X jumps from zero to its maximum at the date of B. Candidates:
>
> - the end cells of the trapezoid convolution at the hard edges of zero-count branch distributions, which jump from zero to their peak at length 0
> - the 300-point messages and their mass windows (their error in the log marginal likelihood is second order: [discretization error of the default grid](../decisions/distribution-tails-and-arithmetic.md#discretization-error-of-the-default-grid))
> - the uniform branch grid ([M-timetree-branch-grid-uniform-resolution.md](M-timetree-branch-grid-uniform-resolution.md))
>
> Refine one at a time and record how the difference to the converged oracle changes.

The mean-field test of the same file ([packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs#L77-L112](../../packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs#L77-L112)) measures up to $1.62 \cdot 10^{-2}$ against the exact Kingman posterior (node X, overlapping clades, strong prior). Halving the step of its 0.02-year oracle changes the oracle CDFs by $1.1 \cdot 10^{-5}$ to $7.1 \cdot 10^{-4}$, and the differences change by less than $5 \cdot 10^{-4}$ against the finer oracle, so the oracle does not explain them. They are the error of the mean-field approximation together with the discretization ([M-timetree-marginal-node-posteriors-miscalibrated-at-low-signal.md](M-timetree-marginal-node-posteriors-miscalibrated-at-low-signal.md)).

## Validation

Each remaining difference falls below `1e-4` against the second-order oracle of the test, and the test passes with a limit of `1e-4` on every case and node. Then its `TODO` lines go away.
