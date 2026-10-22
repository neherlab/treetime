# Time-pass posteriors deviate from the exact posterior by up to 2.9e-3 in CDF

Without a coalescent prior, the marginal time inference is exact message passing on a tree, so its node posteriors should match the exact posterior up to the discretization of the passes. `test_marginal_exact_time_passes_match_the_integrated_posterior_without_a_prior` [packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs#L36-L70](../../packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs#L36-L70) compares them on `((A,B)X,C)R` with tips at 2012, 2010 and 2011, JC69 branches of 1000 sites with known substitution counts, and an oracle that integrates the exact joint density on a 0.005-year grid. Each case asserts its own limit per node.

The time passes run with the production discretization: branch distributions of `fn compute_branch_length_distribution()`, backward and forward messages of 300 points (`GRID_POINTS`) with mass windows at `EPS` $= 5 \cdot 10^{-4}$, and the convolution of two functions by the trapezoid rule over each overlap.

The oracle sums the inner integral of the exact density on its own grid: it leaves out the zero-length branch above X (`j < i` in `fn exact_marginals()`) and gives the sample at a hard end its full weight. Its CDFs therefore converge at first order in its step at most nodes (observed order 1.00 to 1.04 at the finest steps), and at second order at node X with mixed counts and with both tips at zero and at node R with long branches. A converged oracle comes from halving its step seven times (0.005 to $3.9 \cdot 10^{-5}$ years) and a Richardson extrapolation with the observed order; the extrapolation moves the CDFs of the finest step by at most $2.0 \cdot 10^{-5}$.

Largest CDF difference per case and node to the converged oracle and to the 0.005-year oracle of the test, at `df43a205d` (trapezoid rule in the convolutions) and, in parentheses, on `rust` at `f4bd506d0` (rectangle sums):

| Case (counts A, B, X, C) | Node | Converged oracle | Oracle of the test (step 0.005) | Limit |
| --- | --- | --- | --- | --- |
| mixed counts (3, 0, 2, 6) | X | $2.5 \cdot 10^{-5}$ ($2.5 \cdot 10^{-5}$) | $2.3 \cdot 10^{-5}$ | `1e-4` |
| mixed counts (3, 0, 2, 6) | R | $6.9 \cdot 10^{-5}$ ($1.9 \cdot 10^{-3}$) | $8.6 \cdot 10^{-4}$ ($1.1 \cdot 10^{-3}$) | `1e-2` |
| zero counts near the tips (1, 0, 0, 2) | X | $1.3 \cdot 10^{-4}$ ($4.8 \cdot 10^{-4}$) | $2.4 \cdot 10^{-4}$ ($6.4 \cdot 10^{-4}$) | `1e-3` |
| zero counts near the tips (1, 0, 0, 2) | R | $2.0 \cdot 10^{-4}$ ($1.4 \cdot 10^{-2}$) | $2.0 \cdot 10^{-4}$ ($1.4 \cdot 10^{-2}$) | `1e-3` |
| long branches (8, 4, 5, 12) | X | $2.2 \cdot 10^{-4}$ ($2.2 \cdot 10^{-4}$) | $2.2 \cdot 10^{-4}$ | `1e-3` |
| long branches (8, 4, 5, 12) | R | $2.9 \cdot 10^{-4}$ ($2.9 \cdot 10^{-4}$) | $2.9 \cdot 10^{-4}$ | `1e-3` |
| both tips at zero (0, 0, 1, 1) | X | $6.1 \cdot 10^{-5}$ ($6.1 \cdot 10^{-5}$) | $6.2 \cdot 10^{-5}$ | `1e-4` |
| both tips at zero (0, 0, 1, 1) | R | $2.9 \cdot 10^{-3}$ ($3.7 \cdot 10^{-3}$) | $5.3 \cdot 10^{-3}$ ($1.3 \cdot 10^{-3}$) | `1e-2` |

Against the converged oracle the trapezoid rule makes no case worse: it removes most of the error at node R with mixed counts and with zero counts near the tips, and lowers node R with both tips at zero from $3.7 \cdot 10^{-3}$ to $2.9 \cdot 10^{-3}$. Against the oracle of the test, node R with both tips at zero appears to rise from $1.3 \cdot 10^{-3}$ to $5.3 \cdot 10^{-3}$: the first-order error of that oracle is $2.4 \cdot 10^{-3}$ there, and it has the opposite sign of the error of the rectangle sums, so it hid most of that error.

> [!IMPORTANT]
> **Investigation required.** Which part of the discretization of the time passes causes the remaining differences, the largest at node R with both tips at zero ($2.9 \cdot 10^{-3}$), where both tip branches have zero counts and the subtree density of X jumps from zero to its maximum at the date of B. Candidates:
>
> - the end cells of the trapezoid convolution at the hard edges of zero-count branch distributions, which jump from zero to their peak at length 0
> - the 300-point messages and their mass windows ([H-timetree-marginal-likelihood-of-node-times-biased-at-300-grid-points.md](H-timetree-marginal-likelihood-of-node-times-biased-at-300-grid-points.md) shows a first-order bias of the same messages)
> - the uniform branch grid ([M-timetree-branch-grid-uniform-resolution.md](M-timetree-branch-grid-uniform-resolution.md))
>
> Refine one at a time and record how the difference to the converged oracle changes.

The mean-field test of the same file ([packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs#L77-L112](../../packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs#L77-L112)) measures up to $1.62 \cdot 10^{-2}$ against the exact Kingman posterior (node X, overlapping clades, strong prior). Halving the step of its 0.02-year oracle changes the oracle CDFs by $1.1 \cdot 10^{-5}$ to $7.1 \cdot 10^{-4}$, and the differences change by less than $5 \cdot 10^{-4}$ against the finer oracle, so the oracle does not explain them. They are the error of the mean-field approximation together with the discretization ([M-timetree-marginal-node-posteriors-miscalibrated-at-low-signal.md](M-timetree-marginal-node-posteriors-miscalibrated-at-low-signal.md)).

## Validation

The test compares with a converged oracle. The first-order extrapolation $2 C_{h/4} - C_{h/2}$ from the oracle CDFs at the steps 0.0025 and 0.00125 years lies within $3.3 \cdot 10^{-6}$ of the converged oracle for every case and node, and needs about 20 times the arithmetic of the current oracle (the double sum grows with the square of the point count). With it, the test passes at `1e-4` or tighter on every case, and its `TODO` line goes away.
