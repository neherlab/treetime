# Time-pass posteriors on long branches keep a CDF error of 2.5e-4 that more grid points do not remove

Without a coalescent prior, the marginal time inference is exact message passing on a tree, so its node posteriors should match the exact posterior up to the discretization of the passes. `test_marginal_exact_time_passes_match_the_integrated_posterior_without_a_prior` [packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs#L36-L71](../../packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs#L36-L71) compares them on `((A,B)X,C)R` with tips at 2012, 2010 and 2011, JC69 branches of 1000 sites with known substitution counts, and an oracle that integrates the exact joint density on a 0.005-year grid. Each case asserts its own limit per node.

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

Refining the message grid (`GRID_POINTS` 300, 600, 900; 2026-10-24) separates two kinds of difference:

| Case | Node | 300 | 600 | 900 |
| --- | --- | --- | --- | --- |
| both tips at zero | R | $2.93 \cdot 10^{-3}$ | $5.5 \cdot 10^{-4}$ | $2.65 \cdot 10^{-4}$ |
| long branches | X | $2.2 \cdot 10^{-4}$ | $1.5 \cdot 10^{-4}$ | $1.4 \cdot 10^{-4}$ |
| long branches | R | $2.5 \cdot 10^{-4}$ | $2.5 \cdot 10^{-4}$ | $2.5 \cdot 10^{-4}$ |

- **Both tips at zero**: the difference falls by 5.3 and 11 times, at least the factors 4 and 9 of second order. This is the accepted second-order error of the default grid ([discretization error of the default grid](../decisions/distribution-tails-and-arithmetic.md#discretization-error-of-the-default-grid))
- **Long branches**: the difference stays at about $2.5 \cdot 10^{-4}$ at node R. A narrower mass window (`EPS` $5 \cdot 10^{-4}$ to $5 \cdot 10^{-5}$) does not lower it either ($3.1 \cdot 10^{-4}$)

> [!IMPORTANT]
> **Investigation required.** The source of the floor of about $2.5 \cdot 10^{-4}$ on the long-branch case, which neither the message grid nor the mass window sets. Candidates: the range or the uniform spacing of the branch grid ([M-timetree-branch-grid-uniform-resolution.md](M-timetree-branch-grid-uniform-resolution.md)), and the cells at the hard edges of the convolution (anchoring resampled operands at their denser end, 0d1cf2121, leaves this case unchanged: node R $2.865 \cdot 10^{-4}$ before and after; [kb/decisions/distribution-convolution-trapezoid-rule.md](../decisions/distribution-convolution-trapezoid-rule.md#operands-resampled-to-a-common-spacing-keep-their-denser-end-on-a-grid-point)). Refine the branch grid alone and record how the difference changes.

With the mean-field coalescent prior, `test_marginal_exact_mean_field_passes_converge_to_the_exact_mean_field_posterior` ([packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs#L73-L116](../../packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs#L73-L116)) compares the passes on `((A,B)X,(C,D)Y)R` (tips at 2012, 2010.5, 2011.5 and 2010) with the exact posterior of the same factor graph: the JC69 branch likelihoods and the frozen prior factors of `CoalescentPrior`, integrated by the trapezoid rule over the joint times of X, Y and R on a 0.00125-year grid from 1995. The model error of the mean field against the Kingman posterior is in [M-timetree-marginal-node-posteriors-miscalibrated-at-low-signal.md](M-timetree-marginal-node-posteriors-miscalibrated-at-low-signal.md).

- **Oracle error**: the CDFs of the oracle converge at second order (the difference to a 0.000625-year oracle falls 4.2 to 5.0 times per halving of the step from 0.005 years). The Richardson estimate of the error of the 0.00125-year oracle is at most $4.5 \cdot 10^{-6}$ (node Y, short internal branches, $T_c = 1$). A grid start in 2003 cut off $2.5 \cdot 10^{-4}$ of the root mass at $T_c = 10$ on the overlapping clades; from 1995 the cut is below $10^{-20}$
- **Assertions**: the largest CDF difference per node must fall at least 3 times per halving of the spacing of the branch grids (300, 599 and 1197 points), and the oracle error must stay below one tenth of the difference at 1197 points. A prior factor scaled by 1.001 in the passes fails every case
- **Convergence**: on `rust` at f452d6572 every node converges to the oracle, with factors 3.10 to 8.56 per halving. The smallest factors are those of node X of the overlapping clades at $T_c = 10$ (3.10 and 3.46); between 1197 and 2393 points the differences of that case level off near $5 \cdot 10^{-6}$ to $1 \cdot 10^{-5}$

Largest CDF difference to the oracle per case and node at 300, 599 and 1197 points of the branch grid (counts on the branches above A, B, C, D, X and Y):

| Case, $T_c$ | Node | 300 | 599 | 1197 |
| --- | --- | --- | --- | --- |
| short internal branches (3, 0, 2, 1, 1, 2), 1 | X | $3.99 \cdot 10^{-4}$ | $9.01 \cdot 10^{-5}$ | $1.50 \cdot 10^{-5}$ |
| short internal branches (3, 0, 2, 1, 1, 2), 1 | Y | $1.68 \cdot 10^{-2}$ | $2.70 \cdot 10^{-3}$ | $5.12 \cdot 10^{-4}$ |
| short internal branches (3, 0, 2, 1, 1, 2), 1 | R | $3.32 \cdot 10^{-3}$ | $7.58 \cdot 10^{-4}$ | $1.24 \cdot 10^{-4}$ |
| short internal branches (3, 0, 2, 1, 1, 2), 10 | X | $3.17 \cdot 10^{-4}$ | $7.38 \cdot 10^{-5}$ | $1.56 \cdot 10^{-5}$ |
| short internal branches (3, 0, 2, 1, 1, 2), 10 | Y | $6.65 \cdot 10^{-3}$ | $1.01 \cdot 10^{-3}$ | $1.35 \cdot 10^{-4}$ |
| short internal branches (3, 0, 2, 1, 1, 2), 10 | R | $1.49 \cdot 10^{-3}$ | $3.31 \cdot 10^{-4}$ | $5.14 \cdot 10^{-5}$ |
| overlapping clades (6, 3, 4, 3, 3, 4), 1 | X | $3.82 \cdot 10^{-4}$ | $8.11 \cdot 10^{-5}$ | $1.95 \cdot 10^{-5}$ |
| overlapping clades (6, 3, 4, 3, 3, 4), 1 | Y | $6.35 \cdot 10^{-4}$ | $1.05 \cdot 10^{-4}$ | $2.63 \cdot 10^{-5}$ |
| overlapping clades (6, 3, 4, 3, 3, 4), 1 | R | $3.91 \cdot 10^{-4}$ | $1.03 \cdot 10^{-4}$ | $2.39 \cdot 10^{-5}$ |
| overlapping clades (6, 3, 4, 3, 3, 4), 10 | X | $1.08 \cdot 10^{-4}$ | $3.48 \cdot 10^{-5}$ | $1.01 \cdot 10^{-5}$ |
| overlapping clades (6, 3, 4, 3, 3, 4), 10 | Y | $1.47 \cdot 10^{-4}$ | $3.13 \cdot 10^{-5}$ | $6.18 \cdot 10^{-6}$ |
| overlapping clades (6, 3, 4, 3, 3, 4), 10 | R | $1.80 \cdot 10^{-4}$ | $5.09 \cdot 10^{-5}$ | $9.33 \cdot 10^{-6}$ |
| long internal branches (2, 1, 1, 2, 3, 3), 1 | X | $6.22 \cdot 10^{-3}$ | $1.37 \cdot 10^{-3}$ | $1.61 \cdot 10^{-4}$ |
| long internal branches (2, 1, 1, 2, 3, 3), 1 | Y | $6.63 \cdot 10^{-4}$ | $1.00 \cdot 10^{-4}$ | $1.89 \cdot 10^{-5}$ |
| long internal branches (2, 1, 1, 2, 3, 3), 1 | R | $8.08 \cdot 10^{-4}$ | $1.84 \cdot 10^{-4}$ | $3.17 \cdot 10^{-5}$ |

> [!IMPORTANT]
> **Investigation required.** At the default grid the error reaches $1.68 \cdot 10^{-2}$ (node Y, short internal branches, $T_c = 1$) and $6.2 \cdot 10^{-3}$ (node X, long internal branches, $T_c = 1$), above every case without a prior. Both errors fall at second order or faster, so they are discretization error, but their size at 300 points is not explained. Find which spacing sets them (branch grid, message grid, or the sampling of the prior factor on the message grid) and record whether the default grid needs more points where a hard tip date bounds a short branch.

## Validation

The long-branch difference falls with the grid at second order, or its source is found and accepted. The limits of the test then follow the second-order error of the default grid, and its `TODO` lines go away.
