# Convolution drops the operands' extrapolated tails

`fn convolution_function_function_fine()` in [packages/treetime-distribution/src/distribution_ops/convolve.rs](../../packages/treetime-distribution/src/distribution_ops/convolve.rs) convolves grid values only. The operands' `Linear` and `HardApproach` tails beyond their grids are dropped, and `fn convolve_across_edge()` refits the result's soft tail from its outermost grid points.

Branch-length distributions are windowed at $\epsilon = 5 \cdot 10^{-4}$ of the mass per side (`fn rewindow_to_mass()` in `branch_length_likelihood.rs`). For a zero-mutation branch (exponential with rate $\lambda$, grid cut at $\ln(1/\epsilon)/\lambda$) the dropped tail carries about 3% of the variance and 0.4% of the mean, so each convolution adds systematically too little variance.

> [!IMPORTANT]
> **Investigation required.** Whether the variance loss changes results: confidence-interval coverage of node dates on the validation suite for runs dominated by short branches (low mu × L), with the tails included against without.

Investigated 2026-10-19 (prototype extending `Linear` tails onto the grid): no gain in low-signal coverage (0.830 to
0.817 at $L=10^3$, $\mu=3\cdot10^{-4}$), interval widths within 1% at a fixed rate; see [the proposal's
results](../proposals/timetree-forward-pass-grid-support.md#operand-tails-issue-3).

## Fix direction

Options, in order of scope:

- Extend each operand's soft tail onto the grid for a fixed number of decay lengths before the FFT
- Convolve exponential tails by the first-order recursion $y_i = e^{-k\,dx}\,y_{i-1} + w_0 q_{i-1} + w_1 q_i$ (exact for the piecewise-linear interpolant, $O(n)$ at any $k\,dx$)
- Split operands into a gridded core and analytic tail pieces and sum the partial convolutions

Closed-form tail moments (`TailMoments`, `fn SoftTailLaw::moments()`, `fn HardApproachLaw::moments()`, quadrature-tested) exist on the unmerged branch `moment-matched-convolution`.

## Related

- [kb/proposals/timetree-forward-pass-grid-support.md](../proposals/timetree-forward-pass-grid-support.md)
- [kb/proposals/exponential-branch-length-tails.md](../proposals/exponential-branch-length-tails.md)
