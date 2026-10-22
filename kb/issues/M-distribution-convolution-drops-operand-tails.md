# Convolution drops the operands' extrapolated tails

The FFT convolution of two functions, `fn convolution_function_function_fine()` [packages/treetime-distribution/src/distribution_ops/convolve.rs#L260-L323](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L260-L323), convolves grid values only. The operands' `Linear` and `HardApproach` tails beyond their grids are dropped, and the edge convolution (`fn WindowedConvolution::new()` in [packages/treetime-distribution/src/distribution_ops/edge_convolution.rs#L220-L265](../../packages/treetime-distribution/src/distribution_ops/edge_convolution.rs#L220-L265)) refits the result's soft tail from its outermost grid points. The convolutions of a function with a point or a date range are exact and keep the tail laws ([kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md#exact-point-and-range-convolution)), so every edge between two function distributions remains affected.

Branch-length distributions are windowed at $\epsilon = 5 \cdot 10^{-4}$ of the mass per side (`fn compute_branch_length_distribution()` calls `fn rewindow_to_mass()` [packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L63](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L63), defined in [packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L28-L42](../../packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L28-L42); the left side is `Hard` or `HardApproach` and keeps its grid edge, so only the right tail is cut). For a zero-mutation branch (exponential with rate $\lambda$, grid cut at $\ln(1/\epsilon)/\lambda$) the dropped tail carries about 3% of the variance and 0.4% of the mean, so each convolution adds systematically too little variance.

> [!IMPORTANT]
> **Investigation required.** Whether the variance loss changes results: confidence-interval coverage of node dates on the validation suite for runs dominated by short branches (low mu * L), with the tails included against without.

Investigated 2026-10-19 (prototype commit `6d06d211a`, extending `Linear` tails onto the grid): no gain in low-signal coverage (0.830 to 0.817 at $L=10^3$, $\mu=3\cdot10^{-4}$), interval widths within 1% at a fixed rate; see [kb/proposals/timetree-forward-pass-grid-support.md](../proposals/timetree-forward-pass-grid-support.md#operand-tails-issue-3).

## Effect of the trapezoid rule at soft ends

The convolution integrates each overlap by the trapezoid rule ([kb/decisions/distribution-convolution-trapezoid-rule.md](../decisions/distribution-convolution-trapezoid-rule.md)). At a soft grid end the operand continues beyond the grid, and the convolution drops that tail. The rectangle sum gave the end sample an extra half weight, which partly offset the dropped mass; the trapezoid rule removes that offset. In the analytic validation suite of `treetime-validation`, Gaussians cut at a soft grid end are 1% to 26% worse in relative L2 error with the trapezoid rule than with the rectangle sum (identical for every convolution algorithm):

| suite | case | points | domain | relative L2, rectangle | trapezoid | mass error, rectangle | trapezoid |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `conv-gaussian-gaussian` | `python_notebook_case_1` | 101 | $[-5, 5]$ | 3.01e-5 | 3.80e-5 | 2.59e-5 | 3.27e-5 |
| `conv-gaussian-gaussian` | `python_notebook_case_2` | 101 | $[-5, 5]$ | 1.14e-2 | 1.23e-2 | 1.16e-2 | 1.24e-2 |
| `conv-gaussian-gaussian` | `tight_truncation` | 1201 | $[-6, 6]$ | 2.67e-3 | 2.70e-3 | 2.68e-3 | 2.70e-3 |
| `conv-gaussian-gaussian` | `coarse_grid` | 161 | $[-8, 8]$ | 2.6e-14 | 4.1e-14 | 2.5e-14 | 3.8e-14 |
| `conv-gaussian-exponential` | `tight_truncation` | 201 | $[-2, 15]$ | 2.32e-2 | 2.58e-2 | 1.99e-2 | 2.22e-2 |

The other Gaussian cases are at round-off ($10^{-14}$ to $10^{-11}$) and unchanged. Extending the soft tails onto the grid (see Fix direction) removes the dropped mass, and with it the difference between the two rules.

## Fix direction

Options, in order of scope:

- Extend each operand's soft tail onto the grid for a fixed number of decay lengths before the FFT
- Convolve exponential tails by the first-order recursion $y_i = e^{-k\,dx}\,y_{i-1} + w_0 q_{i-1} + w_1 q_i$ (exact for the piecewise-linear interpolant, $O(n)$ at any $k\,dx$)
- Split operands into a gridded core and analytic tail pieces and sum the partial convolutions

Closed-form tail moments (`TailMoments`, `fn SoftTailLaw::moments()`, `fn HardApproachLaw::moments()`, quadrature-tested) exist on the unmerged branch `moment-matched-convolution`.

## Related

- [kb/proposals/timetree-forward-pass-grid-support.md](../proposals/timetree-forward-pass-grid-support.md)
- [kb/proposals/exponential-branch-length-tails.md](../proposals/exponential-branch-length-tails.md)
