# Probability Distributions

## Implemented

- [x] Analytical Gaussian (PDF, product, convolution)
- [x] Analytical Exponential (PDF, convolution, special case a=b)
- [x] Gaussian-Exponential convolution
- [/] Grid distributions (uniform grid and linear interpolation implemented; support-boundary semantics unresolved - [kb/issues/M-distribution-support-boundary-semantics-unresolved.md](../issues/M-distribution-support-boundary-semantics-unresolved.md))
- [x] ScaledArray pattern (normalized values + log-scale factor)

## Partial or Not Implemented

- [ ] Unified Distribution class (v0: wraps scipy.interpolate.interp1d)
- [ ] Delta functions (point masses, v0 `Distribution.delta_function()`)
- [/] Distribution multiplication (implemented; [kb/issues/M-distribution-mixed-support-operations-lose-exact-intersection-boundaries.md](../issues/M-distribution-mixed-support-operations-lose-exact-intersection-boundaries.md), [kb/issues/M-distribution-product-grid-resolution-diverges-from-v0.md](../issues/M-distribution-product-grid-resolution-diverges-from-v0.md))
- [/] Distribution division (implemented; [kb/issues/M-distribution-mixed-support-operations-lose-exact-intersection-boundaries.md](../issues/M-distribution-mixed-support-operations-lose-exact-intersection-boundaries.md))
- [ ] Numerical integration (v0 Simpson's rule, trapezoidal)
- [ ] FFT transform (v0 `Distribution.fft()`)
- [ ] FWHM calculation (full width half maximum)
- [ ] Effective support (threshold-based range)
- [ ] Grid refinement (v0 `_adjust_grid()` adaptive)
- [ ] BranchLenInterpolator (v0: node-specific branch length probability)
  - [ ] Input mode (Gaussian/Poisson approximation)
  - [ ] Marginal mode (from profile pairs)
  - [ ] Joint mode (from compressed state pairs)
  - [ ] Gamma rescaling (for relaxed clock)
  - [ ] Adaptive grid construction (log near zero, quadratic tails)
