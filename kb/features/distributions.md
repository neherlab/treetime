# Probability Distributions

## Implemented

- [x] Analytical Gaussian (PDF, product, convolution)
- [x] Analytical Exponential (PDF, convolution, special case a=b)
- [x] Gaussian-Exponential convolution
- [x] Grid distributions (uniform grid, linear interpolation, and explicit boundary behavior; [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md))
- [x] ScaledArray pattern (normalized values + log-scale factor)

## Partial or Not Implemented

- [ ] Unified Distribution class (v0: wraps scipy.interpolate.interp1d)
- [ ] Delta functions (point masses, v0 `Distribution.delta_function()`)
- [x] Distribution multiplication (uniform-grid divergence documented in [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md))
- [/] Distribution division (implemented; [kb/issues/M-distribution-plain-division-fixed-floor.md](../issues/M-distribution-plain-division-fixed-floor.md))
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
