# Loose tolerances in test_gaussian_product.rs

`packages/treetime-distribution/src/__tests__/test_gaussian_product.rs` uses `assert_relative_eq!` with `epsilon = 0.1` to `1.0` at 10 call sites (lines 98, 115, 136, 151, 203, 216, 256, 269, 288, 292).

Likely justified by the grid-based numerical integration (201-point trapezoidal rule) used to compute Gaussian products, which has inherent discretization error scaling with grid spacing. The tolerances measure agreement between grid-based and analytical Gaussian product parameters (mean, sigma, FWHM).

Tightening requires either finer grids (slower tests) or analytical comparisons bypassing the grid.
