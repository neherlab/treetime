# GTR site-specific interpolation tolerance requires investigation

## Summary

Property tests for GTR site-specific rate variation use 1e-2 tolerance with a `TODO(investigate)` tag. The tolerance is 100x looser than the project standard and requires verification of whether it reflects genuine interpolation error or test gaming.

## Details

`packages/treetime/src/gtr/__tests__/test_prop_gtr_site_specific.rs:283:`

The test compares GTR transition probabilities computed via eigendecomposition against values interpolated from a precomputed 61-point grid. The justification comment (lines 244-270) explains:

- Linear interpolation on a 61-point grid has inherent accuracy limits
- Observed max error is approximately 1.8e-3
- The 1e-2 tolerance provides margin above observed error

## Open questions

1. Was the 1e-2 tolerance measured as the interpolation error bound, or was it widened iteratively until the test passed?
2. Is the 61-point grid density sufficient for the rate variation range being tested?
3. Would a denser grid (e.g., 200 points) reduce interpolation error below 1e-6?
4. Does the error scale with rate magnitude or branch length?

## Required investigation

- Measure actual max interpolation error across the full parameter space (not just observed in random samples)
- Determine whether the error bound is a property of the grid density or of the interpolation method
- If grid-density-limited: compute the grid size needed for 1e-6 accuracy and assess performance tradeoff
- If method-limited: evaluate cubic spline or other higher-order interpolation

## Impact

- The test validates that site-specific GTR is approximately correct but does not enforce tight numerical agreement
- If the tolerance masks a real bug, site-specific rate variation results could be wrong by up to 1%
