# Uniform-spacing validation rejects valid offset grids

`has_uniform_spacing()` compares every adjacent spacing with the first using a fixed 100-ULP tolerance measured at the spacing magnitude. This rejects coordinate arrays produced from exact uniform-grid metadata when the coordinate offset is much larger than the spacing.

For a grid represented by

$$
x_i = x_{\min} + i\Delta x,
$$

floating-point subtraction of adjacent coordinates carries rounding error at the magnitude of $x_i$. Comparing the resulting differences in ULPs of the much smaller $\Delta x$ therefore magnifies cancellation error.

A representable grid with `x_min = 1813.11225188`, `dx = 0.19991935859615637`, and 1000 points regenerates adjacent differences between `0.19991935859593468` and `0.19991935859616206`. [`packages/treetime-utils/src/array/ndarray.rs`](../../packages/treetime-utils/src/array/ndarray.rs) rejects that array even though [`Grid::x_at()`](../../packages/treetime-grid/src/grid.rs) generated every coordinate from one stored spacing.

## Impact

Public array-based constructors can reject valid uniformly spaced data at large coordinate offsets. This affects calendar times, timestamps, and other domains where an offset is large relative to grid resolution. Callers that already hold a typed `Grid` can avoid the failure by preserving its metadata, but external arrays still need a sound validation contract.

## Required behavior

Define uniformity against a numerically justified error bound that accounts for coordinate magnitude and spacing while still rejecting genuinely non-uniform arrays. Cover large positive and negative offsets, small spacing, accumulated index multiplication, and deliberate spacing perturbations on both sides of the acceptance boundary.

This is distinct from [M-grid-invalid-numeric-domain.md](M-grid-invalid-numeric-domain.md), which covers non-finite values and spacing too small to produce strictly increasing represented coordinates.
