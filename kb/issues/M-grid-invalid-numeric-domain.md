# Grid constructors accept non-finite and non-representable spacing

The `Grid` invariant requires positive spacing, but constructors validate it with ordered comparisons that do not reject NaN. Rust floating-point values implement partial ordering, and `fn PartialOrd::partial_cmp()` returns `None` for a NaN comparison [[doc](https://doc.rust-lang.org/std/cmp/trait.PartialOrd.html#tymethod.partial_cmp)]:

- `fn Grid::from_start_dx()` accepts a NaN origin or spacing because `NaN <= 0` is false [packages/treetime-grid/src/grid.rs#L32-L43](../../packages/treetime-grid/src/grid.rs#L32-L43)
- `fn Grid::from_range_n_points()` accepts a NaN endpoint and can calculate NaN spacing [packages/treetime-grid/src/grid.rs#L46-L58](../../packages/treetime-grid/src/grid.rs#L46-L58)

`fn Grid::from_range_dx()` [packages/treetime-grid/src/grid.rs#L60-L76](../../packages/treetime-grid/src/grid.rs#L60-L76) and `fn GridFn::from_arrays_nonuniform()` [packages/treetime-grid/src/grid_fn.rs#L52-L86](../../packages/treetime-grid/src/grid_fn.rs#L52-L86) count their points through `MaxGridPoints::point_count()`, which returns an error for a NaN, negative, or infinite count and for a count above the grid point limit ([kb/decisions/distribution-grid-point-limit.md](../decisions/distribution-grid-point-limit.md)). They no longer panic on non-finite input, but they still accept a finite range whose spacing cannot be represented.

Finite inputs can also violate the effective-spacing invariant. At sufficiently large magnitudes, `x_min + dx` can round back to `x_min`, so adjacent generated coordinates are equal even though stored `dx` is positive. Interpolation and interval lookup then operate on a grid whose represented coordinates are not strictly increasing.

Derived `Deserialize` bypasses constructor-only validation and can materialize invalid fields directly.

## Uniform-grid acceptance silently changes coordinates

`fn Grid::from_array()` [packages/treetime-grid/src/grid.rs#L78-L94](../../packages/treetime-grid/src/grid.rs#L78-L94) uses `fn has_uniform_spacing()` [packages/treetime-utils/src/array/ndarray.rs#L141-L153](../../packages/treetime-utils/src/array/ndarray.rs#L141-L153) whose tolerance `MAX_SPACING_ULPS * endpoint_magnitude * epsilon` is independent of the nominal spacing. This means:

- Arrays with materially unequal intervals pass the uniformity check when endpoint magnitudes are large relative to spacing.
- On acceptance, `from_array()` discards all interior coordinates and reconstructs the grid from `x[0]` and `x[1] - x[0]` [packages/treetime-grid/src/grid.rs#L88-L93](../../packages/treetime-grid/src/grid.rs#L88-L93), producing coordinates that differ from the input array. The rounding error of the first difference grows with the index: the last coordinate `x_min + dx * (n - 1)` [packages/treetime-grid/src/grid.rs#L100-L105](../../packages/treetime-grid/src/grid.rs#L100-L105) can differ from `x[n - 1]` by up to about $(n - 1)\,\lvert x \rvert\,2^{-52}$, about $10^{-10}$ for 300 points near the year 2000. A drift of about $4 \cdot 10^{-11}$ at the grid end was measured during work on the mass model of distributions (2026-10-18). Taking the spacing from the end points, `(x[n - 1] - x[0]) / (n - 1)`, keeps the last coordinate equal to the input
- Non-finite, zero-spacing, and descending grids can also pass because the tolerance check uses ordered comparisons that return `false` for NaN.

## Required behavior

Require finite endpoints and finite positive spacing, use checked point-count conversion, and verify that generated adjacent coordinates are strictly increasing in the represented floating-point type. Validate that the uniform-spacing tolerance scales with the nominal spacing, not only with endpoint magnitude. Apply the invariant at constructors, deserialization, and every public grid-producing boundary. Failures must return contextual errors rather than panic. Count conversions already cannot fail: `fn from_count()` converts exactly and `fn Grid::find_interval_index()` falls back to the last interval [packages/treetime-grid/src/grid.rs#L126-L139](../../packages/treetime-grid/src/grid.rs#L126-L139).

## Fix

- Reject NaN and positive or negative infinity for origins, endpoints, and spacing
- Reject zero and negative spacing
- Reject ranges whose generated adjacent coordinates are not strictly increasing in the represented type
- Take the spacing of `fn Grid::from_array()` from the end points of the array
- Validate derived deserialization, or replace it with a validated implementation, so serialized input cannot bypass the invariant

## Validation

Cover every constructor, deserialization, and the `GridFn` boundary with these cases:

- NaN and positive or negative infinity in each input
- Reversed range and zero spacing
- Extreme finite span
- Repeated representable coordinates (large magnitude where `x_min + dx` rounds to `x_min`)
