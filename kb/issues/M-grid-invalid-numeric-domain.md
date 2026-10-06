# Grid constructors accept non-finite and non-representable spacing

The `Grid` invariant requires positive spacing, but constructors validate it with ordered comparisons that do not reject NaN. Rust floating-point values implement partial ordering, and `fn PartialOrd::partial_cmp()` returns `None` for a NaN comparison [[doc](https://doc.rust-lang.org/std/cmp/trait.PartialOrd.html#tymethod.partial_cmp)]:

- `from_start_dx()` accepts a NaN origin or spacing because `NaN <= 0` is false. [`packages/treetime-grid/src/grid.rs#L33-L44`](../../packages/treetime-grid/src/grid.rs#L33-L44)
- `from_range_n_points()` accepts a NaN endpoint and can calculate NaN spacing. [`packages/treetime-grid/src/grid.rs#L50-L62`](../../packages/treetime-grid/src/grid.rs#L50-L62)

`from_range_dx()` [`packages/treetime-grid/src/grid.rs#L64-L80`](../../packages/treetime-grid/src/grid.rs#L64-L80) and `GridFn::from_arrays_nonuniform()` [`packages/treetime-grid/src/grid_fn.rs#L62-L96`](../../packages/treetime-grid/src/grid_fn.rs#L62-L96) count their points through `MaxGridPoints::point_count()`, which returns an error for a NaN, negative, or infinite count and for a count above the grid point limit ([kb/decisions/distribution-grid-point-limit.md](../decisions/distribution-grid-point-limit.md)). They no longer panic on non-finite input, but they still accept a finite range whose spacing cannot be represented.

Finite inputs can also violate the effective-spacing invariant. At sufficiently large magnitudes, `x_min + dx` can round back to `x_min`, so adjacent generated coordinates are equal even though stored `dx` is positive. Interpolation and interval lookup then operate on a grid whose represented coordinates are not strictly increasing.

Derived `Deserialize` bypasses constructor-only validation and can materialize invalid fields directly.

## Uniform-grid acceptance silently changes coordinates

`fn Grid::from_array()` [`packages/treetime-grid/src/grid.rs#L82-L98`](../../packages/treetime-grid/src/grid.rs#L82-L98) uses `fn has_uniform_spacing()` [`packages/treetime-utils/src/array/ndarray.rs#L267-L279`](../../packages/treetime-utils/src/array/ndarray.rs#L267-L279) whose tolerance `MAX_SPACING_ULPS * endpoint_magnitude * epsilon` is independent of the nominal spacing. This means:

- Arrays with materially unequal intervals pass the uniformity check when endpoint magnitudes are large relative to spacing.
- On acceptance, `from_array()` discards all interior coordinates and reconstructs the grid from `x[0]` and `x[1] - x[0]`, producing coordinates that differ from the input array.
- Non-finite, zero-spacing, and descending grids can also pass because the tolerance check uses ordered comparisons that return `false` for NaN.

## Required behavior

Require finite endpoints and finite positive spacing, use checked point-count conversion, and verify that generated adjacent coordinates are strictly increasing in the represented floating-point type. Validate that the uniform-spacing tolerance scales with the nominal spacing, not only with endpoint magnitude. Apply the invariant at constructors, deserialization, and every public grid-producing boundary. Failures must return contextual errors rather than panic.

## Fix

- Reject NaN and positive or negative infinity for origins, endpoints, and spacing
- Reject zero and negative spacing
- Reject ranges whose generated adjacent coordinates are not strictly increasing in the represented type
- Replace the fallible numeric-conversion `unwrap()` calls with contextual errors: the `T::from()` conversions [`packages/treetime-grid/src/grid.rs#L60`](../../packages/treetime-grid/src/grid.rs#L60), [`#L112`](../../packages/treetime-grid/src/grid.rs#L112), [`#L135`](../../packages/treetime-grid/src/grid.rs#L135), [`packages/treetime-grid/src/grid_fn.rs#L91`](../../packages/treetime-grid/src/grid_fn.rs#L91), and the float-to-`usize` interval index [`packages/treetime-grid/src/grid.rs#L150`](../../packages/treetime-grid/src/grid.rs#L150)
- Validate derived deserialization, or replace it with a validated implementation, so serialized input cannot bypass the invariant

## Validation

Cover every constructor, deserialization, and the `GridFn` boundary with these cases:

- NaN and positive or negative infinity in each input
- Reversed range and zero spacing
- Extreme finite span
- Repeated representable coordinates (large magnitude where `x_min + dx` rounds to `x_min`)
