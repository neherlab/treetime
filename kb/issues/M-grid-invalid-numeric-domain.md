# Grid constructors accept non-finite and non-representable spacing

The `Grid` invariant requires positive spacing, but constructors validate it with ordered comparisons that do not reject NaN. Rust floating-point values implement partial ordering, and `fn PartialOrd::partial_cmp()` returns `None` for a NaN comparison [[doc](https://doc.rust-lang.org/std/cmp/trait.PartialOrd.html#tymethod.partial_cmp)]:

- `from_start_dx()` accepts a NaN origin or spacing because `NaN <= 0` is false. [`packages/treetime-grid/src/grid.rs#L31-L42`](../../packages/treetime-grid/src/grid.rs#L31-L42)
- `from_range_n_points()` accepts a NaN endpoint and can calculate NaN spacing. [`packages/treetime-grid/src/grid.rs#L44-L56`](../../packages/treetime-grid/src/grid.rs#L44-L56)
- `from_range_dx()` accepts non-finite inputs and unwraps a fallible float-to-`usize` conversion. [`packages/treetime-grid/src/grid.rs#L58-L73`](../../packages/treetime-grid/src/grid.rs#L58-L73)

Finite inputs can also violate the effective-spacing invariant. At sufficiently large magnitudes, `x_min + dx` can round back to `x_min`, so adjacent generated coordinates are equal even though stored `dx` is positive. Interpolation and interval lookup then operate on a grid whose represented coordinates are not strictly increasing.

Two public boundaries bypass constructor-only validation: derived `Deserialize` can materialize invalid fields directly, and `fn GridFn::from_arrays_nonuniform()` [`packages/treetime-grid/src/grid_fn.rs#L88-L118`](../../packages/treetime-grid/src/grid_fn.rs#L88-L118) performs a float-to-`usize` conversion before delegating to `Grid`.

## Required behavior

Require finite endpoints and finite positive spacing, use checked point-count conversion, and verify that generated adjacent coordinates are strictly increasing in the represented floating-point type. Apply the invariant at constructors, deserialization, and every public grid-producing boundary. Failures must return contextual errors rather than panic.

## Related tickets

- [kb/tickets/grid-validate-finite-representable-spacing.md](../tickets/grid-validate-finite-representable-spacing.md)
