# Grid resampling collects its values into a vector without its known length

`fn GridFn::resample_with_extrap()` [packages/treetime-grid/src/grid_fn.rs#L388-L402](../../packages/treetime-grid/src/grid_fn.rs#L388-L402) and `fn GridFn::resample_range_dx_clamped()` [packages/treetime-grid/src/grid_fn.rs#L404-L415](../../packages/treetime-grid/src/grid_fn.rs#L404-L415) map the indices of the new grid through a fallible interpolation and collect with `collect::<Result<Vec<T>, Report>>()`, then convert with `Array1::from_vec`. The point count is known before the loop. Collecting into a `Result` passes a lower size hint of zero to the vector, so the vector can grow by repeated reallocation and hold up to twice its final capacity for a moment. The project rules also ask for `ndarray` constructors over `Array1::from(vec)`.

Resampling runs on every regrid of the time inference, with grids of up to millions of points.

> [!IMPORTANT]
> **Investigation required.** Whether these collects reallocate in the release build, and their share of the run time and of the peak memory on `data/sc2/5100` and `data/mpox/clade-ii/1000`.

## Fix direction

Allocate the result once at its known length, for example with `Array1::from_shape_fn` and a fallible fill, or `Array1::zeros(n)` filled through `Zip` with an early error.

## Validation

The resampling tests of `treetime-grid` pass; a profile of `data/sc2/5100` shows no reallocation in these functions.
