# GridFn::from_grid_array resets both boundaries to Error and drops the fitted law

`GridFn::from_grid_array()` builds a grid function with both tails set to `BoundaryBehavior::default()` (which is `Error`), discarding any policy and any fitted `SoftTailLaw` or `HardApproachLaw` the source carried [`packages/treetime-grid/src/grid_fn.rs#L49-L63`](../../packages/treetime-grid/src/grid_fn.rs#L49-L63). Every rebuild that goes through this constructor without re-applying the tails therefore loses the boundary declaration.

## Mechanism

`resample()` works around the reset by re-applying `self.left_extrap` and `self.right_extrap` after the rebuild [`packages/treetime-grid/src/grid_fn.rs#L513-L527`](../../packages/treetime-grid/src/grid_fn.rs#L513-L527). Callers that construct through `from_grid_array` directly do not. Under the log-space proposal the tail is not decoration: it is the domain declaration (hard versus soft) and the fitted approach or decay law that multiplication and re-windowing depend on. A path that regrids and loses the tail silently converts a declared soft or hard boundary into `Error`.

The backward child fold normalizes the accumulator after every multiply. If normalization or any `scale_y` path routes through `from_grid_array` without carrying the tail, the boundary law is discarded before the next child multiplies in, and the hard/soft multiplication rule never fires from the second child onward.

## Required behavior

- Carry the boundary policy and the fitted law across every regrid, not only through `resample()`.
- Make `from_grid_array` (or its callers `scale_y` and `normalize`) preserve `left_extrap` and `right_extrap`, so a regridded function keeps its declared domain and law.
- A round trip through regridding must leave a hard boundary hard, a soft boundary soft, and the fitted coefficients intact.

## Impact

This blocks the integrable `Linear` soft tail in the timetree passes ([N-timetree-passes-omit-integrable-linear-soft-tail.md](N-timetree-passes-omit-integrable-linear-soft-tail.md)): a fitted `Linear(Some(..))` tail attached before a normalize step is reset to `Error` and cannot survive to the next operation. It is the prerequisite for the soft-tail half of the log-space work.

## Related

- [kb/proposals/distribution-log-space-and-hard-soft-boundaries.md](../proposals/distribution-log-space-and-hard-soft-boundaries.md): Part B, "Policy must survive regridding".
