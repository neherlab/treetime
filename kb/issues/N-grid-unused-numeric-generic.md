# Grid numeric generic has only one supported type

`InterpElem` has one implementation, `f64`, but its type parameter propagates through `Grid<T>`, `GridFn<T>`, and `DistributionFunction<T, Y>` [`packages/treetime-grid/src/lib.rs#L16`](../../packages/treetime-grid/src/lib.rs#L16) [`packages/treetime-grid/src/grid.rs#L26`](../../packages/treetime-grid/src/grid.rs#L26) [`packages/treetime-grid/src/grid_fn.rs#L20`](../../packages/treetime-grid/src/grid_fn.rs#L20) [`packages/treetime-distribution/src/distribution_core/function.rs#L13`](../../packages/treetime-distribution/src/distribution_core/function.rs#L13).

## Cost in the current API

- Numeric trait bounds repeat across grid construction, interpolation, serialization, and distribution composition.
- `InterpElem` requires `From<f64> + Into<f64>` so that grid code converts to and from `f64` without fallible `NumCast` calls; only `f64` can meet this bound in practice.
- Generic serialization bounds and error messages obscure the actual wire and numerical contract.
- Reviewers must reason about unsupported numeric types even though no second implementation can instantiate the abstraction.

`YAxisPolicy` is different: `Plain` and `NegLog` are both real production representations and must remain explicit.

## Decision axis

- O1. Make grid coordinates concrete `f64`, removing `InterpElem` and the unused type parameters. This states the actual numerical contract.
- O2. Retain generic coordinates and add at least one required non-`f64` implementation with tests defining conversion, interpolation, and serialization semantics.

No second coordinate type or requirement is currently identified, so O1 is the supported direction. This is an API simplification only; interpolation formulas and floating-point behavior must remain unchanged.

## Validation

- Semantic impact analysis covers `treetime-grid` and `treetime-distribution` public types.
- Grid, interpolation, distribution, serialization, and property tests produce identical `f64` results.
- The `From<f64> + Into<f64>` bounds of `InterpElem` disappear with the type parameter.
