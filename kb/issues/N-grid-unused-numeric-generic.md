# Grid numeric generic has only one supported type

`InterpElem` has one implementation, `f64`, but its type parameter propagates through `Grid<T>`, `GridFn<T>`, and `DistributionFunction<T, Y>` [`packages/treetime-grid/src/lib.rs#L12`](../../packages/treetime-grid/src/lib.rs#L12) [`packages/treetime-grid/src/grid.rs#L24`](../../packages/treetime-grid/src/grid.rs#L24) [`packages/treetime-grid/src/grid_fn.rs#L31`](../../packages/treetime-grid/src/grid_fn.rs#L31) [`packages/treetime-distribution/src/distribution_core/function.rs#L12`](../../packages/treetime-distribution/src/distribution_core/function.rs#L12).

## Cost in the current API

- Numeric trait bounds repeat across grid construction, interpolation, serialization, and distribution composition.
- `NumCast` conversions introduce fallible-looking conversion points whose production instantiation is always `f64`.
- Generic serde bounds and error messages obscure the actual wire and numerical contract.
- Reviewers must reason about unsupported numeric types even though no second implementation can instantiate the abstraction.

`YAxisPolicy` is different: `Plain` and `NegLog` are both real production representations and must remain explicit.

## Decision axis

- O1. Make grid coordinates concrete `f64`, removing `InterpElem` and the unused type parameters. This states the actual numerical contract.
- O2. Retain generic coordinates and add at least one required non-`f64` implementation with tests defining conversion, interpolation, and serde semantics.

No second coordinate type or requirement is currently identified, so O1 is the supported direction. This is an API simplification only; interpolation formulas and floating-point behavior must remain unchanged.

## Validation

- Semantic impact analysis covers `treetime-grid` and `treetime-distribution` public types.
- Grid, interpolation, distribution, serde, and property tests produce identical `f64` results.
- No conversion unwrap remains solely to recover `f64` from the removed abstraction.
