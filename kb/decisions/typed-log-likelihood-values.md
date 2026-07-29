# Type aggregate log-likelihood values

Aggregate and component log-likelihoods cross module boundaries as `LogLh`, a transparent scalar newtype over `f64`. This makes log-domain values distinct from probabilities, linear likelihoods, branch lengths, dates, and optimizer coordinates at compile time.

This is a representation-safety decision, with no numerical or scientific divergence from TreeTime v0. The wrapped value remains an IEEE 754 `f64`, including support for negative infinity and positive log-likelihoods from continuous densities.

## Boundary

`LogLh` is used for scalar fields, parameters, and return values that represent complete or component log-likelihoods, including:

- dense, discrete, and sparse marginal likelihood totals;
- graph and partition root likelihoods;
- coalescent edge and tree totals;
- branch-optimization likelihood results and convergence histories;
- GTR refinement and mugration totals;
- timetree convergence components and totals.

Probability arrays, transition matrices, normalized profiles, derivative coefficients, local logarithms, optimizer coordinates, and distribution grids retain native `f64` and `ndarray` representations. Explicit `LogLh::new()` and `LogLh::value()` calls mark crossings between numerical kernels and typed aggregate boundaries.

## Arithmetic

The type supports only operations with domain-defined results:

- `LogLh + LogLh` and `Sum<LogLh>` combine independent log-likelihood components;
- `LogLh - LogLh` returns an `f64` difference used by convergence checks;
- `-LogLh` returns an `f64` optimizer cost.

It does not implement dereferencing, mixed arithmetic with `f64`, or implicit conversion from `f64`. Those operations would erase the distinction the type provides.

## Representation and serialization

`#[repr(transparent)]` preserves the size and alignment of `f64`. Serde's transparent representation preserves existing JSON and CSV scalar schemas: a `LogLh` serializes as a number rather than an object.

The implementation is defined in [`packages/treetime-primitives/src/log_lh.rs`](../../packages/treetime-primitives/src/log_lh.rs). Unit and compile-fail tests cover arithmetic, layout, serialization, and rejected mixed-domain expressions.
