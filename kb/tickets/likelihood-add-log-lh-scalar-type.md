# Add the LogLh scalar type

Introduce a transparent `LogLh` newtype and use it for aggregate log-likelihood values throughout the core TreeTime library. Preserve native `f64` values inside probability arrays and numerical kernels.

## Type contract

Add a public scalar type with this representation:

```rust
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
#[repr(transparent)]
#[serde(transparent)]
pub struct LogLh(f64);
```

Provide:

- `new(value: f64) -> LogLh`
- `value(self) -> f64`
- additive identity and impossible-likelihood constants where they improve existing call sites
- `Add<LogLh, Output = LogLh>` and `AddAssign<LogLh>`
- `Sum<LogLh, Output = LogLh>` and `Sum<&LogLh, Output = LogLh>`
- `Sub<LogLh, Output = f64>` for log-likelihood differences
- `Neg<Output = f64>` for numerical optimizer costs

Do not implement `Deref`, `DerefMut`, mixed arithmetic with `f64`, implicit `From<f64>`, or a finite/non-positive constructor restriction. `f64::NEG_INFINITY` is a valid impossible likelihood, and continuous densities can produce positive log-likelihoods.

Place the type in a core module that is dependency-safe for partition, coalescent, optimization, mugration, and timetree code. Export one canonical type; do not add aliases or compatibility wrappers because v1 has not shipped.

## Propagation scope

Use structural and semantic search to inventory every production field, parameter, and return value representing an aggregate or component log-likelihood. Propagate `LogLh` through all such boundaries, including:

- `DenseSeqDistribution::log_lh` and `SparseSeqDistribution::log_lh`
- `HasLogLh::get_log_lh()` implementations and `graph_log_lh()`
- dense, discrete, and sparse marginal inference totals
- coalescent edge and whole-tree log-likelihoods
- branch optimization likelihood result fields and objective boundaries
- GTR refinement likelihood totals
- mugration result log-likelihoods
- timetree convergence component and total fields

Keep local logarithms, normalized probability profiles, transition matrices, distribution grids, derivative coefficients, and optimizer coordinates as `f64` when they do not cross a log-likelihood boundary. Convert explicitly with `LogLh::new()` and `.value()` at those boundaries.

Preserve existing field and function names such as `log_lh`, `log_lh_total`, and `compute_sequence_log_lh()`. The type supplies compile-time semantics; names continue to communicate the domain in serialized schemas and source.

## Tests

- Unit-test `LogLh` construction, extraction, arithmetic, comparisons, constants, size, alignment, and transparent JSON serialization.
- Add compile-fail coverage showing that `LogLh + f64`, `f64 + LogLh`, and implicit raw-float assignment do not compile.
- Update existing tests to compare typed log-likelihoods without weakening floating-point tolerances or deriving expected values from the implementation under test.
- Verify unchanged JSON and convergence CSV scalar schemas.
- Run all project tests without filtering.

## Performance validation

Compare the parent revision and implementation revision using the project's performance workflow. Cover representative dense marginal, sparse marginal, and branch optimization workloads. Report benchmark distributions and profiles; investigate any statistically supported regression before completion.

The implementation must retain native ndarray element types and operations. Do not introduce per-element wrapping, validation branches in hot loops, allocation, dynamic dispatch, or conversion through `Vec`.

## Knowledge base updates

After full resolution:

- delete the source issue and its index entry
- update relevant feature and algorithm documentation to show typed aggregate log-likelihood boundaries
- move the selected proposal outcome to the appropriate implemented design documentation without adding a scientific divergence from v0

## Related issues

- Source: [kb/issues/M-likelihood-log-lh-values-use-raw-f64.md](../issues/M-likelihood-log-lh-values-use-raw-f64.md) -- delete after full resolution
- Proposal: [kb/proposals/typed-log-likelihood-values.md](../proposals/typed-log-likelihood-values.md)
