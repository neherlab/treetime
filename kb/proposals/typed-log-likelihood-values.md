# Type log-likelihood values at scalar boundaries

TreeTime represents aggregate log-likelihoods as `f64`. Names such as `log_lh` communicate the representation to readers, but the compiler cannot distinguish a log-likelihood from a probability, likelihood, branch length, calendar date, or optimization parameter. A misplaced raw value can therefore compile and produce plausible numerical output under the wrong algebra.

The problem resembles units of measurement: the machine representation is identical while valid operations depend on the quantity. Rust's newtype idiom provides compile-time distinction, and `#[repr(transparent)]` preserves the wrapped scalar's layout and ABI. TreeTime already uses this pattern for `CalendarTime` and uses phantom policy types to distinguish `Distribution<Plain>` from `Distribution<NegLog>`.

## Constraints

- Scientific scalar types must prevent accidental domain mixing at module boundaries.
- Numerical kernels must retain native `f64` and `ndarray` operations.
- The representation must not allocate or introduce dynamic dispatch.
- Serialization must remain a JSON and CSV number rather than a wrapper object.
- Valid log-likelihoods include `f64::NEG_INFINITY`; continuous densities can also yield positive log-likelihoods. A global finite or non-positive invariant would reject valid domain values.
- The API must remain convenient for summing independent components and computing likelihood differences.

## Ecosystem patterns

- Rust documents the [newtype idiom](https://doc.rust-lang.org/rust-by-example/generics/new_types.html) as a compile-time guarantee that the correct kind of value is supplied.
- Rust's [`transparent` representation](https://doc.rust-lang.org/reference/type-layout.html#the-transparent-representation) gives a one-field wrapper the layout and ABI of its non-zero-sized field.
- [`bio::stats::LogProb`](https://docs.rs/bio/latest/bio/stats/probs/struct.LogProb.html) demonstrates a dedicated log-probability type and stable log-space operations. It is not suitable as TreeTime's general log-likelihood type because log-probabilities and log-likelihoods have different domains, its tuple field is public, and its probability conversion uses Rust-Bio's approximate `FastExp` implementation.
- [`uom`](https://docs.rs/uom/latest/uom/) provides generic dimensional analysis. Log representation changes algebra rather than physical dimension, so a full unit system does not match this problem.
- [`nutype`](https://docs.rs/nutype/latest/nutype/) and [`typed_floats`](https://docs.rs/typed_floats/latest/typed_floats/) enforce numeric invariants. They do not encode the distinction between probability, likelihood, and their logarithms.

## Design axes

### Public type vocabulary

1. **`LogLh` (selected)**: short, explicit, and consistent with established `log_lh` identifiers.
2. `LogLikelihood`: maximally descriptive but repetitive in numerical code.
3. `LogValue<Likelihood>`: reusable machinery at the cost of generic parameters and less direct diagnostics.

Public vocabulary remains domain-specific. A general `LogValue` or `TypedFloat` would permit semantically meaningless tags and obscure which algebra applies.

### Typing boundary

1. **Aggregate scalar boundaries (selected)**: use `LogLh` for fields, return values, and inter-module parameters that represent complete or component log-likelihoods.
2. Every probability element: replace `f64` inside profiles and matrices.
3. Output boundaries only: type serialized metrics while retaining raw internal APIs.

Aggregate scalar typing catches cross-domain mistakes without disrupting `ndarray`, matrix multiplication, interpolation, or vectorized profile operations. Whole-array domain structs such as `DenseSeqDistribution` continue to hold native `Array2<f64>` values.

### Arithmetic

1. **Domain-directed operators (selected)**:
   - `LogLh + LogLh -> LogLh` combines independent log-likelihood components.
   - `Sum<LogLh> -> LogLh` uses `0.0` as the additive identity.
   - `LogLh - LogLh -> f64` produces a log-likelihood difference rather than claiming the result is itself a likelihood.
   - `-LogLh -> f64` exposes an optimizer cost at the numerical optimizer boundary.
2. No operators: require named methods for every operation.
3. Transparent `f64` arithmetic: implement `Deref<Target = f64>` and mixed scalar operators.

Mixed `LogLh`/`f64` arithmetic and `Deref` would recreate the ambiguity the type is intended to prevent. Explicit `new()` and `value()` methods mark representation boundaries.

### Validation

1. **Semantic type without a global range restriction (selected)**: preserve all IEEE values accepted by current algorithms and validate finite/range requirements where an algorithm requires them.
2. Reject non-finite values in the constructor.
3. Restrict values to `[-inf, 0]` as log-probabilities.

The selected design separates representation safety from algorithm-specific validity. It does not normalize, exponentiate, approximate, clamp, or otherwise change numerical behavior.

## Selected design

```rust
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
#[repr(transparent)]
#[serde(transparent)]
pub struct LogLh(f64);
```

`LogLh` provides `new()`, `value()`, domain-directed arithmetic, and constants for the additive identity and impossible likelihood where they improve call-site clarity. It does not implement `Deref`, `DerefMut`, mixed arithmetic with `f64`, or implicit conversion from `f64`.

The type applies to aggregate likelihood values such as partition root log-likelihoods, graph totals, coalescent totals, optimizer likelihood results, mugration results, and convergence metrics. Local intermediate logarithms inside numerical kernels remain `f64` until they cross a typed boundary.

## Performance contract

- `LogLh` has the size and alignment of `f64`.
- Release code contains no validation branch, allocation, dynamic dispatch, or virtual call introduced by wrapping and unwrapping.
- Representative dense marginal, sparse marginal, and branch optimization benchmarks show no statistically supported runtime regression against the parent revision.
- Probability profiles and transition matrices retain their native ndarray element types and operations.

## Validation

- Unit tests verify construction, extraction, addition, summation, subtraction, negation, comparison, and transparent Serde behavior.
- Compile-fail coverage verifies that raw `f64` and `LogLh` cannot be mixed implicitly.
- Existing dense, sparse, coalescent, optimizer, mugration, timetree, JSON, and CSV tests verify unchanged numerical and serialized output.
- Full lint, format, build, and test validation passes.
- Before/after performance measurements cover representative inference and optimization workloads.

## Actionable issue

- [kb/issues/M-likelihood-log-lh-values-use-raw-f64.md](../issues/M-likelihood-log-lh-values-use-raw-f64.md)
