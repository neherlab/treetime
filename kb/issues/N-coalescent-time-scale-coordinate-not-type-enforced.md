# Coalescent time-scale coordinate is not type-enforced

> [!IMPORTANT]
> **Decision required.** The API keeps the documented calendar-year contract on `Distribution` (O1). A dedicated coalescent time-scale type (O2) remains an open alternative, and its representation and serialization contracts are undecided. The options are under "Design axis: time-scale representation".

The coalescent API represents $T_c$ with the generic `Distribution` type. Constant values are coordinate-independent, while nonconstant values require an external convention for their horizontal axis. The coalescent contribution refactor interprets that axis as decimal calendar years, but the type system cannot prevent a caller from supplying a time-before-present distribution.

## Design axis: time-scale representation

- O1. Keep `Distribution` and document the calendar-year contract at every coalescent boundary.
- O2. Introduce a dedicated coalescent time-scale type whose variants encode constant and calendar-domain piecewise values, making invalid coordinate states unrepresentable.

The coalescent API retains O1: nonconstant values use decimal calendar years. O2 remains an open alternative, and its implementation waits until its representation and serialization contracts are decided.

## Related issues

- [kb/algo/coalescent-contribution-refactor.md](../algo/coalescent-contribution-refactor.md)

## Proposed implementation plan

The problem is current: `CoalescentPrior::new` (`packages/treetime/src/coalescent/mean_field.rs` L42) and `LineageExpectation::merger_rate_schedule` (L469) take $T_c$ as a plain `&Distribution`. The plan depends on the open decision.

- **O1**: name the calendar-year axis in the parameter types and the schema docs of every coalescent boundary that takes $T_c$; no code change of behavior; then delete this issue
- **O2**: introduce a type such as `enum CoalescentTimescale { Constant(f64), Calendar(PiecewiseConstantFn) }` in `packages/treetime/src/coalescent/`, built only from the skyline fit or a constant, and replace `&Distribution` in `CoalescentPrior::new`, `merger_rate_schedule`, and the timescale field of the lineages; decide its serialization in the coalescent output first. Tests: construction from a skyline fit evaluates to the fitted segment values at calendar times (oracle: the segment values), and a constant evaluates to itself at any time
- Validation: `just check`, `just test-rs`, and `just generated-check` when a schema changes
