# Coalescent time-scale coordinate is not type-enforced

The coalescent API represents $T_c$ with the generic `Distribution` type. Constant values are coordinate-independent, while nonconstant values require an external convention for their horizontal axis. The coalescent contribution refactor interprets that axis as decimal calendar years, but the type system cannot prevent a caller from supplying a time-before-present distribution.

## Design axis: time-scale representation

- O1. Keep `Distribution` and document the calendar-year contract at every coalescent boundary.
- O2. Introduce a dedicated coalescent time-scale type whose variants encode constant and calendar-domain piecewise values, making invalid coordinate states unrepresentable.

The coalescent API retains O1: nonconstant values use decimal calendar years. O2 remains an open alternative and has no implementation ticket until its representation and serialization contracts are decided.

## Related issues

- [kb/algo/coalescent-contribution-refactor.md](../algo/coalescent-contribution-refactor.md)
- [N-coalescent-skyline-grid-validation-incomplete.md](N-coalescent-skyline-grid-validation-incomplete.md)
