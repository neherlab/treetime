# Validate skyline grids and test the production objective

Make invalid skyline array shapes and domains unrepresentable before endpoint indexing.

## Required changes

- Add a validated skyline-grid type that requires nonempty equal-length time and log-$T_c$ arrays, finite values, and strictly increasing times.
- Construct `SkylineCostFunction` and $T_c$ distributions only from the validated type.
- Propagate contextual validation errors without panics or silent coercion.
- Replace formula-copy objective tests with calls to the production `CostFunction::cost` implementation.

## Validation

- Parameterized rejection tests for empty, unequal-length, non-finite, duplicate, and descending grids.
- Whole-value valid-grid construction tests.
- Production-objective comparisons against independently derived constant-$T_c$ cases.
- Full lint and test suite.

## Related issues

- Source: [kb/issues/N-coalescent-skyline-grid-validation-incomplete.md](../issues/N-coalescent-skyline-grid-validation-incomplete.md)
