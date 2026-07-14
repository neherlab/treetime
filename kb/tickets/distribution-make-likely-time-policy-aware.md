# Make distribution likely-time selection policy-aware

`Distribution<Y>::likely_time()` selects the maximum stored ordinate for `Function` and `Formula` variants. That is correct for `Plain` and inverted for `NegLog`, where the minimum ordinate represents the highest likelihood.

## Implementation

- Make likely-time selection use the active `YAxisPolicy` semantics for both sampled functions and formulas.
- Keep the public result and the `Empty`, `Point`, and `Range` contracts unchanged.
- Reuse one policy-aware extremum rule across `DistributionFunction` and `DistributionFormula`; do not duplicate plain/neg-log branching in callers.

## Validation

- Add parameterized unit cases for `Function` and `Formula` variants.
- Construct matching `DistributionPlain` and `DistributionNegLog` values whose peak and trough occur at the same time, then assert that both return that complete expected value.
- Include an asymmetric grid or formula so selecting the opposite extremum cannot pass accidentally.
- Confirm each new test fails when the extremum direction is inverted, then restore the implementation and run the complete test suite.

## Related issues

- Source: [kb/issues/M-distribution-neglog-likely-time-selects-maximum.md](../issues/M-distribution-neglog-likely-time-selects-maximum.md)
