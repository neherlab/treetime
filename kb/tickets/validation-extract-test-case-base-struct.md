# Extract TestCaseBase struct for shared field accessors

Five TestCase implementations each define 5 identical accessor methods returning &self.field.

## Current state

`GaussianTestCase`, `ExponentialTestCase`, `GaussianExponentialTestCase`, `GaussianPairwiseMultiplicationTestCase`, and `GaussianChainMultiplicationTestCase` each implement `name()`, `description()`, `stress_type()`, `analytical_caution()`, `slowness()` identically.

## Target state

A `TestCaseBase` struct holds the common fields (`name`, `description`, `stress_type`, `analytical_caution`, `slowness`) with accessors implemented once. Each test case embeds `TestCaseBase` and forwards via delegation.

## Implementation

1. Define `TestCaseBase` struct with the 5 common fields
2. Implement the 5 accessor methods on `TestCaseBase`
3. Add a `base()` method to the `TestCase` trait returning `&TestCaseBase`
4. Provide default implementations for all 5 accessor methods via `self.base().field()`
5. Update all 5 test suites to embed `TestCaseBase` and implement `base()`
6. Remove the 25 explicit accessor implementations

## Related issues

Source: `kb/issues/N-validation-test-case-accessor-boilerplate.md`
