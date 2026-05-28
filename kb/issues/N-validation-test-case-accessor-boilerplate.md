# TestCase trait field accessors duplicated across five test suites

Five `TestCase` implementations each define 5 accessor methods that return `&self.field` identically (25 identical method bodies total).

## Locations

Test suite files in `packages/treetime-validation/src/testing/test_suites/`:

- `gaussian.rs:75:` `GaussianTestCase`
- `exponential.rs:72:` `ExponentialTestCase`
- `gaussian_exponential.rs:71:` `GaussianExponentialTestCase`
- `gaussian_pairwise_multiplication.rs:100:` `GaussianPairwiseMultiplicationTestCase`
- `gaussian_chain_multiplication.rs:76:` `GaussianChainMultiplicationTestCase`

Duplicated accessors: `name()`, `description()`, `stress_type()`, `analytical_caution()`, `slowness()` - all return `&self.field`.

## Impact

Pure maintainability. Adding a test suite requires copying 5 trivial accessors.

## Action

Extract a `TestCaseBase` struct holding the common fields and implement the accessors once. Each test suite embeds `TestCaseBase` and delegates via a `Deref` impl or explicit forwarding. Alternative: a derive macro for the trait.
