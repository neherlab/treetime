# Make sparse-profile tolerance explicit

Remove variable tolerance plumbing from sparse-profile normalization assertions.

## Required changes

- Put an explicit literal `max_ulps` value at each scientific assertion site and remove the variable tolerance parameter.
- Preserve expected-first ordering and whole-profile comparison.
- See red by temporarily perturbing one expected profile, then restore it.

## Validation

- Run every sparse marginal normalization test and the full suite.
- Confirm the test fails under the injected fault before restoring it.

## Related issues

- Source: [kb/issues/N-test-floating-tolerance-hidden-by-variable.md](../issues/N-test-floating-tolerance-hidden-by-variable.md)
