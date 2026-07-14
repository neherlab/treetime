# Return errors for unsupported Formula distribution operations

Replace the Formula panic arms in division, mapping, convolution, and scalar multiplication with actionable `Err` results that identify the operation and unsupported operand combination.

## Acceptance criteria

- Exercise `Formula` against every other distribution variant in each operand position supported by the operation's arity.
- Assert that public operation entry points do not unwind for any distribution-variant combination.
- Preserve existing implemented behavior for non-Formula variants.
- Do not add an analytic Formula implementation without a separately approved numerical contract.

## Related issues

Source: [kb/issues/H-distribution-result-api-panics-on-formula.md](../issues/H-distribution-result-api-panics-on-formula.md)
