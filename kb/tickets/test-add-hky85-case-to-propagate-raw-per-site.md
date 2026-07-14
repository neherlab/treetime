# Add HKY85 test case to propagate_raw_per_site tests

## Problem

The `propagate_raw_per_site` forward and backward tests both use JC69, which has a symmetric transition matrix ($P(t) = P(t)^T$). This means forward and backward tests produce identical expected values and cannot distinguish whether the `transpose` flag is handled correctly.

## Task

Add an HKY85 (or other non-symmetric model) test case where forward and backward produce different results, verifying the transpose code path is exercised.

## Location

`packages/treetime/src/partition/marginal_helpers.rs:209,248:`

## Notes

- HKY85 with unequal base frequencies (e.g. pi = [0.4, 0.1, 0.1, 0.4]) produces an asymmetric P(t)
- Expected values must be derived from HKY85 closed-form or a verified reference implementation, not from `gtr.expQt_with_rate()` (that would reintroduce the circular oracle)
- The v0 Python GTR can serve as oracle for non-JC69 models

## Related issues

- Source: [kb/issues/N-test-quality-deficiencies.md](../issues/N-test-quality-deficiencies.md)
