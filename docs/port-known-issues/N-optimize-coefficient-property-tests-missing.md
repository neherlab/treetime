# Coefficient extraction lacks invariant property tests

The eigenvalue-space coefficient extraction (`get_coefficients()` in `optimize_dense.rs` and `optimize_sparse.rs`) computes `k_c = (msg_child . V) * (msg_parent . V_inv^T)`. These coefficients have mathematical invariants that are tested only by manual examples, not by property-based tests.

## Invariants not tested

- **Non-negative site likelihood at t=0**: `sum_c k_{ic} >= 0` for valid probability messages (coefficients sum to a probability, which is non-negative)
- **Multiplicity linearity**: sparse contribution with multiplicity `m` must produce `m *` the single-site log-likelihood, derivative, and second derivative
- **Dense-sparse equivalence**: `m` sparse sites with multiplicity 1 must match one dense contribution with `m` rows
- **Coefficient additivity**: coefficients from independent sites sum independently in the log-likelihood

## Why it matters

The Hessian multiplicity bug (sparse second derivative squaring multiplicity incorrectly) would have been caught by the multiplicity linearity property test. Property tests over random coefficients and branch lengths provide coverage that manual examples miss.

## Impact

Medium. Missing property tests leave algebraic invariants unverified. The `proptest` infrastructure already exists in this module (`test_coefficient_extraction_sparse/` has generators).

## Fix

Add `proptest` tests for the listed invariants, using the existing generator infrastructure.
