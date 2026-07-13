# Extract scalar coalescent merger-rate evaluation

Internal coalescent contribution formulas have scalar lineage-count and Tc values, but they allocate two one-element arrays and two result arrays through `compute_merger_rates()`. Skyline likelihood evaluation separately repeats the same clamp and total-rate arithmetic.

## Implementation

In [`packages/treetime/src/coalescent/integration.rs`](../../packages/treetime/src/coalescent/integration.rs), define a named generic result type instead of returning a same-type tuple:

```rust
pub struct MergerRates<T> {
  pub per_lineage: T,
  pub total: T,
}
```

Extract a reusable scalar calculation that implements the existing v0/v1 formula exactly:

```text
nlineages = max(0.5, k - 1)
per_lineage = nlineages / (2 * Tc)
total = nlineages * (nlineages + 1) / (2 * Tc)
```

Return `MergerRates<f64>` from the scalar function. Change the array-oriented `compute_merger_rates()` to return `MergerRates<Array1<f64>>` and build its values through the scalar calculation, keeping one source of truth for the clamp and arithmetic.

Use the scalar function in:

- [`compute_internal_contribution_single()`](../../packages/treetime/src/coalescent/contributions.rs), replacing the two one-element `Array1::from_vec` allocations and `lambda_t[0]`; and
- [`compute_total_neg_log_lh()`](../../packages/treetime/src/coalescent/skyline.rs), replacing its duplicate `nlineages` and `lambda` arithmetic.

Update current array callers and tests to use the named fields. Do not replace the array API: integration and golden-master tests still require vectorized results.

## Explicit non-goals

Do not change:

- calendar-time/TBP conversion or event-breakpoint sidedness;
- contribution signs, multiplicity, or multiplication order;
- Tc validation and non-positive-Tc behavior;
- Formula construction or per-node closure architecture;
- integral-merger-rate integration; or
- leaf/root coalescent contributions.

The scalar extraction must be behavior-preserving and must not absorb the unresolved coalescent refactor.

## Tests

Add parameterized scalar unit cases below, at, and above the clamp boundary `k = 1.5`, with representative positive Tc values. Assert exact results when the arithmetic is exactly representable and ULP equivalence otherwise.

For arrays containing the same cases, compare the entire `MergerRates<Array1<f64>>` fields against scalar-derived expected arrays. Include differing `k` and Tc values so an accidental broadcast or field swap fails.

Run the existing coalescent golden-master fixtures in [`test_gm_coalescent.rs`](../../packages/treetime/src/coalescent/__tests__/test_gm_coalescent.rs) and existing skyline tests unchanged in meaning. They must show that scalar reuse does not alter merger-rate arrays, node contributions, or skyline likelihoods beyond the project's `1e-6` parity threshold. Existing looser tolerances do not authorize a larger divergence. Validate the new tests by temporarily perturbing an expected scalar rate and confirming failure before restoring it.

## Acceptance criteria

- Scalar contribution evaluation performs no one-element ndarray allocation.
- Array and scalar paths share the exact clamp and rate formulas.
- Contributions and skyline use the scalar API.
- Existing golden fixtures remain unchanged and pass within the project's `1e-6` parity threshold; scalar unit cases use exact or ULP assertions as appropriate.
- No scientific, time-coordinate, validation, or Formula behavior changes.
- The full formatter, linter/build, and test suite pass:

  ```bash
  ./dev/docker/run ./dev/dev f
  ./dev/docker/run ./dev/dev l
  ./dev/docker/run ./dev/dev t
  ```

## Related issues

- Source: [kb/issues/N-coalescent-contribution-evaluation-allocation-overhead.md](../issues/N-coalescent-contribution-evaluation-allocation-overhead.md)
