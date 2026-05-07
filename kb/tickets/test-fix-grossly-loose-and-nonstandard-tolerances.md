# Fix grossly loose, non-standard, and variable test tolerances

## Grossly loose (hard failure per project rules)

- `commands/timetree/coalescent/__tests__/test_integration.rs:157:` `epsilon = 10.0`
- `commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs:80:` `epsilon = 1e0` (= 1.0). Comment: "Tolerance increased from 0.9 to 1.0"

## Non-standard format

- `commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_poisson.rs:46:` `epsilon = 3e-1`
- `commands/mugration/__tests__/test_gm_mugration.rs:92:` `epsilon = 2e-2`

## Variable tolerances (7 instances)

- `commands/timetree/inference/__tests__/test_branch_length_likelihood.rs:80,129,162:` `GRID_SPACING_BL`, `expected_epsilon` (computed)
- `commands/timetree/output/__tests__/test_confidence_extract.rs:232:` `epsilon = dx` (computed)
- `commands/timetree/coalescent/__tests__/test_gm_coalescent.rs:102:` `worst_err < PASS_THRESHOLD` (named constant)
- `commands/ancestral/__tests__/test_marginal_dense.rs:364:` `let epsilon = 1e-6;` (variable)
- `gtr/__tests__/test_prop_gtr_site_specific.rs:283:` `1e-2` with `TODO(investigate)` tag. Justification: linear interpolation on 61-point grid has inherent accuracy limits; observed max error ~1.8e-3.

## Sparse root-invariance proptest tolerance 1e-1

`packages/treetime/src/commands/ancestral/__tests__/test_marginal_root_invariance_prop.rs:57:`

Hides >2-orders-of-magnitude pulley-principle violation. Dense uses 1e-6. The 5-order gap indicates real algorithmic divergence. Related: `M-ancestral-sparse-root-invariance.md`.

## Loose tolerance 1e-6 in div tests

`packages/treetime/src/seq/__tests__/test_div.rs`

6 assertions use `epsilon = 1e-6`, 1 uses `epsilon = 1e-9`. Whether 1e-6 is the tightest passing tolerance has not been measured.

## Related issues

- Source: [N-test-quality-deficiencies.md](../issues/N-test-quality-deficiencies.md) -- delete after full resolution
