# Per-edge optimization: test rewrite for 6 methods

## Problem

The current test suite targets 3 methods. After the preceding method selection issues are implemented, 6 methods exist but only 3 are tested. The golden master uses `BranchOptMethod::Newton`, which has a known Hessian dominance limitation -- it should use `BrentSqrt` (matching v0) for meaningful comparison against v0 reference outputs.

28 call sites across 12 test files hardcode `BranchOptMethod::Newton`. All code paths need coverage.

## Current state

- `test_optimize_method.rs`: 12 test functions / 40 cases targeting 3 methods
- `test_gm_optimize.rs`: golden master uses `BranchOptMethod::Newton` (line 216)
- 10 other test files pass a single `BranchOptMethod` variant to `run_optimize_mixed`
- Smoke tests (`dev/run-smoke-tests`): no `--opt-method` entries

## Approach

### Rewrite `test_optimize_method.rs`

The existing chain rule sqrt tests are kept. New chain rule log tests (from [M-optimize-method-newton-log](M-optimize-method-newton-log.md)) are added. The method-level tests are rewritten for 6 methods.

**Verification criteria** (each is a separate test or test group):

C1 -- Local optimality: the combined (substitution + indel) log-likelihood at the reported optimum exceeds the log-likelihood at nearby points. Use 3 delta scales ($\pm 0.1\%$, $\pm 1\%$, $\pm 10\%$ of optimal $t$). Parameterize across all 6 methods.

C2 -- Stationarity (Newton-specific): at the reported optimum, the implied Newton step $|\ell'/\ell''|$ is small relative to the Newton tolerance. A 10x factor is acceptable (the step from the final position can be larger than the convergence threshold due to nonlinearity). Parameterize across the 3 Newton variants.

C3 -- Cross-method log-likelihood agreement: run all 6 methods on the same input, verify all produce combined log-likelihood within tolerance of each other. Compare log-likelihood not branch lengths (the objective can be flat near the optimum). Tolerance around `1e-3` in log-likelihood.

C4 -- Bracket validity (Brent-specific): the log-likelihood at the optimum exceeds the log-likelihood at both bracket endpoints. Parameterize across all 3 Brent variants.

C5 -- Cross-conditioning (Newton ordering): on indel-bearing edges, `lh_newton_log >= lh_newton_sqrt >= lh_newton` within tolerance. Better-conditioned parameterizations find equal or better optima. This documents the known Newton-t limitation: it is a correct implementation of a limited algorithm, not a bug. If Newton-t produces better log-likelihood than Newton-log, the implementation is wrong.

**Additional tests:**

- Positivity and finiteness with indels ($k = 1, 2, 4$), parameterized across all 6 methods
- Positivity and finiteness without indels, parameterized across all 6
- Brent variants agree with best Newton variant within tolerance

The existing test helpers (`setup_with_indels`, `eval_combined_first_edge`, `eval_metrics_first_edge`) are reusable. They capture the indel rate before optimization (same rate the optimizer uses) for consistent evaluation.

### Parameterize 10 other test files

Each test function that calls `run_optimize_mixed` with a hardcoded method gains an `#[rstest]` parameterization across all 6 methods.

Files to parameterize (call count in parentheses):

| Test file                                          | Calls |
| :------------------------------------------------- | :---- |
| `test_convergence/test_convergence_edge_cases.rs`  | 3     |
| `test_convergence/test_convergence_idempotence.rs` | 3     |
| `test_convergence/test_convergence_iterations.rs`  | 3     |
| `test_damping.rs`                                  | 2     |
| `test_topology_cleanup.rs`                         | 3     |
| `test_optimize_indel.rs`                           | 3     |
| `test_dense_sparse_equivalence/..._bounds.rs`      | 4     |
| `test_dense_sparse_equivalence/..._convergence.rs` | 2     |
| `test_dense_sparse_equivalence/..._validity.rs`    | 2     |

Files NOT parameterized:

- `test_eval_zero_branch_mismatch.rs` (1 call) -- tests prelude guard before dispatch, method-independent
- `test_optimize_zero_sequence_length.rs` (1 call) -- tests error check, method-independent
- `test_dispatch_zero_boundary.rs` -- already covers all 6 methods in `test_dispatch_zero_boundary_k80_identical_sequences` and 4 methods that cannot evaluate exactly at zero in `test_dispatch_zero_boundary_topology_cleanup_collects_k80_internal_edges`. Other tests in the file target specific inner solvers (`newton_inner`, `newton_sqrt_inner`) or the `reconcile_zero_boundary` helper directly; those are intentionally method-specific.

### Golden master switch

Change `test_gm_optimize.rs` from `BranchOptMethod::Newton` to `BranchOptMethod::BrentSqrt`.

`brent-sqrt` matches v0's algorithm and parameterization. Newton-t finds different optima on indel-bearing edges (Hessian dominance), making it unsuitable as a golden master method.

The reference data may need recapture if branch lengths differ between Newton and Brent on this dataset. Run the test first; if it fails, capture new reference with `brent-sqrt` and update fixtures. Only `flu_h3n2_20` is currently enabled; other datasets are commented out.

### Smoke test entries

Add one entry per non-default method to `OPTIMIZE_BATCH_RUNS` in `dev/run-smoke-tests`, all on `flu/h3n2/20`. The default `brent-sqrt` is covered by existing entries.

## Verification

`./dev/docker/run ./dev/dev t` -- all tests pass.

Smoke tests: `./dev/docker/run ./dev/dev br && ./dev/docker/run ./dev/run-smoke-tests .build/docker/release/treetime`

## Dependencies

- Depends on: nothing
- Depended on by: nothing (terminal)

## Cross-references

- Current `test_optimize_method.rs`: `packages/treetime/src/commands/optimize/__tests__/test_optimize_method.rs`
- Current `test_gm_optimize.rs`: `packages/treetime/src/commands/optimize/__tests__/test_gm_optimize.rs`
- Smoke tests: `dev/run-smoke-tests`
- Test support module: `packages/treetime/src/commands/optimize/__tests__/test_convergence/test_convergence_support.rs`
