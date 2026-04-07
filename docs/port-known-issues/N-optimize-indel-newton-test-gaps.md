# Missing end-to-end test coverage for indel-aware Newton optimizer

Three test gaps in the indel-aware Newton optimization path prevent verification of scientific correctness. The gaps interact with a production convergence issue ([M-optimize-newton-indel-hessian-dominance](M-optimize-newton-indel-hessian-dominance.md)).

## Gap 1: No end-to-end optimality verification

Existing tests call `run_optimize_mixed()` on indel-bearing edges but assert only positivity and finiteness:

- `test_optimize_indel_run_optimize_nonzero_with_indels` at [packages/treetime/src/commands/optimize/**tests**/test_optimize_indel.rs#L200](../../packages/treetime/src/commands/optimize/__tests__/test_optimize_indel.rs#L200): asserts `bl > 0.0 && bl.is_finite()`
- `test_optimize_indel_zero_bl_pipeline_escapes_zero` at [packages/treetime/src/commands/optimize/**tests**/test_optimize_indel.rs#L232](../../packages/treetime/src/commands/optimize/__tests__/test_optimize_indel.rs#L232): same assertions

No test verifies that the reported branch length is at a local maximum of the combined (substitution + indel) log-likelihood. An attempted local-maximum test (evaluating combined log-likelihood at perturbed points around the optimum) correctly failed, exposing the premature convergence bug in [M-optimize-newton-indel-hessian-dominance](M-optimize-newton-indel-hessian-dominance.md).

This gap cannot be closed until the convergence bug is fixed. Once fixed, the test should:

1. Run `run_optimize_mixed()` on an edge with indels
2. Evaluate the combined log-likelihood at the reported optimum and at perturbed points ($t \pm \delta$)
3. Assert the optimum has the highest combined log-likelihood

The helper-level test `test_optimize_indel_newton_converges_to_poisson_mle` at [packages/treetime/src/commands/optimize/**tests**/test_optimize_indel.rs#L441](../../packages/treetime/src/commands/optimize/__tests__/test_optimize_indel.rs#L441) verifies only that the Poisson helper function has the expected analytical properties (derivative zero at MLE, Newton step converges). It does not exercise `run_optimize_mixed()`.

## Gap 2: Zero-branch mismatch precondition not pinned

The regression test `test_eval_zero_branch_mismatch_no_nan` at [packages/treetime/src/commands/optimize/**tests**/test_eval_zero_branch_mismatch.rs#L25](../../packages/treetime/src/commands/optimize/__tests__/test_eval_zero_branch_mismatch.rs#L25) exercises the `branch_length == 0.0` path through `run_optimize_mixed()` and asserts finite, non-negative output. It does not verify that any edge contribution has `all_sites_valid_at_zero() == false` before optimization.

Investigation showed that marginal reconstruction with JC69 at zero branch length spreads probability mass across states via the transition matrix, so even fully mismatched leaf sequences (A vs C vs G vs T) produce overlapping profiles after the marginal backward-forward passes. The predicate `!all_sites_valid_at_zero()` may not be reachable through the standard marginal reconstruction pipeline with JC69.

The production guard at [packages/treetime/src/commands/optimize/optimize_unified.rs#L394](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L394) protects against this condition:

```rust
if branch_length == 0.0 && !contributions.iter().all(|c| c.all_sites_valid_at_zero()) {
  branch_length = one_mutation;
}
```

If the condition is unreachable through marginal reconstruction, the guard is defense-in-depth for non-standard callers (direct API users constructing contributions without marginal reconstruction). To test it, either:

- Construct `OptimizationContribution` directly with hand-crafted coefficients that sum to zero, bypassing marginal reconstruction
- Determine whether the guard is unreachable in all production paths and document it as defense-in-depth with a rationale comment

## Gap 3: Finite-difference Hessian tolerance lacks error analysis

The proptest finite-difference second derivative check at [packages/treetime/src/commands/optimize/**tests**/test_coefficient_extraction_dense/test_coefficient_extraction_dense_prop_invariants.rs#L240-L257](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_prop_invariants.rs#L240-L257) uses `max_relative = 1e-3`. Binary search showed `1e-3` passes and `1e-4` fails across the proptest input space, but the tolerance is not justified by error analysis.

For the central difference approximation $f''(x) \approx (f(x+h) - 2f(x) + f(x-h)) / h^2$, the total error is $O(h^2) + O(\epsilon / h^2)$ where $\epsilon \approx 10^{-16}$ is machine epsilon. The optimal step size minimizing total error is $h_{\text{opt}} = \epsilon^{1/4} \approx 1.2 \times 10^{-4}$, giving theoretical best relative error $O(\epsilon^{1/2}) \approx 10^{-8}$.

The current test uses $h = t \times 10^{-3}$ with $t \in [0.01, 1.0]$, so $h \in [10^{-5}, 10^{-3}]$. At $h = 10^{-5}$, the round-off term $\epsilon / h^2 \approx 10^{-6}$ and the truncation term $h^2 \approx 10^{-10}$, so the theoretical relative error is $\sim 10^{-6}$. The `1e-3` tolerance is 1000x looser than this bound.

The gap between theoretical and observed error suggests either: (a) the analytical second derivative has a subtle numerical issue that amplifies error for specific coefficient/branch-length combinations, or (b) the step size is suboptimal for the curvature of the eigenvalue-space log-likelihood.

To close this gap:

1. Use the theoretically optimal step $h = (\epsilon \cdot t)^{1/4}$ instead of $h = t \times 10^{-3}$
2. Measure actual relative errors across the proptest input space
3. Set tolerance to 10x the measured worst case
4. If the measured error exceeds the theoretical bound, investigate whether the analytical second derivative has a numerical issue (e.g. catastrophic cancellation in the Hessian formula)

## Impact

Negligible individually. Gap 1 blocks verification of the optimizer's core scientific claim (branch lengths at local maxima of the combined likelihood). Gaps 2 and 3 are test quality issues.

## Related

- [M-optimize-newton-indel-hessian-dominance](M-optimize-newton-indel-hessian-dominance.md) -- Gap 1 is blocked by this production convergence bug
- [Dense-sparse duplication in OptimizationContribution enum dispatch](L-optimize-eval-dense-sparse-duplication.md) -- the enum dispatch duplication affects the same `all_sites_valid_at_zero()` method referenced in Gap 2
