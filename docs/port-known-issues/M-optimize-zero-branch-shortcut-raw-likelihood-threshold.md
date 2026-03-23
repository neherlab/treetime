# Zero-branch shortcut underflows and depends on a fixed likelihood threshold

The unified optimizer decides whether zero branch length is optimal by multiplying raw zero-length likelihood factors and then checking whether the product is greater than `0.01`. This is numerically unstable and scientifically scale-dependent. The same likelihood shape can pass or fail the shortcut solely because the alignment has more sites or more ambiguous positions.

## Problem

The unified zero-branch shortcut is implemented in [`packages/treetime/src/commands/optimize/optimize_unified.rs#L92-L122`](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L92-L122) and [`packages/treetime/src/commands/optimize/optimize_unified.rs#L153-L173`](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L153-L173):

```rust
pub fn zero_branch_length_likelihood(&self) -> f64 {
  match self {
    OptimizationContribution::Dense(contribution) => contribution.coefficients.sum_axis(Axis(1)).product(),
    OptimizationContribution::Sparse(contribution) => contribution
      .site_contributions
      .iter()
      .map(|coeff| coeff.coefficients.sum().powf(coeff.multiplicity))
      .product(),
  }
}

pub fn is_zero_branch_optimal(contributions: &[OptimizationContribution]) -> bool {
  let zero_branch_length_lh: f64 = contributions
    .iter()
    .map(|contrib| contrib.zero_branch_length_likelihood())
    .product();

  if zero_branch_length_lh > 0.01 {
    let zero_branch_length_derivative: f64 = contributions
      .iter()
      .map(|contrib| contrib.zero_branch_length_derivative())
      .sum();

    if zero_branch_length_derivative < 0.0 {
      return true;
    }
  }

  false
}
```

Scientifically, the local criterion for whether `t = 0` is an optimum is the behavior of the log-likelihood near `t = 0`. For independent sites, log-likelihood is additive. Raw likelihood is multiplicative. Multiplying many small site terms is therefore the least stable representation of this decision problem.

Two distinct problems follow from the current implementation:

- numerical underflow: multiplying many small site likelihoods collapses to `0.0` in `f64`
- non-invariant decision rule: the cutoff `0.01` depends on sequence length and normalization scale, not on whether zero branch length is locally optimal

## Scientific context

For one branch, the optimizer works with per-site or per-pattern terms of the form:

`L_i(t) = sum_c k_ic exp(lambda_c t)`

The total log-likelihood is:

`log L(t) = sum_i log L_i(t)`

At `t = 0`, the meaningful local question is whether the derivative of `log L(t)` at zero is negative, zero, or positive. A fixed cutoff on the product `prod_i L_i(0)` is not part of that criterion. It only reflects how many terms were multiplied and how each term happened to be normalized.

## Examples

### Example 1. Underflow from many ambiguous sites

Suppose 1000 nucleotide sites are fully ambiguous at zero branch length, so each site contributes `L_i(0) = 0.25`.

The raw product is:

`0.25^1000 ≈ 9.3e-603`

That value underflows to `0.0` in `f64`. The current shortcut therefore skips the derivative check entirely, even though the only thing that changed was the number of multiplied sites.

### Example 2. Same local shape, different number of sites

Suppose one site has `L_1(0) = 0.5` and a negative zero-length derivative.

- with 1 site, the product is `0.5`, so the shortcut checks the derivative
- with 20 identical independent sites, the product is `0.5^20 ≈ 9.5e-7`, so the shortcut skips the derivative

The local shape of the log-likelihood did not change. Only the arbitrary raw scale changed. The decision rule should not flip because the same evidence was repeated more times.

## Blast radius

### Affected production code

- `zero_branch_length_likelihood()` in [`packages/treetime/src/commands/optimize/optimize_unified.rs#L92-L101`](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L92-L101)
- `zero_branch_length_derivative()` in [`packages/treetime/src/commands/optimize/optimize_unified.rs#L103-L122`](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L103-L122)
- `is_zero_branch_optimal()` in [`packages/treetime/src/commands/optimize/optimize_unified.rs#L153-L173`](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L153-L173)
- `run_optimize_mixed()` in [`packages/treetime/src/commands/optimize/optimize_unified.rs#L180-L250`](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L180-L250)
- `run_optimize()` in [`packages/treetime/src/commands/optimize/run.rs#L34-L150`](../../packages/treetime/src/commands/optimize/run.rs#L34-L150)
- the `treetime optimize` CLI entry point in [`packages/treetime-cli/src/bin/treetime.rs#L21-L70`](../../packages/treetime-cli/src/bin/treetime.rs#L21-L70)

### User-visible effect

The `optimize` command can fail to collapse a branch to zero even when the derivative criterion favors zero, or it can make the zero-branch decision depend on alignment size rather than on the branch likelihood shape. This affects convergence behavior and final branch lengths for edges near zero.

No current production timetree path calls `run_optimize_mixed()`. The blast radius is therefore narrower than the issue name suggests. The main production impact today is the standalone `optimize` command.

The thresholded raw-product shortcut in [`packages/treetime/src/commands/optimize/optimize_unified.rs#L153-L173`](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L153-L173) is still present in the current `dev` branch, so this issue remains current there.

### Affected and related tests

The current unit tests encode the heuristic instead of challenging it:

- [`packages/treetime/src/commands/optimize/__tests__/test_is_zero_branch_optimal.rs#L85-L95`](../../packages/treetime/src/commands/optimize/__tests__/test_is_zero_branch_optimal.rs#L85-L95) expects low raw zero-length likelihood to skip the derivative check
- [`packages/treetime/src/commands/optimize/__tests__/test_is_zero_branch_optimal.rs#L197-L218`](../../packages/treetime/src/commands/optimize/__tests__/test_is_zero_branch_optimal.rs#L197-L218) tests the `0.01` threshold boundary directly

Convergence tests such as [`packages/treetime/src/commands/optimize/__tests__/test_convergence/test_convergence_edge_cases.rs#L11-L40`](../../packages/treetime/src/commands/optimize/__tests__/test_convergence/test_convergence_edge_cases.rs#L11-L40) exercise zero-length inputs, but they do not test scale invariance or underflow resistance of the shortcut itself.

## Related code and exhaustive search

### Exact production occurrence

The production zero-branch shortcut used by the CLI exists in one path:

- [`packages/treetime/src/commands/optimize/optimize_unified.rs#L153-L173`](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L153-L173)

### Exact or closely similar sibling occurrences

Exhaustive search found two standalone helper implementations with the same raw-product pattern:

- standalone dense optimizer in [`packages/treetime/src/commands/optimize/optimize_dense.rs#L104-L123`](../../packages/treetime/src/commands/optimize/optimize_dense.rs#L104-L123)
- standalone sparse optimizer in [`packages/treetime/src/commands/optimize/optimize_sparse.rs#L160-L175`](../../packages/treetime/src/commands/optimize/optimize_sparse.rs#L160-L175)

The standalone sparse path is already tracked separately in [Standalone sparse optimizer skips derivative check for zero branches](N-optimize-sparse-zero-branch-no-derivative.md). The dense helper uses the same raw likelihood product and the same `0.01` threshold as the unified path, so it shares the underflow and scale-dependence defect even though it does check the derivative sign.

### Related documentation and tests

- [`docs/port-intentional-changes/optimize-newton-raphson-per-edge.md#L33`](../port-intentional-changes/optimize-newton-raphson-per-edge.md#L33) documents the current short-circuit as a property of the optimizer. The current behavior is real, but the threshold-based decision itself is not scientifically justified.
- [`packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_bounds.rs#L50-L57`](../../packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_bounds.rs#L50-L57) already acknowledges that different zero-branch thresholds perturb dense-sparse equivalence tests.
- [`packages/treetime/src/commands/optimize/__tests__/test_is_zero_branch_optimal.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_is_zero_branch_optimal.rs) still encodes the threshold rule in the current `dev` branch. The tests do not yet assert scale invariance or underflow resistance.
- [Timetree crashes on zero-length branches with grid spacing error](H-timetree-crash-grid-zero-branch.md) is a separate issue in the timetree distribution code. It is related by symptom and domain, not by direct call path.
- [Standalone helper optimizers duplicate production logic and drift from unified behavior](N-optimize-standalone-helper-optimizers-duplicate-production-logic.md) tracks the duplication that let the same zero-branch design defect persist in multiple optimizers.

## Proposed solutions

### S1. Use derivative sign with per-site safety checks and log-space aggregation

For each site or pattern at `t = 0`:

- compute `site_lh = L_i(0)`
- if any `site_lh` is non-finite or smaller than a small epsilon, do not use the shortcut
- otherwise, accumulate `log_zero_lh += ln(site_lh)` for diagnostics only
- decide zero optimality from `zero_branch_length_derivative() < 0.0`

Why this works:

- it removes the unstable raw product
- it keeps the decision criterion tied to the local slope of the log-likelihood
- it makes the shortcut scale-invariant with respect to the number of sites

Pros:

- numerically stable
- easy to reason about
- preserves the current derivative machinery
- cleanly falls back to full optimization when `t = 0` is too degenerate to evaluate safely

Cons:

- requires choosing an epsilon for “too close to zero”
- the shortcut becomes conservative for highly degenerate cases because it may refuse to decide and fall back to full optimization

### S2. Rescale each coefficient row before evaluating the zero-length shortcut

Before computing `L_i(0)`, divide each coefficient row by a positive scale such as its row maximum. Positive row scaling changes the raw likelihood scale but does not change:

- the sign of the zero-length derivative ratio
- the position of the branch-length optimum

Then remove the `0.01` cutoff and use derivative sign after confirming the rescaled site terms are finite.

Pros:

- suppresses most underflow at the source
- keeps the per-row arithmetic in ordinary floating-point space

Cons:

- more invasive than S1
- requires careful documentation because the rescaled raw likelihood no longer has the old interpretation

### S3. Remove the zero-branch shortcut entirely

Delete the pre-check and always rely on Newton or grid search.

Pros:

- simplest semantics
- no special-case threshold logic

Cons:

- loses a useful fast path when the branch is obviously zero
- may increase runtime for large trees with many zero or near-zero branches

### Recommended direction

S1 is the best fix. It removes the broken threshold rule while preserving the useful optimization shortcut and the existing derivative calculations.

## Implementation notes

- The main evaluation path already uses log likelihood accumulation in [`packages/treetime/src/commands/optimize/optimize_dense_eval.rs#L31-L42`](../../packages/treetime/src/commands/optimize/optimize_dense_eval.rs#L31-L42) and [`packages/treetime/src/commands/optimize/optimize_sparse_eval.rs#L27-L39`](../../packages/treetime/src/commands/optimize/optimize_sparse_eval.rs#L27-L39). This issue is specifically about the zero-branch pre-check, not the whole optimizer.
- A pure derivative-only shortcut without safety checks is risky because the current derivative formulas divide by zero-length site sums in [`packages/treetime/src/commands/optimize/optimize_unified.rs#L107-L122`](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L107-L122). Any fix that removes the raw-product threshold should replace it with explicit per-site validity checks or a deterministic fallback when a denominator is zero or non-finite.
- If the standalone dense and sparse helper optimizers remain, they should share the corrected zero-branch helper instead of carrying their own inline copies.

## Proposed tests and verification

### T1. Scale-invariance test

Create two synthetic contributions with the same per-site zero-length derivative behavior:

- one with a single site
- one with many duplicated copies of that site

Expected result:

- `is_zero_branch_optimal()` should return the same answer for both
- current code fails this when the raw product crosses the `0.01` threshold

### T2. Underflow-resistance test

Create many zero-length site terms equal to `0.25` or another small positive value that would underflow when multiplied directly.

Expected result:

- the new implementation must not underflow to `0.0`
- if the shortcut declines to decide because a site is too close to zero, it must do so explicitly and deterministically

### T3. Replace threshold-boundary tests

Rewrite [`packages/treetime/src/commands/optimize/__tests__/test_is_zero_branch_optimal.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_is_zero_branch_optimal.rs) so it validates:

- derivative sign handling
- invariance under repeated independent sites
- behavior when one site likelihood is near zero

The tests that assert the `0.01` boundary itself should be removed because they enforce the bug.

### T4. Align helper optimizers or mark them test-only

The standalone dense and sparse helpers should be verified against the same scale-invariant rule. If they are kept, they should either:

- adopt the same corrected zero-branch logic as the unified path
- or be isolated as test-only helpers so their divergent shortcuts do not confuse equivalence tests

On `dev`, this matters for the current helper-based tests too. If the helpers keep the old threshold logic while the unified path is fixed, dense-sparse equivalence tests will need to state explicitly which optimizer they are validating.
