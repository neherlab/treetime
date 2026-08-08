# Frozen coalescent lineage counts are an unsmoothed step function

The coalescent prior is built from lineage counts frozen before the refinement loop; see
[timetree-frozen-lineage-counts-for-coalescent-prior.md](../decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md).
`compute_lineage_counts` returns a `PiecewiseConstantFn` that steps at every merger, and
[`CoalescentModel`](../../packages/treetime/src/coalescent/coalescent.rs) evaluates $k(t)$ from it
pointwise.

Freezing removed the across-round instability the prior used to cause. It does not smooth the
prior's log-density *within* a pass: a node whose posterior straddles a step in $k(t)$ sees a
discontinuous change in the merger rate as it moves across it.

## Evidence

Weak. On `data/ebola/20` with a fixed $T_c$, an earlier measurement left one node of 39
alternating between two positions with `max/rms = 6.24 = sqrt(39)` — exactly one node moving while
the rest of the tree was still. Step-function $k(t)$ is a plausible cause of that local
bistability, but it was not isolated, and the configuration is no longer reproducible now that
fixed-$T_c$ runs converge.

## Impact

None demonstrated. Filed so the hypothesis is not rediscovered.

## Options

- Smooth $k(t)$ before building the prior, e.g. by convolution or by fitting a monotone
  interpolant through the step midpoints. The statistic role must keep the exact step function,
  since the $T_c$ solve integrates against it analytically.
- Leave it, and revisit only if a bistable node reappears with a `max/rms` ratio near
  $\sqrt{n_{\text{nodes}}}$.
