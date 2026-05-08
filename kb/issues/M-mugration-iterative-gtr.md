# Mugration golden master parity with v0

v1 implements iterative GTR inference for mugration, matching v0's `reconstruct_discrete_traits()` ([packages/legacy/treetime/treetime/wrappers.py#L653-L811](../../packages/legacy/treetime/treetime/wrappers.py#L653-L811)) algorithm structure. Golden master trait assignments match v0 for 2/7 datasets (zika, zika with weights). The remaining 5 datasets diverge at ambiguous internal nodes due to two intentional v1 improvements over v0.

## Intentional v1 improvements (D1, D2)

1. Pseudo-count smoothing on initial pi ([packages/treetime/src/commands/mugration/run.rs#L221](../../packages/treetime/src/commands/mugration/run.rs#L221)): v1 applies `apply_pseudo_counts(pi, pc)` before the initial GTR model, giving a smoother prior for the first reconstruction. v0 passes raw weights directly to `GTR.custom()` and reserves `pc` for `infer_gtr()` regularization only.

2. Root state uniform-threshold filtering ([packages/treetime/src/commands/mugration/gtr_refinement.rs#L144-L150](../../packages/treetime/src/commands/mugration/gtr_refinement.rs#L144-L150)): v1 skips the root state contribution to pi estimation when the root posterior is near-uniform (max probability at or below `1/n_states + 1e-10`). v0 always converts the root MAP state to a one-hot count, which injects state-order bias when the root is uninformative. Matches the dense GTR inference path at [packages/treetime/src/gtr/infer_gtr/dense.rs#L158-L166](../../packages/treetime/src/gtr/infer_gtr/dense.rs#L158-L166).

Both produce scientifically defensible results but shift posterior probabilities at ambiguous internal nodes.

## Affected golden master tests

- `test_gm_mugration_outputs`: 2/7 passing (zika, zika_weights), 5/7 ignored (D1/D2 divergence)
- `test_gm_mugration_confidence_outputs`: 7/7 ignored (confidence profiles reflect the same differences)

## Related

- [Full forward-backward reconstruction proposal](../proposals/mugration-full-reconstruction-per-iteration.md)
