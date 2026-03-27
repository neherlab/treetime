# Root state uniform-threshold filtering in mugration GTR inference

v1 skips the root state contribution to equilibrium frequency estimation when the root posterior profile is near-uniform at [packages/treetime/src/commands/mugration/gtr_refinement.rs#L138-L145](../../packages/treetime/src/commands/mugration/gtr_refinement.rs#L138-L145). A root profile with max probability at or below `1/n_states + 1e-10` carries no phylogenetic signal and is excluded.

v0 always converts the root MAP state to a one-hot count via `self.tree.root.cseq` at [packages/legacy/treetime/treetime/treeanc.py#L1608-L1613](../../packages/legacy/treetime/treetime/treeanc.py#L1608-L1613), regardless of posterior confidence.

## Impact

When the root is uninformative (posterior near-uniform across states), v0 injects a full count for whichever state happens to be first in the argmax (state-order dependent). v1 contributes a zero vector instead, leaving pi estimation to transition counts and pseudo-counts alone. The effect is strongest for datasets with many states and weak root signal.

## Rationale

Matches the dense GTR inference path's approach at [packages/treetime/src/gtr/infer_gtr/dense.rs#L158-L166](../../packages/treetime/src/gtr/infer_gtr/dense.rs#L158-L166), which already filters uninformative root profiles. Eliminates state-order bias in the `root_state` vector that feeds into `infer_gtr_impl()` pi estimation.
