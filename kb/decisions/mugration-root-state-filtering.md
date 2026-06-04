# Root state uniform-threshold filtering in mugration GTR inference

By default v1 matches v0: the root MAP state is always folded into the equilibrium-frequency prior, regardless of posterior confidence. v0 converts the root MAP state to a one-hot count via `self.tree.root.cseq` at [packages/legacy/treetime/treetime/treeanc.py#L1608-L1613](../../packages/legacy/treetime/treetime/treeanc.py#L1608-L1613).

v1 offers an opt-in `--filter-uninformative-root` flag (mugration only). When set, it skips the root state contribution to equilibrium frequency estimation for root positions whose posterior profile is near-uniform (max probability at or below `1/n_states + 1e-10`), at [packages/treetime/src/partition/marginal_core.rs#L249](../../packages/treetime/src/partition/marginal_core.rs#L249). The policy is carried on `MarginalData.filter_uninformative_root`, set per partition: `false` for mugration's `PartitionMarginalDiscrete` ([packages/treetime/src/partition/marginal_discrete.rs#L30](../../packages/treetime/src/partition/marginal_discrete.rs#L30)) unless the flag is given.

The nucleotide dense ancestral path keeps filtering enabled (`filter_uninformative_root: true` at [packages/treetime/src/partition/marginal_dense.rs#L52](../../packages/treetime/src/partition/marginal_dense.rs#L52)); that behavior is unchanged by this flag and is documented separately in [gtr-uninformative-root-state-filtering.md](gtr-uninformative-root-state-filtering.md).

## Impact

When the root is uninformative (posterior near-uniform across states) and the flag is set, v1 contributes a zero vector instead of v0's full count for whichever state is first in the argmax (state-order dependent), leaving pi estimation to transition counts and pseudo-counts alone. The effect is strongest for datasets with many states and weak root signal. With the flag off (the default), v1 reproduces v0's root-state contribution.

## Rationale

Filtering eliminates a state-order-dependent bias in the `root_state` vector that feeds `infer_gtr_impl()` pi estimation when the root carries no phylogenetic signal. It is opt-in because exact v0 parity is the default; the divergence occurs only when the flag is set.
