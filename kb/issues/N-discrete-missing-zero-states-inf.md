# Discrete partition constructor accepts zero states

`PartitionMarginalDiscrete::new()` [packages/treetime/src/partition/marginal/discrete/partition.rs#L21-L52](../../packages/treetime/src/partition/marginal/discrete/partition.rs#L21-L52) accepts a `DiscreteStates` set with no states. Its leaf profiles come from `fn missing_trait_profile()` and `fn one_hot_profile()` [packages/treetime/src/partition/marginal/discrete/input.rs#L10-L18](../../packages/treetime/src/partition/marginal/discrete/input.rs#L10-L18), which then have shape `(1, 0)` and contain no elements, so any assertion over all elements would pass vacuously. The defect is acceptance of a zero-state model, not the contents of the empty array.

## Impact

The constructor is `pub(crate)`, and its only production caller is the mugration pipeline [packages/treetime/src/mugration/pipeline.rs#L106](../../packages/treetime/src/mugration/pipeline.rs#L106), which validates `n_states < 2` first at [packages/treetime/src/mugration/pipeline.rs#L74-L78](../../packages/treetime/src/mugration/pipeline.rs#L74-L78) and returns an error, so this path is not reachable through mugration. The data-model constructor nevertheless admits an invalid state space and creates a partition on which later probability operations have no meaningful domain.

## Affected code

- Constructor: [packages/treetime/src/partition/marginal/discrete/partition.rs#L21-L52](../../packages/treetime/src/partition/marginal/discrete/partition.rs#L21-L52)
- Profile helpers: [packages/treetime/src/partition/marginal/discrete/input.rs#L10-L18](../../packages/treetime/src/partition/marginal/discrete/input.rs#L10-L18)
- Guard: [packages/treetime/src/mugration/pipeline.rs#L74-L78](../../packages/treetime/src/mugration/pipeline.rs#L74-L78)

## Fix

Reject a zero-state set at the constructor boundary with an actionable error. A unit test must assert constructor rejection and must not characterize the empty array with an `all(...)` assertion, which is vacuously true for a `(1, 0)` array.
