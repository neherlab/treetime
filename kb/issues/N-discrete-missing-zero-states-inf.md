# Discrete partition constructor accepts zero states

`PartitionMarginalDiscrete::new()` [packages/treetime/src/partition/marginal/discrete/partition.rs#L27-L61](../../packages/treetime/src/partition/marginal/discrete/partition.rs#L27-L61) accepts a `DiscreteStates` set with no states. Its leaf profiles come from `fn missing_trait_profile()` and `fn one_hot_profile()` [packages/treetime/src/partition/marginal/discrete/input.rs#L12-L20](../../packages/treetime/src/partition/marginal/discrete/input.rs#L12-L20), which then have shape `(1, 0)` and contain no elements, so any assertion over all elements would pass vacuously. The defect is acceptance of a zero-state model, not the contents of the empty array.

## Impact

Production code validates `n_states < 2` at [packages/treetime/src/mugration/pipeline.rs#L72-L76](../../packages/treetime/src/mugration/pipeline.rs#L72-L76) and returns an error, so this path is not reachable through mugration. The data-model constructor nevertheless admits an invalid state space and creates a partition on which later probability operations have no meaningful domain.

## Affected code

- Constructor: [packages/treetime/src/partition/marginal/discrete/partition.rs#L27-L61](../../packages/treetime/src/partition/marginal/discrete/partition.rs#L27-L61)
- Profile helpers: [packages/treetime/src/partition/marginal/discrete/input.rs#L12-L20](../../packages/treetime/src/partition/marginal/discrete/input.rs#L12-L20)
- Guard: [packages/treetime/src/mugration/pipeline.rs#L72-L76](../../packages/treetime/src/mugration/pipeline.rs#L72-L76)

## Fix

Reject a zero-state set at the constructor boundary with an actionable error. A unit test must assert constructor rejection and must not characterize the empty array with an `all(...)` assertion, which is vacuously true for a `(1, 0)` array.
