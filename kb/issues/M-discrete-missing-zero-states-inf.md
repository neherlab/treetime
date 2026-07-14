# Discrete partition constructor accepts zero states

`uniform_profile(n_states)` at [packages/treetime/src/partition/marginal_discrete.rs#L247](../../packages/treetime/src/partition/marginal_discrete.rs#L247) accepts `n_states == 0`. The resulting profile has shape `(1, 0)` and contains no elements, so any assertion that all elements are infinite would pass vacuously. The defect is acceptance of a zero-state model, not the contents of the empty array.

## Impact

Production code validates `n_states < 2` at [packages/treetime/src/mugration/mugration.rs#L118-L122](../../packages/treetime/src/mugration/mugration.rs#L118-L122) and returns an error, so this path is not reachable through mugration. The data-model constructor nevertheless admits an invalid state space and creates a partition on which later probability operations have no meaningful domain.

## Affected code

- Constructor: [packages/treetime/src/partition/marginal_discrete.rs#L247](../../packages/treetime/src/partition/marginal_discrete.rs#L247)
- Guard: [packages/treetime/src/mugration/mugration.rs#L118-L122](../../packages/treetime/src/mugration/mugration.rs#L118-L122)

## Fix

Reject `n_states == 0` at the constructor boundary. Tests must assert constructor rejection and must not characterize the empty array with an `all(...)` assertion.

## Related tickets

- [kb/tickets/discrete-guard-zero-states-in-missing-constructor.md](../tickets/discrete-guard-zero-states-in-missing-constructor.md)
