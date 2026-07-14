# Reject zero-state discrete partitions

`uniform_profile(n_states)` at [packages/treetime/src/partition/marginal_discrete.rs#L247](../../packages/treetime/src/partition/marginal_discrete.rs#L247) accepts `n_states == 0` and creates an empty profile with shape `(1, 0)`. A zero-state alphabet is invalid, and the constructor does not reject it.

## Impact

Production mugration code validates `n_states < 2` at [packages/treetime/src/mugration/mugration.rs#L118-L122](../../packages/treetime/src/mugration/mugration.rs#L118-L122). The lower-level constructor still needs to enforce its own state-space invariant.

## Affected code

- Constructor: [packages/treetime/src/partition/marginal_discrete.rs#L247](../../packages/treetime/src/partition/marginal_discrete.rs#L247)
- Guard: [packages/treetime/src/mugration/mugration.rs#L118-L122](../../packages/treetime/src/mugration/mugration.rs#L118-L122)

## Fix

- Make the constructor return an actionable error when `n_states == 0`.
- Add a unit test that asserts constructor rejection.
- Do not test the empty profile with `iter().all(...)`; that assertion is vacuously true for a `(1, 0)` array.

## Related issues

- Source: [kb/issues/M-discrete-missing-zero-states-inf.md](../issues/M-discrete-missing-zero-states-inf.md)
