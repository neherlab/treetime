# Test mugration transition counting against hand-computed values

Exercise `fn MarginalData::count_transitions()` [`packages/treetime/src/partition/marginal_core.rs#L483-L521`](../../packages/treetime/src/partition/marginal_core.rs#L483-L521) through the discrete partition on a small discrete-state tree, complementing the existing shared contract tests.

## Acceptance criteria

- A hand-computed two-state fixture verifies the full expected transition matrix, dwell-time vector, and root state.
- The hand oracle distinguishes near-uniform root filtering and verifies discrete delegation; existing tests remain the oracle for parent/child orientation, diagonal zeroing, branch-length weighting, and dense/sparse equality.
- The oracle does not call the production accumulation helpers under test.

## Related issues

- Source: [kb/issues/N-mugration-count-transitions-untested.md](../issues/N-mugration-count-transitions-untested.md)
