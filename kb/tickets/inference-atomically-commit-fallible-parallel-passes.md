# Atomically commit fallible parallel inference passes

Make marginal, optimize, and timetree parallel passes preserve their complete input state on error.

## Required changes

- Return keyed immutable deltas from workers instead of mutating shared payloads.
- Collect all deltas and errors before publication.
- Commit in deterministic graph-key order only after every worker succeeds.
- Run required post-pass/root corrections as part of the same transaction.

## Validation

- Inject failures at the first, middle, and last logical key under several worker counts.
- Compare whole graph/partition state before and after each failure.
- Assert successful results are identical across worker counts and repeated runs.
- Full lint and test suite.

## Related issues

- Source: [kb/issues/M-inference-fallible-parallel-passes-partially-commit.md](../issues/M-inference-fallible-parallel-passes-partially-commit.md)
