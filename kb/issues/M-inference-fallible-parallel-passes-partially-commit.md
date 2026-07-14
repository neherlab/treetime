# Fallible parallel inference passes partially commit state

Marginal reconstruction, branch optimization, and timetree branch-distribution passes mutate graph or partition state inside fallible parallel loops. If one worker fails, completed siblings remain committed; which siblings ran depends on scheduling.

Examples include fallible Rayon callbacks in the backward and forward marginal frontiers [packages/treetime/src/partition/marginal_core.rs#L99-L126](../../packages/treetime/src/partition/marginal_core.rs#L99-L126) [packages/treetime/src/partition/marginal_core.rs#L210-L240](../../packages/treetime/src/partition/marginal_core.rs#L210-L240).

## Impact

The command returns an error but leaves a nondeterministic hybrid of old and new scientific state. Retrying or inspecting the graph after failure is unsafe, and cleanup such as root correction may never run.

## Potential solutions

- O1. Compute immutable keyed deltas in parallel and commit after global success.
- O2. Mutate transaction-local graph copies and replace the original after success. This preserves atomicity but duplicates complete state.

## Recommendation

Workers compute immutable keyed deltas. Collect every result, return immediately on any error without mutation, then commit successful deltas in deterministic key order. This separates parallel computation from state publication.

## Related issues

- [N-ancestral-parallel-sparse-leaf-error-atomicity-unverified.md](N-ancestral-parallel-sparse-leaf-error-atomicity-unverified.md)
