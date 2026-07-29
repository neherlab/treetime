# Timetree refinement errors can leave partial state

`Refinement::run()` mutates the graph, partitions, and clock model through a scientifically ordered sequence. Missing or non-finite times encountered during polytomy scoring are rejected before topology mutation, but failures after the first successful mutation do not roll back earlier transitions.

## Impact

The command pipeline propagates the error and emits no result. An in-memory caller that catches the error can still observe a mixture of old and new refinement state and must discard the graph, partitions, and clock model rather than retrying with them.

This includes failures during node naming, partition reconciliation, marginal reconstruction, timetree inference, and clock-model estimation. Fallible parallel inference passes can publish a nondeterministic subset of worker results; that broader defect is tracked in [M-inference-fallible-parallel-passes-partially-commit.md](M-inference-fallible-parallel-passes-partially-commit.md).

## Required contract

Refinement needs an explicit publication boundary:

- Compute graph, partition, and clock-model deltas without publishing them, then commit all deltas after every transition succeeds.
- Alternatively, build complete transaction-local state and replace the caller's state only after success.

The first option preserves direct ownership without duplicating the complete graph and partition state. Tests must inject failures after topology mutation, partition reconciliation, inference, and clock estimation, then verify that every caller-visible input remains unchanged.
