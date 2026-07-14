# Add property tests for ClockSet algebraic identities and Fitch parsimony invariants

## No property tests for ClockSet algebraic identities

`packages/treetime/src/payload/clock_set.rs:53-172:`

`+`, `-`, `+=`, `-=`, `fn propagate_averages` lack property test coverage for algebraic identities (associativity, commutativity, identity element).

## No property tests for Fitch parsimony invariants

Score invariant under rerooting, state-set subset relation between parent and child Fitch sets.

## Related issues

- Source: [kb/issues/N-test-coverage-gaps.md](../issues/N-test-coverage-gaps.md) -- delete after full resolution
