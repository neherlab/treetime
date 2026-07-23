# Timetree refinement iteration mixes independent state transitions

`fn run_refinement_iteration()` applies relaxed-clock updates, resolves topology, reconciles every partition, selects inference order from topology dirtiness, compares ancestral states, and re-estimates the clock model [`packages/treetime/src/timetree/refinement.rs#L27`](../../packages/treetime/src/timetree/refinement.rs#L27).

## Hidden workflow

The operation mutates both graph and clock model through a fixed scientifically meaningful sequence. It accepts the pipeline's internal collaborators directly and returns an anonymous `(usize, usize)` pair for sequence differences and resolved nodes. Callers must recover those meanings through local destructuring names.

Topology mutation changes which partition reconciliation and inference operations are valid. The current function encodes that transition as local branching rather than as a state or result that subsequent steps must handle exhaustively.

## Required boundary

- Give each scientifically distinct transition a named operation over a complete refinement context.
- Represent topology reconciliation requirements explicitly after a topology-changing step.
- Return a named iteration outcome containing sequence differences, resolved-node count, topology status, and clock diagnostics required by convergence.
- Preserve the exact operation order, fallback behavior, and numerical policy unless a separate parity decision approves a change.

## Validation

- Golden-master comparisons cover refinement with and without topology changes.
- Unit tests verify each transition's preconditions and complete output state.
- Injected errors prove the selected graph/partition/clock atomicity contract.
- Callers no longer interpret tuple positions or infer topology status from counts.

## Related issues

- [H-core-command-module-shared-ops-entanglement.md](H-core-command-module-shared-ops-entanglement.md)
- [M-timetree-convergence-metric-deficiencies.md](M-timetree-convergence-metric-deficiencies.md)
- [M-inference-fallible-parallel-passes-partially-commit.md](M-inference-fallible-parallel-passes-partially-commit.md)
