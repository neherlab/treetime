# Mutation rendering and format projection are partially inconsistent

Tree output uses a shared zero-based `MutationEvent` representation, canonical substitution rendering, checked UShER MAT coordinate conversion, and shared tree writers. Amino-acid node-data duplicates part of that rendering contract, and two format projections remain undecided.

## Current contract

- `MutationEvent` represents substitutions, aligned insertions, and aligned deletions in one domain type [`packages/treetime/src/seq/mutation.rs#L51`](../../packages/treetime/src/seq/mutation.rs#L51).
- `fn mutation_event_strings()` performs checked one-based rendering for every event variant [`packages/treetime/src/seq/mutation.rs#L57`](../../packages/treetime/src/seq/mutation.rs#L57).
- UShER MAT conversion uses checked addition, checked `i32` narrowing, reference bounds checks, and explicit rejection of unsupported tracks and indels [`packages/treetime/src/commands/shared/tree_output.rs#L977`](../../packages/treetime/src/commands/shared/tree_output.rs#L977).

## Problems

- Amino-acid node-data has a private second implementation of insertion/deletion string expansion using unchecked arithmetic instead of the canonical helper [`packages/treetime/src/commands/ancestral/aa_node_data.rs#L308`](../../packages/treetime/src/commands/ancestral/aa_node_data.rs#L308).
- Amino-acid indel representation and coordinate parity with augur remain undecided in [N-amino-acid-mutation-indel-representation-undecided.md](N-amino-acid-mutation-indel-representation-undecided.md).
- PhyloXML mutation projection remains blocked on its format contract in [N-io-phyloxml-mutation-property-contract-undecided.md](N-io-phyloxml-mutation-property-contract-undecided.md).

## Required invariant

Let $p_0$ be the internal zero-based position and $p_1$ the one-based position required by an external format:

$$p_1 = p_0 + 1$$

The inverse is defined only for $p_1 \ge 1$:

$$p_0 = p_1 - 1$$

Every narrowing conversion must prove the target range. Shared spelling uses the canonical mutation API; format-specific projection remains explicit where wire contracts differ.

## Validation

- Replace the amino-acid duplicate with the canonical checked renderer after the amino-acid output policy is approved.
- Retain round-trip and boundary rejection tests for zero, overflow, and unrepresentable mutations.
- Compare every approved format projection with an independent reference fixture.

## Related issues

- [N-amino-acid-mutation-indel-representation-undecided.md](N-amino-acid-mutation-indel-representation-undecided.md)
- [N-io-phyloxml-mutation-property-contract-undecided.md](N-io-phyloxml-mutation-property-contract-undecided.md)
