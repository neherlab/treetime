# Leaf and tip vocabulary drifts across phylogenetic operations

The generic graph layer consistently names terminal nodes leaves, and CLI concepts name sampled taxa tips. Core phylogenetic operations mix both terms for the same concept: `find_tip_group_root()` reports `no dated leaves`, while `complete_alignment_for_leaves()` reports missing tips [`packages/treetime/src/clock/reroot.rs#L190`](../../packages/treetime/src/clock/reroot.rs#L190) [`packages/treetime/src/ancestral/attach.rs#L12`](../../packages/treetime/src/ancestral/attach.rs#L12).

## Why the distinction matters

A graph leaf is defined by topology: it has no children. A phylogenetic tip is a sampled taxon and can carry names, dates, sequences, and metadata. These sets usually coincide, but algorithms and diagnostics should not rely on that equivalence implicitly.

Mixed names make it unclear whether an operation requires terminal topology or sampled biological data. That ambiguity is especially relevant to pruning, incomplete metadata, and attachment workflows.

## Required vocabulary

- Use `leaf` for graph-topology predicates and traversal results.
- Use `tip` for sampled taxa and user-facing phylogenetic concepts.
- At boundary operations, name both when the contract requires a leaf to have tip data.

## Validation

- Classify identifiers and diagnostics by graph or phylogenetic ownership before renaming.
- Preserve CLI flag names and externally documented terms unless the public contract is intentionally changed.
- Add explicit tests for terminal nodes missing sample metadata where an operation depends on the distinction.
