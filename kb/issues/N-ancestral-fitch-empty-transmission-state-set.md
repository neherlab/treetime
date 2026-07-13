# Fitch transmission filtering can produce an empty state set

## Problem

The Fitch backward recurrence excludes a child's substitution state when the position lies in that edge's `transmission` ranges. If every child is excluded at one candidate position, both the intersection and union of the remaining child profiles are empty. The empty set is stored as a variable Fitch state and a later forward pass calls `StateSet::get_one()`, which has no valid state to select.

No production code currently assigns `SparseEdgePartition::transmission`, so the invalid state is dormant. Defining the correct result requires specifying whether a transmitted position contributes no evidence, inherits an external state, or should be excluded from the node entirely.

## Evidence

- Filtering and empty intersection/union: [`packages/treetime/src/ancestral/fitch_sub.rs`](../../packages/treetime/src/ancestral/fitch_sub.rs)
- Forward state selection: [`packages/treetime/src/ancestral/fitch_sub.rs`](../../packages/treetime/src/ancestral/fitch_sub.rs) (`resolve_root_forward`, `resolve_nonroot_substitutions_forward`)
- Transmission storage with no current writers: [`packages/treetime/src/partition/sparse.rs`](../../packages/treetime/src/partition/sparse.rs) (`SparseEdgePartition::transmission`)

## Required decision

Specify transmission semantics for Fitch substitution states, then reject or represent the all-children-filtered case without constructing an empty `StateSet`.
