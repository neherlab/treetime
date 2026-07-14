# Fitch transmission filtering can produce an empty state set

The Fitch backward recurrence excludes a child's substitution state when the position lies in that edge's `transmission` ranges [packages/treetime/src/ancestral/fitch_sub.rs#L43-L51](../../packages/treetime/src/ancestral/fitch_sub.rs#L43-L51). If every child is excluded at one candidate position, both the intersection and union of the remaining profiles are empty, so the backward pass stores an empty variable state [packages/treetime/src/ancestral/fitch_sub.rs#L61-L72](../../packages/treetime/src/ancestral/fitch_sub.rs#L61-L72). A later forward pass calls `StateSet::get_one()` without a valid state to select [packages/treetime/src/ancestral/fitch_sub.rs#L110-L119](../../packages/treetime/src/ancestral/fitch_sub.rs#L110-L119).

No production code currently assigns `SparseEdgePartition::transmission`, so the invalid state is dormant [packages/treetime/src/partition/sparse.rs#L34-L47](../../packages/treetime/src/partition/sparse.rs#L34-L47). Activating the field without first defining its evidence semantics would turn the dormant state into a failure or arbitrary reconstruction.

## Potential solutions

- O1. Treat a transmitted position as absent evidence from that child and represent the all-children-filtered case explicitly as missing evidence.
- O2. Require a transmitted position to inherit a state supplied by another partition or boundary object; reject construction when that state is unavailable.
- O3. Exclude the position from the node's Fitch candidate set when every child is filtered. This is valid only if exclusion means the node has no state obligation at that position.

## Recommendation

Define the biological and partition-boundary meaning of `transmission` before selecting an option. Then make the all-children-filtered state representable or reject it with a typed error; never construct an empty `StateSet`. No implementation ticket is ready while those semantics remain undecided.

## Related issues

- [M-ancestral-fitch-polytomy-recurrence-not-minimum.md](M-ancestral-fitch-polytomy-recurrence-not-minimum.md)
- [N-ancestral-sparse-remove-insert-pattern.md](N-ancestral-sparse-remove-insert-pattern.md)
