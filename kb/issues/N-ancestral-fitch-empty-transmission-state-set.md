# Fitch transmission filtering can produce an empty state set

> [!WARNING]
> **Needs review.** The backward pass no longer stores an empty variable state. When every child is filtered, `fn resolve_variable_positions_backward()` now skips the position (`if child_profiles.is_empty() { continue; }`, [packages/treetime/src/partition/fitch/sub.rs#L52-L54](../../packages/treetime/src/partition/fitch/sub.rs#L52-L54), commit `219a014a`). The discovery pass ignores `transmission` and can already have written `VARIABLE_CHAR` at that position [packages/treetime/src/partition/fitch/sub.rs#L99-L107](../../packages/treetime/src/partition/fitch/sub.rs#L99-L107), so the node can keep the sentinel in its sequence with no `variable` entry. The forward passes iterate only `variable` [packages/treetime/src/partition/fitch/sub.rs#L118-L129](../../packages/treetime/src/partition/fitch/sub.rs#L118-L129). The effect of that leftover sentinel is not verified.

The Fitch backward recurrence excludes a child's substitution state when the position lies in that edge's `transmission` ranges [packages/treetime/src/partition/fitch/sub.rs#L36-L40](../../packages/treetime/src/partition/fitch/sub.rs#L36-L40). If every child is excluded at one candidate position, no child state set remains to intersect or unite [packages/treetime/src/partition/fitch/sub.rs#L52-L74](../../packages/treetime/src/partition/fitch/sub.rs#L52-L74), and the position has no valid state for the forward pass to select [packages/treetime/src/partition/fitch/sub.rs#L118-L184](../../packages/treetime/src/partition/fitch/sub.rs#L118-L184).

No production code currently assigns `SparseEdgeObs::transmission`, so the invalid state is dormant [packages/treetime/src/partition/storage/sparse.rs#L99](../../packages/treetime/src/partition/storage/sparse.rs#L99). Activating the field without first defining its evidence semantics would turn the dormant state into a failure or arbitrary reconstruction.

## Potential solutions

- O1. Treat a transmitted position as absent evidence from that child and represent the all-children-filtered case explicitly as missing evidence.
- O2. Require a transmitted position to inherit a state supplied by another partition or boundary object; reject construction when that state is unavailable.
- O3. Exclude the position from the node's Fitch candidate set when every child is filtered. This is valid only if exclusion means the node has no state obligation at that position.

## Recommendation

Define the biological and partition-boundary meaning of `transmission` before selecting an option. Then make the all-children-filtered state representable or reject it with a typed error; never construct an empty `StateSet` or leave a sentinel without a state set. Implementation waits until those semantics are decided.

## Related

- [kb/decisions/ancestral-fitch-plurality-on-multifurcations.md](../decisions/ancestral-fitch-plurality-on-multifurcations.md)
