# Coalescent events are incomplete after topology changes

A topology-changing refinement pass clears internal time distributions. Coalescent event collection then silently filters nodes without `likely_time`, so the next inference step can omit internal events and continue with an incomplete prior.

`fn collect_tree_events()` adds an event only when both `time_distribution` and `likely_time()` are present [packages/treetime/src/coalescent/events.rs#L14-L43](../../packages/treetime/src/coalescent/events.rs#L14-L43); it checks only that the final event vector is nonempty, not that every required node contributed.

## Potential solutions

- O1. Rebuild complete time state immediately after topology mutation and make event collection reject missing state.
- O2. Derive event times directly during collection from distributions. This couples event construction to inference and still requires an explicit missing-state policy.

## Recommendation

Define complete node-time state as a precondition for event collection. After topology mutation, rebuild required time distributions and likely times before computing events. A missing required time returns a contextual error instead of removing the node from the event set.

## Validation

- Topology-changing and topology-preserving refinement cases.
- Assert the event multiset contains every required leaf, internal node, and root contribution.
- Inject one missing internal time and require an error naming the node.
- Compare the full coalescent likelihood against an independent direct calculation.
