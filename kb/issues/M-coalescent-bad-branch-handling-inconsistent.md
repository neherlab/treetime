# Coalescent bad_branch handling is inconsistent across code paths

Bad branches are excluded from some coalescent computations but included in others. V0 excludes bad branches from lineage events. The inconsistency means the lineage count function and edge likelihood see bad-branch contributions that the backward pass explicitly filtered out.

## Handling by code path

| Code path                                 | bad_branch handling               | Location               |
| ----------------------------------------- | --------------------------------- | ---------------------- |
| Backward pass: child message combination  | Excluded (`continue`)             | `backward_pass.rs:64`  |
| Backward pass: outgoing message to parent | Excluded (guarded)                | `backward_pass.rs:114` |
| Forward pass: leaf early return           | No check (returns for all leaves) | `forward_pass.rs:73`   |
| Lineage event collection                  | **Included** (no filter)          | `events.rs:27-43`      |
| Edge data collection                      | Excluded (guarded, PR#851)        | `edge_data.rs:55-62`   |
| Total log-likelihood                      | Excluded (via edge data, PR#851)  | `total_lh.rs:24-35`    |

## V0 behavior

V0 `TreeTime.make_time_tree()` builds coalescent events from nodes that pass the `bad_branch` filter. Bad branches are excluded from the lineage count and their edges are excluded from the coalescent likelihood.

## Impact

When bad branches exist in the tree:

- The lineage count $k(t)$ is inflated by bad-branch sampling events, increasing the per-branch merger rate $\kappa(t)$ and the expected merger count $H(t)$
- Bad-branch edges contribute to the total log-likelihood even though their time messages were excluded from the backward pass
- The Tc optimizer sees an objective function that mixes filtered and unfiltered tree structure

## Partial fix

[PR#851](https://github.com/neherlab/treetime/pull/851) added a `bad_branch()` guard to `collect_coalescent_edges`, fixing edge data collection and total log-likelihood. Lineage event collection (`events.rs`) remains unguarded -- bad-branch nodes still inflate $k(t)$.

## Remaining decision

Exclude bad branches from lineage events (`collect_tree_events`), matching v0. The lineage-count invariant (deltas sum to zero) needs adjustment to account for skipped subtrees.

## Related issues

- [H-timetree-coalescent-events-incomplete-after-topology-change.md](H-timetree-coalescent-events-incomplete-after-topology-change.md): event state is incomplete after topology changes (separate but same code path)
