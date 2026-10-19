# Prune chooses every internal edge to collapse before it collapses any

`fn prune_internal_nodes()` [packages/treetime/src/prune/prune.rs#L26-L66](../../packages/treetime/src/prune/prune.rs#L26-L66) first collects every internal edge that `--prune-short`, `--prune-empty` or a named node selects, from the branch lengths and mutation counts of the input, and then collapses the collected edges one by one. `fn collapse_edge()` [packages/treetime/src/optimize/topology/collapse.rs#L7-L42](../../packages/treetime/src/optimize/topology/collapse.rs#L7-L42) adds the length of the collapsed edge to each child edge and chains its mutations onto them. A child edge that was selected as short or empty is therefore collapsed even when the merge makes it longer than the threshold or gives it mutations.

v0 decides per node from the node's own branch length and does not add the removed length to the children [packages/legacy/treetime/treetime/treeanc.py#L1475-L1495](../../packages/legacy/treetime/treetime/treeanc.py#L1475-L1495).

> [!IMPORTANT]
> **Investigation required.** No wrong output has been shown. Find a tree where a selected child edge becomes longer than `--prune-short` or gains mutations after its parent edge collapses, and compare the result with v0 and with a version that checks each edge again after the collapses above it. Decide from that comparison whether the selection must use the merged values.
