# The short-branch collapse adds the collapsed length to the child edges

## Deviation

Before the timetree loop, v1 collapses every non-root internal edge that is shorter than a tenth of a mutation ($0.1/L$, $L$ the number of sites) and whose reconstructed states equal the states of its parent (`fn collapse_short_internal_edges()` in [packages/treetime/src/optimize/topology/collapse.rs](../../packages/treetime/src/optimize/topology/collapse.rs), called from [packages/treetime/src/timetree/pre_loop.rs](../../packages/treetime/src/timetree/pre_loop.rs)). The children of the removed node move to its parent, and `fn collapse_edge()` adds the length of the removed edge to the edge of each child. The `optimize` command collapses with the same function after its loop.

## v0 behavior

`TreeAnc.prune_short_branches()` moves the children to the parent and drops the length of the removed edge: each child keeps its own branch length ([packages/legacy/treetime/treetime/treeanc.py#L1475-L1495](../../packages/legacy/treetime/treetime/treeanc.py#L1475-L1495), the reassignment at [L1492](../../packages/legacy/treetime/treetime/treeanc.py#L1492)). `TreeTime._run()` calls it before the clock filter and the first reroot, and again after the second branch-length optimization in marginal mode ([packages/legacy/treetime/treetime/treetime.py#L235-L266](../../packages/legacy/treetime/treetime/treetime.py#L235-L266)).

## Rationale

The length of the removed edge is part of the distance from the root to every tip below it. v1 keeps every root-to-tip distance unchanged by the collapse. In v0 the root-to-tip distance of each tip below a collapsed edge becomes shorter by the length of that edge, and by the sum of the lengths when several collapsed edges lie on one path. The collapse removes an edge because its length is below the resolution of the alignment and its states do not change, which is a decision about the topology, not evidence that the distance is zero.

## Effect

- Each child edge of a collapsed node is longer than in v0 by less than $0.1/L$, the length of the removed edge
- Input mode: the input branch lengths are the data of the timetree loop (`InputBranchModel`), so the branch likelihoods of these child edges differ from v0
- Both modes: the clock filter, the reroot regressions, and the Gaussian clock estimate at the start of the loop read root-to-tip distances that are longer than in v0 by the sum of the collapsed lengths on each path
- Marginal mode: the branch likelihoods of the loop come from the ancestral reconstruction, not from these lengths, and the ML branch-length optimization after the first reroot optimizes every edge again
