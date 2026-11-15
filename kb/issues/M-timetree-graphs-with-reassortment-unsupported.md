# Timetree cannot run on graphs with reassortment

Time inference requires a tree: the forward and backward time passes call `require_tree()` ([packages/treetime/src/timetree/inference/backward_pass.rs#L36](../../packages/treetime/src/timetree/inference/backward_pass.rs#L36), [packages/treetime/src/timetree/inference/forward_pass.rs#L42](../../packages/treetime/src/timetree/inference/forward_pass.rs#L42)), as do clock regression and rerooting. v0 `treetime arg` dates reassortant segment trees by analysing one focal tree at a time with masked branches ([kb/features/arg.md](../features/arg.md)).

On a graph with reassortment a node has one time but may have several parents, so the time graph has undirected cycles and two message passes no longer give exact marginals. The options and the question of shared timetree parameters are in [kb/proposals/edge-scoped-inheritance.md](../proposals/edge-scoped-inheritance.md#time-trees).

> [!IMPORTANT]
> **Decision required.** Loopy message passing, conditioning on one parent per reassortant node in turn (v0), or exact treatment of the cycles; see the proposal for the trade-offs.

## Depends on

- [M-partition-passes-require-tree-topology.md](M-partition-passes-require-tree-topology.md): branch-length likelihoods per edge come from the scope-aware sequence passes
