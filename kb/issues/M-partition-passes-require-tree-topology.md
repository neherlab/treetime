# Partition passes require tree topology

The Fitch and marginal passes (dense, discrete, sparse) accept only trees: they call `require_tree()` at their start and read only the first parent of a node. A graph that combines several segment trees, with nodes that inherit different segments from different parents, therefore cannot be reconstructed, although the graph and its schedule accept it ([kb/algo/graph.md](../algo/graph.md)).

## Locations

- Tree check: `fn require_tree()` [packages/treetime-graph/src/schedule.rs#L54](../../packages/treetime-graph/src/schedule.rs#L54), called by [packages/treetime/src/partition/fitch/passes.rs#L97](../../packages/treetime/src/partition/fitch/passes.rs#L97), [#L199](../../packages/treetime/src/partition/fitch/passes.rs#L199), [packages/treetime/src/partition/marginal/shared/pass.rs#L98](../../packages/treetime/src/partition/marginal/shared/pass.rs#L98), [#L269](../../packages/treetime/src/partition/marginal/shared/pass.rs#L269), and the sparse backward and forward passes ([packages/treetime/src/partition/marginal/sparse/backward.rs](../../packages/treetime/src/partition/marginal/sparse/backward.rs), [packages/treetime/src/partition/marginal/sparse/forward.rs](../../packages/treetime/src/partition/marginal/sparse/forward.rs))
- First parent only: `parents.first()` and `parents(key).first()` in the same files and in [packages/treetime/src/partition/marginal/sparse/sequences.rs#L23](../../packages/treetime/src/partition/marginal/sparse/sequences.rs#L23)
- Root of all sites at a node without parents: [packages/treetime/src/partition/marginal/shared/pass.rs#L193](../../packages/treetime/src/partition/marginal/shared/pass.rs#L193), [packages/treetime/src/partition/marginal/sparse/backward.rs#L90](../../packages/treetime/src/partition/marginal/sparse/backward.rs#L90)

## Fix direction

Edge scopes as in [kb/proposals/edge-scoped-inheritance.md](../proposals/edge-scoped-inheritance.md): each edge carries a set of sites per partition, every site has at most one parent at every node, and the passes combine and send messages only over the scopes, with root terms at the root sites of each node. Order: Fitch, dense marginal, discrete, sparse; whole-partition scopes first, ranges within a partition after.

## Validation

- Tree inputs (full scopes everywhere): all outputs byte-identical to the current passes
- Toy reassortment graph ([data/toy/reassortment/](../../data/toy/reassortment/)): per partition, states, mutations, and log-likelihood equal to a run on that segment's tree alone; the total log-likelihood equals the sum of the per-segment runs
- Analytic case: two leaves under Jukes-Cantor with one partition rooted at a different node than the other
