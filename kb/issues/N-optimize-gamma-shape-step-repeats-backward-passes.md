# The gamma shape step repeats the backward passes of the marginal update

With `--gamma-categories K`, each optimize iteration runs the shape step (`fn PartitionSet::estimate_gamma_shapes()`, [packages/treetime/src/partition/set.rs#L205](../../packages/treetime/src/partition/set.rs#L205)) and then the marginal update. Every evaluation of the shape step runs the K backward passes of the partition, and the first evaluation is at the current shape. The marginal update after the step runs the K backward passes again at the selected shape. When the shape does not move, as at the lower bound 0.02 on `data/zika/86` and `data/mpox/clade-ii`, those passes repeat passes that the shape step already ran on the same branch lengths and topology.

On `data/mpox/clade-ii/100` with K = 4, the sparse backward passes take 53% of the optimize run, and reusing the passes of the selected shape saves K passes per iteration, about 15% of the run. Reuse needs a cache of the category passes keyed by the shape, the branch lengths and the topology, so that no stale pass reaches the update.

## Validation

- Optimize with `--gamma-categories 4` gives the same outputs, and the number of backward passes per iteration drops by K when the shape stays at a bound
