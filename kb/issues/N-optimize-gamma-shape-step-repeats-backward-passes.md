# The gamma shape step repeats the backward passes of the marginal update

With `--gamma-categories K`, each optimize iteration runs the shape step (`fn PartitionSet::estimate_gamma_shapes()`, [packages/treetime/src/partition/set.rs#L205](../../packages/treetime/src/partition/set.rs#L205)) and then the marginal update. Every evaluation of the shape step runs the K backward passes of the partition, and the first evaluation is at the current shape. The marginal update after the step runs the K backward passes again at the selected shape. When the shape does not move, as at the lower bound 0.02 on `data/zika/86` and `data/mpox/clade-ii`, those passes repeat passes that the shape step already ran on the same branch lengths and topology.

On `data/mpox/clade-ii/100` with K = 4, the sparse backward passes take 53% of the optimize run, and reusing the passes of the selected shape saves K passes per iteration, about 15% of the run. Reuse needs a cache of the category passes keyed by the shape, the branch lengths and the topology, so that no stale pass reaches the update.

## Validation

- Optimize with `--gamma-categories 4` gives the same outputs, and the number of backward passes per iteration drops by K when the shape stays at a bound

## Proposed implementation plan

The shape step still discards its passes: the objective closure of `fn PartitionSet::estimate_gamma_shapes()` calls `mixture_log_lh()` and keeps only the log-likelihood, and `fn compute_iteration()` [packages/treetime/src/optimize/run_loop.rs#L244](../../packages/treetime/src/optimize/run_loop.rs#L244) then calls `PartitionSet::marginal_update()` [packages/treetime/src/partition/set.rs#L133](../../packages/treetime/src/partition/set.rs#L133), which runs the K backward passes again.

1. Profile first (`just profile` on `data/mpox/clade-ii/100` with `--gamma-categories 4`) to confirm the 53% and 15% figures on the current code.
2. Let the shape objective return the backward-pass state of each category together with its log-likelihood, and keep the state of the selected shape (the evaluation at the returned shape) in the result of `optimize_gamma_shape()`.
3. Pass that state into the marginal update of the same iteration as a typed value that carries the shape, the branch-length map, and a topology identity, and let `marginal_update()` run the backward passes only when the key differs. A key mismatch must fall back to a full pass, never reuse the state.
4. Tests: (a) the marginal update with the reused state gives the same profiles and log-likelihood as a fresh update, compared with `assert_ulps_eq!` (oracle: the fresh update on the same inputs, an independent code path for the forward results only when the reuse is a separate path); (b) a changed branch length or shape invalidates the state (the update result equals the fresh one).
5. Validate with `just test-rs`, and compare `just r treetime optimize --gamma-categories 4` outputs on `data/zika/86` before and after (byte-equal outputs), plus a new profile that shows K fewer backward passes per iteration.
