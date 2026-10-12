# Edges in series are optimized as one sum

For a reversible model the likelihood depends only on the sum of edges joined by a node without data of its own and with exactly two incident edges: a bifurcating root or a single-child internal node (Felsenstein's pulley principle). v1 optimizes the sum of each such chain once per iteration and keeps the ratios of the chain's edge lengths. v0 does this for the bifurcating root only.

**Type**: Behavioral change (same estimand at the root; extended to single-child nodes).

**v0 location**: `optimize_tree_marginal()` at [`treeanc.py#L1317-L1337`](https://github.com/neherlab/treetime/blob/1a5c93d4f7929c3cc9302859015a32562dae02fd/packages/legacy/treetime/treetime/treeanc.py#L1317-L1337). Optimizes the root arc from the profiles of both subtrees, damps the total, distributes it by the previous ratio. The block runs once per root child, so the root arc is optimized and damped twice per iteration.

**v1 location**: `run_optimize_mixed_inner()` and `struct SeriesChain` in [`packages/treetime/src/optimize/dispatch.rs`](../../packages/treetime/src/optimize/dispatch.rs).

## Algorithm

1. Chains are the connected groups of edges joined through non-leaf nodes with exactly two incident edges. Leaves are the only nodes with data. A single-child root has one edge, whose length the likelihood does not depend on; it is not part of a chain.
2. Before the per-edge pass, each chain records its length ratios (equal split if the sum is zero) and picks its longest edge.
3. The longest edge is optimized with the other chain edges held at their current lengths, using the indel count of the whole chain. Since the likelihood depends only on the sum, the result plus the held lengths is the optimal sum for the current messages. The held edges are skipped in the per-edge pass.
4. The sum is distributed by the recorded ratios. Per-edge damping preserves the ratios.
5. If the optimized edge's optimum is zero, the sum is only bounded by the held lengths; the chain is set to zero, and the next iteration optimizes the sum exactly from zero.

The indel term is approximate: with a fixed ratio it depends on the sum, but evaluating the combined count on the optimized edge alone does not reproduce it exactly.

## Rationale

Optimizing the edges of a chain independently against the same messages gives $s_{n+1} = 2s^* - s_n$ for the sum: a neutral 2-cycle around the optimum $s^*$. With damping weight $d$ the error factor is $2d - 1$, so the loop converged only through heavy damping. With the joint sum, $s^*$ is reached in one undamped iteration for fixed messages.

## Evidence

Simulated data from `validation/` (fixed topology, perturbed start), damping 0, root sum relative to its converged value over iterations 1-4: independent updates 2-cycle (const rep00 L1000 mu0.01: 0.89, 1.12, 0.89, 1.11; exp rep07 L1000 mu0.01: 1.58, 0.65, 1.36, 0.65) and stop by the `Worsened` rule up to 5 log-lh below the best; joint sums are within 1.3% after one iteration and exact after one to three on all six sets tested. `test_root_preservation.rs` checks ratio preservation, the zero-start split, and convergence of the sum for the root and for a chain through a single-child node.
