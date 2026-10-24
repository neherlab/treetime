# The root search misses positive-rate roots inside an edge

## Summary

With a required positive rate, the root search finds the unconstrained optimum of the root position on each edge, and then rejects the edge when the clock rate at that optimum is negative. It does not search the part of the edge where the rate is positive. A tree can have positive-rate root positions while every edge optimum has a negative rate, and the search then stops with "Clock rate is negative for all root positions".

## Counterexample

Tree `(A:3,B:2,C:1)R;` with dates A 0, B 1, C 2.

- At R, the divergences are A 3, B 2, C 1, so the rate is $-1$
- Root at distance $x$ from R on the edge to A ($0 \le x \le 3$): the divergences are $3 - x$, $2 + x$, $1 + x$. The least-squares rate is $x - 1$ and the residual sum of squares is $\tfrac{2}{3}x^2$, so the optimum of this edge is $x = 0$ (R), and the rate is positive for $x > 1$. At $x = 2$, the divergences are 1, 4, 3 and the rate is $+1$
- Root on the edge to B: the rate is $-1$ for every $x$
- Root on the edge to C: the rate is $-1 - x$ for every $x$

Every edge optimum is R with rate $-1$, so the search rejects all of them. A run of the least-squares root search (`fn reroot_at_best_root()`) on this tree returns the error "Clock rate is negative for all root positions. The data may lack temporal signal. Please specify --clock-rate explicitly." with `force_positive_rate = true`, and keeps R with rate $-1$ with `force_positive_rate = false`.

## Locations

- v1: `fn search_edges()` [packages/treetime/src/reroot/edge_search.rs#L19](../../packages/treetime/src/reroot/edge_search.rs#L19) optimizes each edge with `fn optimize_exact()` and rejects a negative-rate optimum at [packages/treetime/src/reroot/edge_search.rs#L56-L58](../../packages/treetime/src/reroot/edge_search.rs#L56-L58)
- v0: `fn TreeRegression.find_best_root()` [packages/legacy/treetime/treetime/treeregression.py#L321-L352](../../packages/legacy/treetime/treetime/treeregression.py#L321-L352) keeps the optimum of a branch from `_optimal_root_along_branch()` only when its slope is non-negative ([#L347](../../packages/legacy/treetime/treetime/treeregression.py#L347)), so v1 has the same behavior as v0
- [kb/algo/clock.md](../algo/clock.md), "Best Root Search", documents the v1 rule: the rate is checked on the optimum, not imposed as a constraint

> [!IMPORTANT]
> **Decision required.** The root search either keeps parity with v0, or searches the positive-rate part of each edge.
>
> - **Keep v0 parity**: an edge counts only when its unconstrained optimum has a positive rate. Data such as the counterexample stop with an error and need `--clock-rate`. This is the current behavior
> - **Search the positive-rate part of each edge**: on each edge, minimize the objective over the positions with a positive rate (the rate is a ratio of functions of the position, so the positive-rate set is a union of intervals), and take the best constrained optimum over all edges. The counterexample then gets a root inside the edge to A. This is a deviation from v0 that needs approval. When the constrained optimum lies at the boundary where the rate is zero, the result has a rate near zero, which [kb/algo/clock.md](../algo/clock.md) names as the reason for the current rule
>
> Evidence: the counterexample above, checked by hand and by a run of the v1 root search.
