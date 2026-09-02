# Optimize loop resolves reversion-driven polytomies

The `optimize` loop adds a topology move with no v0 counterpart: when a polytomy forces a substitution onto an internal edge and then reverts it on a child, the move relocates the child so the reversion cancels. This removes homoplasy that an arbitrary binary resolution of a polytomy introduces, and that v0's loop never targets.

## What v0 does

v0's branch-length optimization loop (`TreeAnc.optimize_tree_marginal`) collapses short internal branches (`prune_short_branches`) but has no step that inspects or removes reversions. A reversion is only ever dropped as an incidental side effect when an edge is collapsed for an unrelated reason and its substitutions compose with the child's. v0 has no shared-mutation merge and no reversion-aware polytomy handling.

## What v1 does

Once per iteration, [`prune_and_merge_in_loop()`](../../packages/treetime/src/optimize/run_loop.rs#L319) runs a single per-polytomy routine, [`resolve_polytomies()`](../../packages/treetime/src/optimize/topology/resolve_polytomy.rs#L31), after the zero-optimal edge collapse. Sparse partitions only; dense partitions carry no per-edge mutation lists. The routine applies three steps at each polytomy and iterates to a fixpoint:

- Merge siblings that share substitutions under a helper node ([`merge_single_polytomy()`](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs#L64)). This canonicalizes several children reverting the same position into one reverting child.
- Hoist the reverting child with the largest reversion count ([`hoist_reverting_child()`](../../packages/treetime/src/optimize/topology/hoist_reversions.rs#L74)), ties broken by edge key.
- Retire helper nodes: collapse mutation-free edges to nodes the routine created in this pass.

The hoist inserts a new node $N$ between parent $u$ and node $v$, grouping $v$ with one reverting child $c$. The parent-edge substitutions $M_v$ partition by position against the child substitutions $M_c$ into three disjoint sets: $T$ (untouched by the child), $H$ (chained, $a \to b$ then $b \to d$), and $R$ (reverted, $a \to b$ then $b \to a$). The resulting edges carry:

| edge      | substitutions                                                                                          |
| --------- | ------------------------------------------------------------------------------------------------------ |
| $u \to N$ | $T$                                                                                                    |
| $N \to v$ | $H \cup R$                                                                                             |
| $N \to c$ | $\mathrm{compose}(M_v, M_c)$ at $M_c$ positions ($H' \cup D$, where $D$ are the child's own positions) |

The net substitution change is $\Delta = -\lvert R\rvert$: one mutation removed per reverted position, none added. The two relocated edges keep their keys via [`Graph::reparent_edge()`](../../packages/treetime-graph/src/graph_ops.rs#L123), so their partition entries stay valid and only their content is rewritten.

Branch lengths split in proportion to substitution count so that root-to-$v$ and root-to-$c$ distances are preserved exactly and re-fit on the next iteration: $b(u \to N) = b(u \to v)\cdot\lvert T\rvert / \lvert M_v\rvert$, $b(N \to v) = b(u \to v) - b(u \to N)$, $b(N \to c) = b(N \to v) + b(v \to c)$. This deliberately differs from the Jukes-Cantor recomputation used by the `prune` merge ([prune-merge-jukes-cantor-branch-length.md](prune-merge-jukes-cantor-branch-length.md)), which would discard the converged branch lengths inside the loop.

Indels use an all-or-nothing rule, because `compose_indels` merges overlapping and adjacent ranges rather than being position-keyed. When no child indel overlaps or is adjacent to a parent indel, the parent indels move cleanly above $N$ and $N \to c$ carries the child's own indels; otherwise the parent indels stay on $N \to v$ and $N \to c$ carries the full composition. Both branches preserve the root-to-$v$ and root-to-$c$ indel content, and the substitution gain is unaffected.

## Why v1 differs

An arbitrary binary resolution of a polytomy can spend two mutations (a substitution plus its reversion) where one suffices, and tree builders produce such resolutions routinely. Fitch's forward pass never emits a reversion on a binary node's child edge, so the pattern is inherently a polytomy phenomenon; that is why the move belongs in the per-polytomy routine rather than as a general edge-level rewrite. Inserting a node beats re-attaching the child directly to the parent, which would duplicate $T$ onto the re-attached edge and give $\Delta = \lvert R\rvert - \lvert T\rvert$, almost never a gain because $\lvert T\rvert$ is typically much larger than $\lvert R\rvert$.

The routine runs every iteration, independent of whether a collapse fired, so reversions present in the input and in polytomies formed by earlier iterations are resolved. It cannot oscillate with the zero-optimal collapse: the lexicographic (total fitch mutation count, node count) potential strictly decreases on every applied merge, hoist, or retirement, even though the hoist deliberately relocates mutation-carrying edges that the collapse step refuses to touch.

## Greedy limitation

When every child of a polytomy reverts the same position, the merge reduces them to a single reverting child and the hoist declines (moving that child would leave the node childless). One residual reversion remains rather than dissolving a pre-existing node ([kb/issues/M-optimize-reversion-hoist-single-child-residual.md](../issues/M-optimize-reversion-hoist-single-child-residual.md)). Disagreeing reverting children likewise leave irreducible homoplasy; exact resolution is subtree parsimony re-optimization, which is out of scope.

## Affected commands

- [`optimize`](../../packages/treetime/src/optimize/run_loop.rs#L319) - runs the routine in its per-iteration topology-cleanup step
- `prune` is unchanged: it still calls the shared-mutation merge directly and does not adopt the reversion hoist

## Tests

- Unit (hoist mechanics): [test_hoist_reversions.rs](../../packages/treetime/src/optimize/topology/__tests__/test_hoist_reversions.rs) - no duplication of $T$, chained composition, reversion removal, distance preservation, per-partition splits, indel cases.
- Integration (merge to hoist to retire): [test_resolve_polytomy.rs](../../packages/treetime/src/optimize/topology/__tests__/test_resolve_polytomy.rs) - the worked example reaching the parsimony optimum, incompatible splits stopping at the greedy bound, helper retirement sparing pre-existing internal nodes, root and reversion-free polytomies left untouched.
- Property: [test_prop_resolve_polytomy.rs](../../packages/treetime/src/optimize/topology/__tests__/test_prop_resolve_polytomy.rs) - the potential never rises and strictly falls on any change, and leaves, non-negative branch lengths, and the single-root tree shape are preserved.
