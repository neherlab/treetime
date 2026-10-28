# Optimize loop resolves reversion-driven polytomies

The `optimize` loop adds a topology move with no v0 counterpart: when a polytomy forces a substitution onto an internal edge and then reverts it on a child, the move relocates the child so the reversion cancels. This removes homoplasy that an arbitrary binary resolution of a polytomy introduces, and that v0's loop never targets.

## What v0 does

v0's branch-length optimization loop (`TreeAnc.optimize_tree_marginal`) collapses short internal branches (`prune_short_branches`) but has no step that inspects or removes reversions. A reversion is only ever dropped as an incidental side effect when an edge is collapsed for an unrelated reason and its substitutions compose with the child's. v0 has no shared-mutation merge and no reversion-aware polytomy handling.

## What v1 does

Once per iteration, [`prune_and_merge_in_loop()`](../../packages/treetime/src/optimize/run_loop.rs#L318) runs a single per-polytomy routine, [`resolve_polytomies()`](../../packages/treetime/src/optimize/topology/resolve_polytomy.rs#L22), after the zero-optimal edge collapse. Sparse partitions only; dense partitions carry no per-edge mutation lists. The routine applies three steps at each polytomy and iterates to a fixpoint:

- Merge siblings that share substitutions under a helper node ([`merge_single_polytomy()`](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs#L58)). This canonicalizes several children reverting the same position into one reverting child. A group that contains every child of the node is not merged, because the node would keep a single child.
- Hoist the reverting child with the largest reversion count ([`hoist_reverting_child()`](../../packages/treetime/src/optimize/topology/hoist_reversions.rs#L134)), ties broken by edge key. When the hoist leaves the node with a single child, the edge into the node is collapsed: its substitutions compose with those of the remaining child edge, so a reversion on that child edge cancels.
- Retire helper nodes: collapse mutation-free edges to nodes the routine created in this pass.

The hoist inserts a new node $N$ between parent $u$ and node $v$, grouping $v$ with one reverting child $c$. The parent-edge substitutions $M_v$ partition by position against the child substitutions $M_c$ into three disjoint sets: $T$ (untouched by the child), $H$ (chained, $a \to b$ then $b \to d$), and $R$ (reverted, $a \to b$ then $b \to a$). The resulting edges carry:

| edge      | substitutions                                                                                          |
| --------- | ------------------------------------------------------------------------------------------------------ |
| $u \to N$ | $T$                                                                                                    |
| $N \to v$ | $H \cup R$                                                                                             |
| $N \to c$ | $\mathrm{compose}(M_v, M_c)$ at $M_c$ positions ($H' \cup D$, where $D$ are the child's own positions) |

The net substitution change is $\Delta = -\lvert R\rvert$: one mutation removed per reverted position, none added. The two relocated edges keep their keys via [`Graph::reparent_edge()`](../../packages/treetime-graph/src/graph_ops.rs#L73), so their partition entries stay valid and only their content is rewritten.

Branch lengths split in proportion to substitution count so that root-to-$v$ and root-to-$c$ distances are preserved exactly and re-fit on the next iteration: $b(u \to N) = b(u \to v)\cdot\lvert T\rvert / \lvert M_v\rvert$, $b(N \to v) = b(u \to v) - b(u \to N)$, $b(N \to c) = b(N \to v) + b(v \to c)$. This deliberately differs from the Jukes-Cantor recomputation used by the `prune` merge ([prune-merge-jukes-cantor-branch-length.md](prune-merge-jukes-cantor-branch-length.md)), which would discard the converged branch lengths inside the loop.

Indels use an all-or-nothing rule, because `compose_indels` merges overlapping and adjacent ranges rather than being position-keyed. When no child indel overlaps or is adjacent to a parent indel, the parent indels move cleanly above $N$ and $N \to c$ carries the child's own indels; otherwise the parent indels stay on $N \to v$ and $N \to c$ carries the full composition. Both branches preserve the root-to-$v$ and root-to-$c$ indel content, and the substitution gain is unaffected.

### Reversions across a bifurcating root

A degree-2 root is a reversible pass-through: its two child edges $R \to v$ and $R \to s$ are the two halves of the single unrooted edge $v\,\text{--}\,s$, and the substitution at a site is free to sit on either half. At a site where $v$ equals the root but the sibling $s$ differs, the distinguishing substitution $\text{root} \to s$ sits on the sibling edge, and $v$'s own parent edge is empty there. The per-node reversion scan then sees nothing to revert on $v$'s parent edge, so a child of $v$ that makes the same change is missed, and the tree keeps a substitution and its reversion that a different rooting would avoid.

[`try_hoist_reverting_child()`](../../packages/treetime/src/optimize/topology/resolve_polytomy.rs#L102) closes this gap for a node directly under a bifurcating root. Reversion scoring augments the parent edge with the inverted sibling substitutions at parent-empty sites ([`count_child_reversions()`](../../packages/treetime/src/optimize/topology/hoist_reversions.rs#L20)), so a cross-root reversion is counted when picking the best child. Before hoisting the winning child, [`slide_bifurcating_root_for_child()`](../../packages/treetime/src/optimize/topology/hoist_reversions.rs#L60) re-roots exactly the reverted sites onto the sibling's state: it sets the root character to $s$, drops the sibling substitution, and adds its inverse to $v$'s parent edge, moving the root clamp ($\texttt{root\_sequence}$ and the root node MAP) in lock-step so the next iteration's marginal reconstruction stays consistent. The reversion is then local to $v$'s parent edge and the standard hoist above removes it unchanged. The re-rooting is substitution-count neutral on its own, so it is applied only for the positions the immediately following hoist consumes, preserving the strict decrease of the potential. The net effect is that the routine reaches the same mutation count regardless of where a bifurcating root is placed, matching the result obtained when the same reversion sits on a genuine internal edge.

This applies only to a degree-2 root. A root with three or more children is a genuine unrooted vertex whose edges are real tree edges, and every non-root internal node has a real parent edge, so the artifact does not arise elsewhere.

### Mutations that the steps read

The steps decide only from mutations that enter the likelihood:

- Substitutions: the merge index and the reversion count ([`count_child_reversions()`](../../packages/treetime/src/optimize/topology/hoist_reversions.rs#L20)) read the Fitch substitutions of an edge that its marginal estimate of the current iteration also contains ([`SparseReconstruction::supported_fitch_subs()`](../../packages/treetime/src/partition/marginal/reconstruction.rs#L352)). Fitch assigns an arbitrary state at an ambiguous site, so a Fitch substitution that the marginal reconstruction does not place on the edge is no evidence for a shared ancestor. An edge without a marginal estimate, as in the `prune` command, keeps all its Fitch substitutions; collapse, merge, and hoist compose or split the marginal estimates of the edges they rewrite.
- Indels: the merge index includes indels and helper retirement keeps an edge with indels only under `IndelPolicy::Include`. Under `--no-indels` (`IndelPolicy::Exclude`) indels are not part of the likelihood, so an indel-only merge would create an edge with no likelihood support. The moves still compose and split indels, so the sparse data stay consistent.

### Merges that the likelihood rejects

A merge nests the old tree: with the new edge at length zero, the merged tree has the likelihood of the unmerged one. The merge helps only when the branch-length optimizer gives the new edge a positive length. A substitution that is the most likely one on each sibling edge does not make the grouped edge likely, so the next branch-length step can drive it to zero, the zero-optimal collapse removes it, and the merge would be proposed again in the following iteration.

`struct MergeLedger` (in [`merge_shared_mutations.rs`](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs#L105)) holds, for one optimize run, the set of child nodes of every helper node that a merge created. When the collapse removes the edge into such a helper, the ledger records its child set as rejected before polytomy resolution runs, and `merge_single_polytomy()` skips every later group with the same children, also under the grandparent after the parent edge collapsed. The `prune` command merges with an empty ledger.

## Why v1 differs

An arbitrary binary resolution of a polytomy can spend two mutations (a substitution plus its reversion) where one suffices, and tree builders produce such resolutions routinely. Fitch's forward pass never emits a reversion on a binary node's child edge, so the pattern is inherently a polytomy phenomenon; that is why the move belongs in the per-polytomy routine rather than as a general edge-level rewrite. Inserting a node beats re-attaching the child directly to the parent, which would duplicate $T$ onto the re-attached edge and give $\Delta = \lvert R\rvert - \lvert T\rvert$, almost never a gain because $\lvert T\rvert$ is typically much larger than $\lvert R\rvert$.

The routine runs every iteration, independent of whether a collapse fired, so reversions present in the input and in polytomies formed by earlier iterations are resolved. Within one pass it reaches a fixpoint: the lexicographic (total Fitch mutation count, node count) potential strictly decreases on every applied merge, hoist, retirement, or single-child collapse. Across iterations the routine alternates with the zero-optimal collapse, which removes every internal edge of length zero, also edges that the routine created. A merge whose new edge collapses would then repeat in every iteration; the merge ledger above stops this cycle. The ledger tracks merges only.

## Greedy limitation

When every child of a polytomy reverts the same position, the merge does not group all of them, and the hoist moves the reverting children out one at a time. When one child is left, the collapse of the edge into the node cancels its reversion, so this case reaches the parsimony optimum; the collapse removes the input node, which has no split of its own once it has a single child. Disagreeing reverting children leave irreducible homoplasy; exact resolution is subtree parsimony re-optimization, which is out of scope.

## Affected commands

- [`optimize`](../../packages/treetime/src/optimize/run_loop.rs#L318) - runs the routine in its per-iteration topology-cleanup step
- `prune` calls the shared-mutation merge directly, with `IndelPolicy::Include` and an empty ledger, and does not adopt the reversion hoist. It shares the rule that no group of all children is merged

## Tests

- Unit (hoist mechanics): [test_hoist_reversions.rs](../../packages/treetime/src/optimize/topology/__tests__/test_hoist_reversions.rs) - no duplication of $T$, chained composition, reversion removal, distance preservation, per-partition splits, indel cases, and the root slide (sibling substitution moved onto the parent edge and root clamp, reversion removed by the following hoist).
- Unit (merge evidence): [test_merge_shared_mutations.rs](../../packages/treetime/src/optimize/topology/__tests__/test_merge_shared_mutations.rs) - a Fitch substitution counts only where the marginal estimate has it.
- Unit (step rules): [test_resolve_polytomy_rules.rs](../../packages/treetime/src/optimize/topology/__tests__/test_resolve_polytomy_rules.rs) - indel-only merges and helper retirement under both indel policies, no merge of all children, no single-child node after a hoist.
- Unit (ledger): [test_topology_cleanup_merge_ledger.rs](../../packages/treetime/src/optimize/__tests__/test_topology_cleanup_merge_ledger.rs) - a rejected group is skipped, also under the grandparent after the parent collapsed.
- Sparse against dense: [test_topology_cleanup_sparse_dense.rs](../../packages/treetime/src/optimize/__tests__/test_topology_cleanup_sparse_dense.rs) - root log-likelihood and edge objectives on a gapped alignment after a merge and after a hoist.
- Integration (merge to hoist to retire): [test_resolve_polytomy.rs](../../packages/treetime/src/optimize/topology/__tests__/test_resolve_polytomy.rs) - the worked example reaching the parsimony optimum, the bifurcating-root worked example reaching the one-mutation bipartition cost, incompatible splits stopping at the greedy bound, helper retirement sparing pre-existing internal nodes, root and reversion-free polytomies left untouched.
- Property: [test_prop_resolve_polytomy.rs](../../packages/treetime/src/optimize/topology/__tests__/test_prop_resolve_polytomy.rs) - the potential never rises and strictly falls on any change, leaves, non-negative branch lengths, and the single-root tree shape are preserved, and a bifurcating-root site settles at its reroot-invariant bipartition cost regardless of how many children carry the substitution.
