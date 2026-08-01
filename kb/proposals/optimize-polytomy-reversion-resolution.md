# Proposal: reversion-driven polytomy resolution in the optimize loop

## Summary

Extend the optimize loop's topology cleanup with a third local move: when a child edge of
node $v$ carries the reversion of a substitution mapped to $v$'s parent edge, insert a new
node above $v$ that groups $v$ and that child, moving the non-reverted mutations above the
new node. The move removes exactly one mutation per reverted position and never adds any.

Combined with the existing shared-mutation merge and a helper-node cleanup step, this turns
per-polytomy topology cleanup into a single three-step routine: **merge shared mutations →
hoist reverting child → retire helper nodes**.

## Current state

`fn run_optimize_loop()` performs topology cleanup once per iteration through
`fn prune_and_merge_in_loop()`
[`packages/treetime/src/optimize/run_loop.rs#L309-L353`](../../packages/treetime/src/optimize/run_loop.rs#L309-L353):

1. `fn find_zero_optimal_internal_edges()` collects internal edges the per-edge optimizer
   drove to exactly zero, excluding edges that carry substitutions or indels.
2. `fn collapse_edge()` contracts each such edge, producing polytomies.
3. `fn merge_shared_mutation_branches()` groups siblings in those polytomies that share
   identical substitutions under a new internal node.

Step 3 is gated on step 2 having fired: `prune_and_merge_in_loop` returns early when the
zero-optimal edge list is empty.

Nothing in the loop acts on reversions. A reversion is currently only ever removed as a side
effect of `fn compose_substitutions()`
[`packages/treetime/src/seq/mutation.rs#L176-L227`](../../packages/treetime/src/seq/mutation.rs#L176-L227)
when an edge happens to be collapsed for an unrelated reason.

## Motivation

An arbitrary binary resolution of a polytomy can force a substitution onto an internal edge
and then require its reversion on one child. The pattern costs two mutations where a
different resolution costs one. Tree builders produce such resolutions routinely, and the
optimize loop already creates fresh polytomies by collapsing zero-length edges, so the
pattern arises both in the input and during the run.

Fitch's forward pass resolves an internal node toward the parent state when its children
conflict, so a *binary* node never emerges from reconstruction with a reversion on one child
edge — the reconstruction would place a single mutation on the other child instead. The
pattern is therefore inherently a **polytomy** phenomenon, which is why the move belongs in
the per-polytomy routine rather than as a general edge-level rewrite.

## Move definition

Let $v$ be an internal non-root node, $u$ its parent, $e_p = u \to v$ carrying substitutions
$M_v$, and $e_c = v \to c$ a child edge carrying $M_c$. Because every edge holds at most one
substitution per position (asserted in `fn compose_substitutions()`), $M_v$ partitions by
position against $M_c$ into three disjoint sets:

- $R$ — positions where $M_c$ holds the exact inverse (**reversions**)
- $H$ — positions where $M_c$ holds a different substitution (**chains**: $a \to b$ then $b \to d$)
- $T$ — the remainder of $M_v$, untouched by the child

Let $D$ be the positions of $M_c$ absent from $M_v$, and let $H'$ denote the composed
chains $\{a \to d\}$ for $p \in H$.

The move inserts a new node $N$ between $u$ and $v$:

| edge | substitutions | count |
| --- | --- | --- |
| $u \to N$ | $T$ | $\lvert T\rvert$ |
| $N \to v$ | $H \cup R$ | $\lvert H\rvert + \lvert R\rvert$ |
| $N \to c$ | $H' \cup D$ | $\lvert H\rvert + \lvert D\rvert$ |

Before the move the two edges carry $(\lvert T\rvert + \lvert H\rvert + \lvert R\rvert) +
(\lvert R\rvert + \lvert H\rvert + \lvert D\rvert)$. After, they carry
$\lvert T\rvert + 2\lvert H\rvert + \lvert R\rvert + \lvert D\rvert$.

$$\Delta = -\lvert R\rvert$$

**Trigger: $R \neq \emptyset$. Gain: $\lvert R\rvert$, independent of $\lvert T\rvert$.**

$N \to c$ is exactly $\mathrm{compose}(M_v, M_c) \setminus T$, so it can be produced with the
existing `fn chain_fitch_subs()` rather than new composition logic.

### Why not re-attach the child directly to the parent

The simpler move — detach $c$ and make it a child of $u$ with
$\mathrm{compose}(M_v, M_c) = T \cup H' \cup D$ — gives
$\Delta = \lvert R\rvert - \lvert T\rvert$, because it duplicates $T$ onto the re-attached
edge. On real data $\lvert T\rvert$ is typically much larger than $\lvert R\rvert$, so that
form would almost never fire. It is rejected.

`fn merge_shared_mutation_branches()` would in principle recover the node-insertion state
from the re-attached state on a later iteration: the shared set between $u \to v = M_v$ and
$u \to c = T \cup H' \cup D$ is precisely $T$. But it cannot when $u$ is left with two
children — `fn find_polytomy_nodes()` requires `degree_out > 2` and
`fn merge_single_polytomy()` breaks at `child_edges.len() <= 2`
[`packages/treetime/src/optimize/topology/merge_shared_mutations.rs#L64-L100`](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs#L64-L100).
That is a wrong fixed point, not a delayed one. Inserting $N$ directly avoids it and also
avoids optimizing branch lengths against the inflated intermediate tree for one iteration.

### Branch lengths

Distance-preserving split, proportional to substitution count:

$$b(u \to N) = b(u \to v)\cdot\frac{\lvert T\rvert}{\lvert M_v\rvert},\quad
b(N \to v) = b(u \to v) - b(u \to N),\quad
b(N \to c) = b(N \to v) + b(v \to c)$$

Root-to-$v$ and root-to-$c$ distances are both preserved exactly; the next iteration re-fits.

This deliberately differs from `fn merge_sibling_group()`, which recomputes branch lengths
from Jukes-Cantor distance
[`packages/treetime/src/optimize/topology/merge_shared_mutations.rs#L312-L331`](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs#L312-L331).
That is appropriate for `prune`, which has no optimizer state, but inside the optimize loop it
would discard converged branch lengths.

### Indels

`fn compose_indels()` merges overlapping and adjacent ranges rather than being position-keyed,
so an exact three-way split is not always definable. Conservative rule: hoist $e_p$'s indels
to $u \to N$ only when the child's indels neither cancel nor overlap them; otherwise leave
them on $N \to v$ and let $N \to c$ carry $\mathrm{compose}(I_v, I_c)$. Always correct, and
leaves at most $\lvert I_v\rvert$ unclaimed. The substitution gain $\lvert R\rvert$ is
unaffected either way.

## Interaction with the shared-mutation merge

Running the merge **first** canonicalizes the input to the hoist. When several siblings revert
the same position, their edges literally share the mutation $b \to a$ at $p$, so the merge
groups them under a helper node. The hoist then sees a single reverting child.

Worked example: $M_v = \{p, q\}$, children $c_1$ and $c_2$ both reverting $p$, $c_3$ keeping it.

| step | edges | total |
| --- | --- | --- |
| start | $u\to v=\{p,q\}$, $v\to c_1=\{\lnot p\}$, $v\to c_2=\{\lnot p\}$, $v\to c_3=\emptyset$ | 4 |
| 1. merge | $u\to v=\{p,q\}$, $v\to N'=\{\lnot p\}$, $N'\to c_1=\emptyset$, $N'\to c_2=\emptyset$ | 3 |
| 2. hoist $N'$ | $u\to N=\{q\}$, $N\to v=\{p\}$, $N\to N'=\emptyset$, $v\to c_3=\emptyset$ | 2 |
| 3. retire $N'$ | $u\to N=\{q\}$, $N\to v=\{p\}$, $N\to c_1=\emptyset$, $N\to c_2=\emptyset$ | 2 |

The parsimony optimum is 2 ($q$ on every lineage, $p$ only on $c_3$'s). The routine reaches it.

Step 3 is not cosmetic. The hoist produces an **empty** $N \to c$ edge exactly when the child
edge consisted of nothing but reversions, which is the normal case after a merge. Collapsing
it is what dissolves the helper node and makes $c_1$ and $c_2$ genuine siblings of $v$. The
merge creates the helper, the hoist relocates it, the collapse retires it.

The merge cannot undo a hoist: after a merge places shared substitution $s$ on the group's
parent edge, $s$ is removed from every group member's edge, and since an edge holds at most
one substitution per position, no member can then carry $s$'s reversion.

## Greedy limits

When reverting children disagree on *which* positions they revert, one hoist is often already
optimal. With $M_v = \{p_1, p_2, q\}$, $c_1$ reverting $p_1$, $c_2$ reverting $p_2$, $c_3$
keeping both, hoisting $c_1$ takes the total from 5 to 4. No tree does better: $p_1$ is needed
by $\{c_2, c_3\}$ and $p_2$ by $\{c_1, c_3\}$, which are incompatible splits, so the residual
reversion on $v \to c_2$ is irreducible homoplasy.

Genuine greedy failures exist. With $c_1$ reverting $\{p_1, p_2\}$ and $c_2$ reverting
$\{p_2, p_3\}$, greedy reaches 5 against an optimum of 4. Exact resolution is subtree
parsimony re-optimization and is out of scope; `fn merge_shared_mutation_branches()` carries
the same greedy caveat via `fn greedy_disjoint_group_matching()`.

Rule adopted: **one reverting child per node per round**, choosing the child with the largest
$\lvert R\rvert$, ties broken by edge key for determinism. Iterate the node to a fixpoint.

## Termination

Potential function: (total fitch mutation count, node count), lexicographic.

- hoist: strictly decreases the first component by $\lvert R\rvert \geq 1$
- merge: strictly decreases the first component by $\text{total\_shared}\cdot(k-1)$
- helper retirement and zero-optimal collapse: first component unchanged, second decreases

No cycling. This has to be stated explicitly because the new move deliberately relocates
mutation-carrying edges, which `fn find_zero_optimal_internal_edges()` documents as previously
forbidden precisely to avoid a merge/collapse oscillation
[`packages/treetime/src/optimize/run_loop.rs#L270-L272`](../../packages/treetime/src/optimize/run_loop.rs#L270-L272).

## Design

### New module

`packages/treetime/src/optimize/topology/resolve_polytomy.rs`

```rust
/// Merge → hoist → retire, at one polytomy. Returns whether anything changed.
fn resolve_one(
  graph: &mut GraphAncestral,
  sparse_partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  dense_partitions: &[Arc<RwLock<PartitionMarginalDense>>],
  node_key: GraphNodeKey,
) -> Result<bool, Report>;

/// Drive `resolve_one` over all polytomies to a fixpoint. Returns moves applied.
pub fn resolve_polytomies(
  graph: &mut GraphAncestral,
  sparse_partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  dense_partitions: &[Arc<RwLock<PartitionMarginalDense>>],
) -> Result<usize, Report>;
```

`packages/treetime/src/optimize/topology/hoist_reversions.rs` holds the $(T, H, R)$ split and
the node insertion. The split is a single merge-walk over two position-sorted `Vec<Sub>`,
the same shape as `fn compose_substitutions()`.

Sparse-only, matching `fn merge_shared_mutation_branches()`: dense partitions carry no
mutation lists. `optimize` builds either a sparse or a dense partition, never both
[`packages/treetime/src/optimize/pipeline.rs#L93-L102`](../../packages/treetime/src/optimize/pipeline.rs#L93-L102),
so the routine is inert under `--dense`.

### Helper retirement scope

Step 3 collapses a mutation-free internal edge **only when its target is a node created
during this invocation of `resolve_one`**. Track created keys in a local set.

The restriction matters: a merge's $N' \to c_i$ edges are frequently mutation-free too, but
their targets are pre-existing subtree roots. Collapsing those would flatten topology the
input asserted and discard its branch lengths.

Within that scope the rule relaxes `fn find_zero_optimal_internal_edges()`, which additionally
requires `bl == 0`. Justified because these edges were synthesized moments earlier by the
routine itself, so there is no optimizer decision to override.

### Partition bookkeeping

Node insertion follows `fn merge_sibling_group()`: `add_node`, `add_edge`, `remove_edge`, then
rewrite the sparse `partition.edges` entries. No `treetime-graph` change is required.

Reparenting via a new `Graph::reparent_edge()` primitive (mutating `set_source`, as
`Graph::collapse_edge()` already does internally) would avoid edge-key churn and preserve
dense edge state. Deferred: `fn merge_sibling_group()` already establishes remove-and-add as
the accepted pattern, and matching it keeps the two routines uniform.

Register the new node and edge keys in dense partitions, mirroring
`PartitionMarginalDense::apply_reroot()`
[`packages/treetime/src/partition/marginal_dense.rs#L96-L120`](../../packages/treetime/src/partition/marginal_dense.rs#L96-L120).
`fn merge_sibling_group()` omits this; it is harmless today only because dense and sparse
partitions never coexist in `optimize`. Not worth adding a second instance of the same latent
trap.

Stale `msg_to_parent` / `msg_to_child` / `msg_from_child` on rewritten edges need no handling:
`fn marginal_backward()` overwrites them wholesale
[`packages/treetime/src/partition/marginal_passes.rs#L386`](../../packages/treetime/src/partition/marginal_passes.rs#L386),
which is what the collapse path already relies on.

### Loop integration

`fn prune_and_merge_in_loop()` becomes:

1. zero-optimal collapse (existing, optimizer-driven — it is what creates polytomies)
2. `resolve_polytomies()` — **every iteration**, no longer gated on step 1 having fired
3. `graph.build()` + `assign_node_names()`

`fn merge_shared_mutation_branches()` is no longer called directly from the loop;
`fn merge_single_polytomy()` is promoted to `pub(crate)` for reuse by `resolve_one`. The
public wrapper stays for `prune`
[`packages/treetime/src/prune/pipeline.rs#L79`](../../packages/treetime/src/prune/pipeline.rs#L79).

Each applied move sets `topology_changed`, which resets `best_lh`. The monotone potential
bounds how many iterations can keep firing.

## Open axes

- **`prune` adoption.** `prune` already calls the merge and would benefit from the same
  routine, but adopting it changes `prune`'s output. Separate decision.
- **Accepted wart.** In the common case the merge creates a helper node that step 3 destroys
  moments later. Fusing the two passes to skip the transient would entangle them for no
  parsimony gain.
- **Per-iteration cost.** The pass is $O(\text{edges} \times \text{mutations per edge})$ per
  round. `Graph::remove_edge()` scans all nodes, so node insertion is $O(V)$ — inherited from
  `fn merge_sibling_group()`. Measure before optimizing.

## Validation plan

Unit tests in `packages/treetime/src/optimize/topology/__tests__/`, reusing the fixture
helpers from `test_collapse_edge.rs` (`fn make_sparse_partition()`,
`fn populate_test_nodes()`, `fn nwk_read_str()`):

- polytomy with large $\lvert T\rvert$ — asserts $T$ is *not* duplicated, distinguishing the
  move from re-attachment to the parent
- chained ($H$) positions composed correctly on both output edges
- the worked merge → hoist → retire example above, end to end, asserting the final mutation
  count is 2
- incompatible-splits example: total goes 5 → 4 and stops
- helper retirement does not dissolve pre-existing internal nodes on mutation-free edges
- indel cancellation, and the overlapping-indel fallback
- root skipped; distance preservation of root-to-node paths
- multi-partition

Loop-level:

- planted reversion converges with the reversion removed, no oscillation, iteration count not
  regressed
- property test in the style of `test_prop_merge_shared_mutations.rs`: the potential
  (mutation count, node count) decreases lexicographically on every applied move

Datasets:

- ebola/20 and sc2/4500: count moves applied, total mutation count before and after, final
  log-likelihood, iterations to convergence
- confirm no change on trees with no reversions

## Related

- [kb/features/optimize.md](../features/optimize.md) — needs a topology section; the existing
  shared-mutation merge is undocumented there
- [kb/proposals/optimize-short-branch-pruning.md](optimize-short-branch-pruning.md) —
  complementary post-loop cleanup on branch length rather than mutation content
- [kb/decisions/prune-merge-jukes-cantor-branch-length.md](../decisions/prune-merge-jukes-cantor-branch-length.md) —
  the JC branch-length convention this proposal declines to reuse inside the loop
- [kb/decisions/timetree-no-zero-branch-collapse-in-loop.md](../decisions/timetree-no-zero-branch-collapse-in-loop.md)
