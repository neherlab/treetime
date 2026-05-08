# Chapter 7: Polytomy resolution

[Back to index](README.md) | Previous: [Chapter 6: Zero-length branches](6-zero-length-branches.md) | Next: [Chapter 8: Initial estimation](8-initial-estimation.md)

## What is a polytomy?

A **polytomy** (multifurcation) is an internal node with three or more children. Most tree-building algorithms produce binary (bifurcating) trees where every internal node has exactly two children. A polytomy means three or more lineages diverge from the same ancestor:

```
Bifurcation (normal):     Polytomy (3-way):       Polytomy (5-way):

      P                        P                        P
     / \                     / | \                   / / | \ \
    A   B                   A  B  C                 A B  C  D  E
```

### Soft vs hard polytomies

A **soft polytomy** results from insufficient data. The lineages did diverge at different times, but the alignment lacks mutations to determine the order. With more data, a soft polytomy would resolve into binary splits. Most polytomies in molecular phylogenetics are soft.

A **hard polytomy** represents genuine simultaneous divergence -- three or more lineages arose at the same time. This can happen during rapid radiation events. Hard polytomies are rare.

The distinction matters: soft polytomies should be resolved (the true tree is binary), while hard polytomies should be preserved.

### Why tree builders produce polytomies

Tree builders produce fully bifurcating trees. When the data does not support a particular split, the builder picks one of many possible binary resolutions and assigns near-zero length to the unsupported internal branch. The number of possible binary resolutions grows rapidly:

| Polytomy degree | Binary resolutions |
| --------------- | ------------------ |
| 3               | 3                  |
| 4               | 15                 |
| 5               | 105                |
| 10              | 34,459,425         |

The tree builder's choice among these is arbitrary. Pruning the near-zero branches (Chapter 6) reveals the underlying polytomy. This chapter covers the reverse: restructuring a polytomy into supported binary splits.

### Downstream consequences

Arbitrary binary resolutions cause: inflated apparent tree resolution, confounded ancestral reconstruction (mutations assigned to wrong branches), and numerical issues (zero-length branches in rate estimation). See [Chapter 1](1-introduction.md) for the full discussion.

## Strategy 1: Neighbor-joining

The neighbor-joining (NJ) algorithm (Saitou and Nei 1987) starts from a star tree (complete polytomy) and iteratively joins the closest pair of taxa, creating binary structure. At each step, NJ selects the pair `(i,j)` minimizing the total tree length -- not just the two closest taxa, but accounting for distances to all other taxa via the Q-criterion:

```
Q(i,j) = (n-2) * d(i,j) - sum_k d(i,k) - sum_k d(j,k)
```

where n is the number of taxa and d is the distance matrix.

NJ is guaranteed to recover the correct tree from additive distances (Studier and Keppler 1988). Atteson (1999) proved the **l-infinity safety radius**: NJ recovers the correct topology when the maximum entry-wise distance error is less than half the shortest edge length. This bound is tight -- no distance method can do better.

NJ runs in O(n^3) time and is a natural fit for polytomy resolution because it starts from a star topology. Frederick Matsen proposed NJ-based resolution in [treetime issue #109](https://github.com/neherlab/treetime/issues/109) as an O(n^1.5) alternative. Not implemented in v0 or v1.

## Strategy 2: Greedy likelihood-based merging (temporal)

This is what TreeTime uses in the `timetree` command. For each polytomy, the algorithm evaluates every pair of children and computes the likelihood gain from merging them under a new internal node.

Given a polytomy node P with children, merging A and B creates a new node N:

```
Before:              After:

      P                    P
    / | \                / | \
   A  B  C             N  C  ...
                       / \
                      A   B
```

The **cost gain** of merging is:

```
cost_gain = [old_LH(P-A) + old_LH(P-B)] - [new_LH(N-A) + new_LH(N-B)] - penalty(P-N)
```

The penalty `mu * L * dt` discourages long zero-mutation branches: longer branches, higher mutation rates, and longer alignments make zero-mutation branches less probable. The optimal time placement of N is found by Brent's method.

The algorithm greedily picks the best pair above the **resolution threshold** (default 0.05), creates the new node, updates the pairwise gain matrix, and repeats until no pair exceeds the threshold.

**Stretched vs compressed:** v0 separates children into "stretched" (`mutation_length < clock_length` -- fewer mutations than the clock predicts, candidates for artifacts) and "compressed" (at least as many mutations as predicted -- genuine divergence). By default only stretched children are merged. v1 does not implement this distinction.

**Complexity:** O(n^2) pairwise comparisons per polytomy. Acceptable for small polytomies (3-10 children) but impractical for large ones.

v0 code: `_poly()` at [`packages/legacy/treetime/treetime/treetime.py#L713-L870`](../../../packages/legacy/treetime/treetime/treetime.py#L713-L870).

v1 code: `resolve_polytomies()` at [`packages/treetime/src/commands/timetree/optimization/polytomy.rs`](../../../packages/treetime/src/commands/timetree/optimization/polytomy.rs).

## Strategy 3: Stochastic coalescent resolution

For large polytomies, greedy merging produces "caterpillar-like" subtrees (long chains of binary nodes). The stochastic approach simulates a backward-in-time coalescent process (Kingman 1982):

1. Sort children by sampling time (most recent first).
2. Count mutations per child: `n_muts = round(mutation_length * L)`.
3. Simulate backward in time. At each step:
   - Branches with zero remaining mutations are "ready to coalesce."
   - A random event occurs: either a mutation is removed from a branch (making it closer to ready) or two ready branches merge under a new node.
   - Coalescence rate comes from the population model or a dummy rate.
4. Remaining uncoalesced branches stay as direct children of the parent.

The process is stochastic -- different runs produce different topologies. This reflects genuine uncertainty about the branching order within the polytomy.

**Complexity:** O(n) per polytomy.

v0 code: `generate_subtree()` at [`packages/legacy/treetime/treetime/treetime.py#L872-L1010`](../../../packages/legacy/treetime/treetime/treetime.py#L872-L1010). v0 has deprecated greedy mode with a warning recommending stochastic resolution.

v1: not implemented. Tracked as `N-timetree-stochastic-polytomy-unimplemented`.

## Strategy 4: Shared-mutation merging (sequence-based)

This strategy does not use time or likelihood. It examines the substitution sets on child branches and groups siblings sharing identical mutations:

```
Before (A and B share G100T):      After merging:

      P                                  P
    / | \                               / \
   A  B  C                            N   C
   |  |  |                           / \
  G100T  G100T  A200C               A   B
  A200C  T300G
                                    Edge P-N: G100T (shared)
                                    Edge N-A: A200C (unique)
                                    Edge N-B: T300G (unique)
```

New branch length = `#shared_mutations / alignment_length`.

The algorithm finds all polytomy nodes, computes substitution intersections for all child pairs, greedily merges the pair with the most shared mutations, and repeats until no pair shares mutations.

**Complexity:** O(n^2 \* m) where n is polytomy degree and m is mutations per branch.

v0: no formal implementation. The design document describes "ad-hoc scripts" in nextstrain pipelines. See the [nextstrain ad-hoc scripts](#the-nextstrain-ad-hoc-scripts) section below for the script locations.

v1 code: `merge_shared_mutation_branches()` in [`packages/treetime/src/representation/algo/topology_cleanup/merge_shared_mutations.rs`](../../../packages/treetime/src/representation/algo/topology_cleanup/merge_shared_mutations.rs). Available as `--merge-shared-mutations` on the `prune` command and in the `optimize` topology-cleanup loop.

## Comparison

| Aspect           | NJ                   | Greedy temporal          | Stochastic coalescent   | Shared-mutation           |
| ---------------- | -------------------- | ------------------------ | ----------------------- | ------------------------- |
| Input data       | Distance matrix      | Time distributions       | Mutation counts + times | Substitution sets         |
| Criterion        | Minimize tree length | Maximize likelihood gain | Simulate coalescent     | Maximize shared mutations |
| Complexity       | O(n^3)               | O(n^2)                   | O(n)                    | O(n^2 \* m)               |
| Deterministic    | Yes                  | Yes                      | No                      | Yes                       |
| Large polytomies | Good                 | Bad (caterpillar)        | Good                    | Acceptable                |
| Requires time    | No                   | Yes                      | Yes                     | No                        |
| In v0            | No                   | Yes (deprecated)         | Yes (recommended)       | No                        |
| In v1            | No                   | Yes (timetree)           | No                      | Yes (prune only)          |

## The star tree paradox

Steel and Penny (2000) showed that ML always prefers a resolved tree over the true star tree given enough data. This means any ML-based resolution strategy must require a **minimum likelihood gain** (the resolution threshold) before accepting a merge. The default 0.05 in TreeTime is a heuristic balance between over-resolving and under-resolving.

Lewis, Holder, and Holsinger (2005) showed Bayesian inference has the same bias. They proposed placing prior weight on polytomous topologies.

## The nextstrain ad-hoc scripts

Production pipelines use ad-hoc scripts for polytomy handling:

- [**`sanitize_trees.py`**](https://github.com/nextstrain/seasonal-flu/blob/511e3661aaef0ca491648514f7afc9f04aeae824/scripts/sanitize_trees.py) (seasonal-flu): calls `prune_short_branches()` to collapse near-zero branches.
- [**`find_clusters.py`**](https://github.com/nextstrain/ncov/blob/17fc7a64911baf7ef81b1b968b32537c441153a8/scripts/find_clusters.py#L28) (ncov): collapses branches < 1e-5, identifies identical-sequence clusters at polytomies.

These perform the "prune" half. The "resolve" half is what `--merge-shared-mutations` and the timetree polytomy resolution provide.

## References

- Saitou, N., and M. Nei. 1987. "The Neighbor-Joining Method." _Mol. Biol. Evol._ 4(4):406-425. https://doi.org/10.1093/oxfordjournals.molbev.a040454
- Studier, J. A., and K. J. Keppler. 1988. "A Note on the Neighbor-Joining Algorithm." _Mol. Biol. Evol._ 5(6):729-731. https://doi.org/10.1093/oxfordjournals.molbev.a040527
- Atteson, K. 1999. "The Performance of Neighbor-Joining Methods." _Algorithmica_ 25(2-3):251-278. https://doi.org/10.1007/pl00008277
- Kingman, J. F. C. 1982. "On the Genealogy of Large Populations." _J. Appl. Prob._ 19(A):27-43. https://doi.org/10.2307/3213548
- Fitch, W. M. 1971. "Toward Defining the Course of Evolution." _Syst. Zool._ 20(4):406-416. https://doi.org/10.2307/2412116
- Steel, M., and D. Penny. 2000. "Parsimony, Likelihood, and the Role of Models." _Mol. Biol. Evol._ 17(6):839-850. https://doi.org/10.1093/oxfordjournals.molbev.a026364
- Lewis, P. O., M. T. Holder, and K. E. Holsinger. 2005. "Polytomies and Bayesian Phylogenetic Inference." _Syst. Biol._ 54(2):241-253. https://doi.org/10.1080/10635150590924208
- Buneman, P. 1974. "A Note on the Metric Properties of Trees." _J. Combin. Theory B_ 17(1):48-50. https://doi.org/10.1016/0095-8956(74)90047-1
- Sagulenko, P., V. Puller, and R. A. Neher. 2018. "TreeTime." _Virus Evol._ 4(1):vex042. https://doi.org/10.1093/ve/vex042
