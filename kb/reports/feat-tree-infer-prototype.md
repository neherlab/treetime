# Prototype of a parsimony tree builder

This report records what a working prototype of a native tree builder showed about the design in [feat-tree-infer.md](feat-tree-infer.md). The prototype follows the recommended combination of that report: parsimony placement on a mutation-annotated tree (MAT), parsimony SPR, then the existing TreeTime pipeline. It is a research prototype in `packages/treetime/src/build/` with an evaluation program `packages/treetime/examples/build_tree.rs`; it has no CLI command and approves no design.

## Executive summary

- **Exact parsimony placement is cheap at the scale of the example data.** A branch-and-bound search places all 5100 SARS-CoV-2 genomes in 1.1 s and visits 14 to 34 times fewer nodes than a full scan, with the same result. See [Placement search](#placement-search)
- **A greedy walk from local-branching-index (LBI) peaks is not reliable on its own.** It matched the exact search on mpox and on perfect phylogenies, but placed 2 to 85% of samples worse on SARS-CoV-2, RSV, dengue, Lassa, and H3N2, by up to 1900 substitutions for a single sample. As a seed of the exact search it saves little. See [Centroid search](#centroid-search)
- **Reusing the Fitch pass is the largest single improvement.** Refreshing the MAT with `fn create_fitch_partition()` between SPR rounds lowers the score by 0.6 to 7%; SPR on the MAT alone finds few moves. The built trees end 0.0 to 4.3% above the Fitch parsimony of the reference trees in `data/`. See [Tree quality](#tree-quality)
- **The output runs through `clock` and `timetree` without changes, but the dating differs.** On SARS-CoV-2 the least-squares root and the clock rate depend on the order of insertion, while the parsimony score barely changes. A builder needs a dating check in its validation, not only parsimony and split agreement. See [Dating with built trees](#dating-with-built-trees)

## What the prototype does

| Step | Implementation | Reused code |
| --- | --- | --- |
| Reference sequence | Per-site most frequent canonical state of the alignment | |
| Query encoding | Differences from the reference at canonical sites; gaps, `N`, and ambiguity codes are missing data | `fn find_letter_ranges_by()` |
| MAT | `Graph` plus a sorted `Vec<Sub>` per edge, the sample of each leaf, the parsimony score, and an upper bound on the substitutions below each node | `Graph`, `Sub`, `fn compose_substitutions()` |
| Start | The root (the reference) with the first sample as its only child. Later samples split that edge, so the first star tree forms by placement | |
| Placement | Score of node $v$: $c(v) = \lvert D_q(v) \rvert$, the differences of query $q$ from the state of $v$. Moving down edge $e$ changes it by $\sum_{s \in e} [q_s \ne s.\mathrm{qry}] - [q_s \ne s.\mathrm{ref}]$. Splitting $e$ after the substitutions that $q$ shares costs $c(\mathrm{parent}) - \mathrm{shared}(e)$ | |
| Batches | A batch of samples is placed in parallel against a frozen MAT. A sample whose target was changed earlier in the batch is placed again | |
| SPR | Detach each subtree in postorder, collapse a parent left with one child, reinsert at the best placement; ties keep the original place | |
| Refresh | Rebuild the root sequence and all edge substitutions from Fitch on the current topology | `fn create_fitch_partition()` |
| Output | Newick with branch lengths of substitutions per site | `NodeNames`, `fn remove_stem_root()` |
| Evaluation | Fitch parsimony of any topology; informative splits (internal edges with at least one Fitch substitution) shared with a reference tree | `fn create_fitch_partition()`, `fn graph_clade_keys()` |

Tests check the placement against an analytic limit: on random trees with infinite-sites mutations, every search mode recovers the true informative splits, and the parsimony score equals the number of variable mutations.

## Placement search

The exact search is a depth-first branch-and-bound. Each substitution on a path changes the score by at most one, so a node $v$ with score $c(v)$ bounds every placement in its subtree by $c(v) - d(v)$, where $d(v)$ is the largest number of substitutions on a path from $v$ to a leaf. The MAT keeps $d(v)$ as an upper bound: insertion raises it along the path to the root, and detaching leaves it unchanged, which keeps the bound valid. Children are visited in order of their bound, so a good placement is found early.

| Dataset | Samples | Full scan, node visits | Branch-and-bound, node visits | Ratio |
| --- | --- | --- | --- | --- |
| `sc2/4500` | 5100 | 19,232,790 | 779,943 | 25 |
| `rsv/a/2000` | 1999 | 3,554,502 | 255,443 | 14 |
| `mpox/clade-ii/2000` | 1989 | 2,573,712 | 155,266 | 17 |
| `dengue/2000` | 1954 | 3,470,503 | 100,696 | 34 |

The full scan of each placement is the oracle in the tests and in the `compare` option of the example program; the branch-and-bound never differed from it.

## Centroid search

The centroid search computes the LBI on the MAT (branch length = number of substitutions + 0.5, $\tau = 1$), takes the highest local maxima as start nodes together with the root, and walks best-first to neighbors whose score is at most the best score plus a slack. The table compares its placements with the full scan on the same tree (placement only, batch size 64).

| Dataset | Centroids | Slack | Samples placed worse | Total excess | Largest excess | Node visits |
| --- | --- | --- | --- | --- | --- | --- |
| `sc2/4500` | 8 | 1 | 877 of 5099 | 15,873 | 117 | 34,838 |
| `sc2/4500` | 32 | 6 | 110 of 5099 | 963 | 96 | 241,643 |
| `rsv/a/2000` | 8 | 1 | 1193 of 1998 | 143,300 | 806 | 9931 |
| `rsv/a/2000` | 32 | 6 | 275 of 1998 | 24,811 | 595 | 18,373 |
| `mpox/clade-ii/2000` | 8 | 1 | 16 of 1988 | 53 | 9 | 73,934 |
| `tb/149` | 8 | 1 | 0 of 148 | 0 | 0 | 749 |

The greedy walk is exact on tree-like data: the score decreases monotonically along the path to the best placement, and the perfect-phylogeny tests confirm it. Homoplasy and missing data create local minima, and the errors compound: a misplaced sample adds substitutions that attract later samples to the same wrong place. Seeding the branch-and-bound with the centroid placement reduced the visits on `sc2/4500` from 779,943 to 548,796 and was slower in wall time, because each start needs the full state of its node.

## Search from reference nodes

The selection of reference nodes is separate from the search that uses them, and neither depends on the root of the MAT. Two selections exist: the LBI peaks with the time scale $\tau = r / \ln n$, where $r$ is the radius of the tree in substitutions and $n$ the number of samples, and the node of least eccentricity in each component left after cutting the $k$ longest internal branches.

The branch-and-bound search starts at the reference nearest to the query and expands best-first in both directions. A move down to child $u$ is bounded by $c(u)$ minus the longest substitution path below $u$; a move up from $v$ to its parent $p$ by $c(p)$ minus the longest path from $p$ that avoids the subtree of $v$. Both lengths come from one pass over the frozen tree per batch. Every part of the tree is visited or pruned by a valid bound, so the result equals the exhaustive search from any start: the start changes how soon a tight upper bound is found, not the answer, and local minima cannot trap the search.

Placement only, batches of 64; the placement score equals that of the exhaustive search in every run:

| Dataset | Root start: node visits / time (s) | Nearest LBI peak: node visits / time (s) | Range over four reference selections, node visits |
| --- | --- | --- | --- |
| `lassa/L/500` | 21,132 / 0.23 | 18,298 / 0.24 | 18,196 to 18,298 |
| `dengue/2000` | 103,312 / 0.31 | 87,062 / 0.44 | 86,771 to 87,124 |
| `rsv/a/2000` | 251,909 / 0.23 | 199,660 / 0.34 | 199,336 to 199,891 |
| `mpox/clade-ii/2000` | 157,908 / 3.16 | 146,825 / 3.28 | 146,825 to 147,005 |
| `sc2/4500` | 781,762 / 0.77 | 549,692 / 1.53 | 549,692 to 549,925 |

The four selections are LBI peaks with the default and with $\tau = 20$, and 16 and 64 cut branches. Observations:

- The reference start saves 7 to 30% of the node visits, and the saving hardly depends on which references are chosen
- The visits that remain are the nodes whose lower bound does not exceed the optimal cost. That set depends on the query and the bound, not on the start, so the tightness of the bound sets the cost of the search. The current bound counts every substitution below a node, although the query can only gain at sites where it differs from that node
- Wall time rose, because the bounds in both directions are recomputed for every batch and the nearest reference is found by comparing the query with every reference
- An earlier version assigned each node to its nearest reference and accepted a placement inside that cell by the triangle inequality. It was exact as well, but slower, and certified only 1 to 68% of the queries

## Tree quality

Settings: branch-and-bound placement, samples with the fewest missing sites first, batches of 64, exhaustive SPR, and a Fitch refresh before each SPR round and at the end, until a round brings no gain. Fitch parsimony is computed on both trees by the same function, with gaps as indels. The reference trees are the `tree.nwk` files of `data/`, which come from maximum-likelihood builders.

| Dataset | Samples | Sites | Build time (s) | Peak memory | Fitch, built | Fitch, reference | Difference | Shared splits | Only built | Only reference |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `tb/149` | 149 | 216 | 0.01 | 6 MB | 216 | 216 | 0.0% | 44 | 3 | 4 |
| `flu/h3n2/500` | 476 | 1407 | 0.08 | 11 MB | 1515 | 1507 | +0.5% | 188 | 15 | 33 |
| `lassa/L/500` | 414 | 10,402 | 2.6 | 63 MB | 78,478 | 75,263 | +4.3% | 257 | 149 | 151 |
| `rsv/a/2000` | 1999 | 15,225 | 3.5 | 137 MB | 42,470 | 42,268 | +0.5% | 1338 | 166 | 246 |
| `mpox/clade-ii/2000` | 1989 | 197,209 | 18.9 | 1.2 GB | 5599 | 5562 | +0.7% | 537 | 87 | 196 |
| `dengue/2000` | 1954 | 10,723 | 5.6 | 146 MB | 55,673 | 55,018 | +1.2% | 1369 | 229 | 281 |
| `sc2/4500` | 5100 | 29,903 | 11.1 | 504 MB | 26,342 | 26,175 | +0.6% | 1969 | 417 | 516 |

Build times are wall time on 12 cores without reading the alignment; peak memory includes the alignment. MAT score after each stage:

| Dataset | Placement | + SPR | + SPR and Fitch refresh | Time: placement / + SPR / + refresh (s) |
| --- | --- | --- | --- | --- |
| `sc2/4500` | 26,730 | 26,693 | 26,488 | 1.1 / 3.1 / 13.6 |
| `mpox/clade-ii/2000` | 5666 | 5638 | 5605 | 3.9 / 4.1 / 21.3 |
| `rsv/a/2000` | 43,923 | 43,846 | 42,678 | 0.2 / 1.4 / 3.5 |
| `dengue/2000` | 59,076 | 58,876 | 55,904 | 0.3 / 2.5 / 5.4 |

Observations:

- SPR on the MAT alone gains 0.1 to 0.5%. Each subtree is moved with the fixed state of its root, so moves that need a different state at the subtree root are not found; matOptimize uses the Fitch state sets of the subtree root instead
- The Fitch refresh brings most of the gain, because placement fixes the internal states greedily and the reference sequence is the consensus, not an optimized root
- The refresh runs Fitch on full sequences of every node, which dominates time and memory on mpox (197 kb genomes). The MAT itself is small; [kb/proposals/mutation-first-sequence-representation.md](../proposals/mutation-first-sequence-representation.md) describes the representation that would remove this cost
- The MAT treats ambiguity codes as missing data and Fitch treats them as state sets, so an SPR move that gains one substitution on the MAT can be undone by the next refresh. The loop stops when a round brings no net gain
- The remaining excess over the reference trees, and the lower number of informative splits than the reference trees, point to missing moves of sibling groups that share substitutions. The `optimize` command's shared-mutation merge (`fn merge_shared_mutation_branches()`) does such moves on the sparse reconstruction

## Why the maximum-likelihood trees were more parsimonious

Two diagnostics on the real datasets:

- **The builder's SPR started from the reference tree** (`start_from_reference`) improves the Fitch score by at most 12 substitutions (`sc2/4500`: 26,175 to 26,163; `flu/h3n2/500` and `lassa/L/500`: no change). The reference trees are local optima of the builder's moves, and the built trees are worse local optima: the difference comes from the search, not from the objective
- **Per site** (`sites`), the extra substitutions of the built trees sit at sites that already change often: the median minor-allele count of these sites is 36 against 2 for all variable sites on `sc2/4500`, and 124 against 8 on `rsv/a/2000`. On `lassa/L/500` the excess is spread over 2207 sites with 1 to 3 substitutions each, the signature of a different arrangement deep in the tree

The two moves of the `optimize` command, merging siblings that share a substitution and hoisting a child that reverts its parent edge, are single-subtree moves that the exhaustive SPR already finds. What single-subtree SPR cannot do is to change several subtrees and their internal states at once.

## Local rebuilds

For every internal node, the ball of radius $r$ branches is cut out of the MAT. Its uppermost node becomes the local root; the nodes at distance $r$ and the samples inside the ball become pseudo-leaves with their current states; the interior is rebuilt by exhaustive placement in 4 orders, each followed by SPR, and replaces the old interior when it needs fewer substitutions. Balls with more than 24 pseudo-leaves (48 for $r = 5$) are skipped.

### Simulated data with a known tree

The data of `validation/` (v0 simulator: Kingman coalescent, random GTR without rate variation among sites, 10,000 sites, replicates 0 and 1 of the constant-size schedule). `data` has 400 tips, `data_n30` 1200; the messy sets add gaps, ambiguity codes, and runs of N. The divergence from the consensus is 0.6, 2.0, 5.7, and 15.6% for $\mu$ = 0.0003, 0.001, 0.003, and 0.01. Recall and precision are those of the informative splits of the true tree (internal edges with at least one Fitch substitution); the score difference is that of the Fitch scores of the built and the true tree.

| Set | $\mu$ | Without rebuild: score difference / recall / precision | Rebuild $r = 3$: score difference / recall / precision / time (s) |
| --- | --- | --- | --- |
| `data` | 0.0003 | 0.00% / 1.000 / 1.000 | 0.00% / 1.000 / 1.000 / 0.3 |
| `data` | 0.003 | +1.57% / 0.991 / 0.991 | −0.01% / 1.000 / 0.997 / 1.2 |
| `data` | 0.01 | +0.18% / 0.980 / 0.979 | 0.00% / 0.989 / 0.983 / 2.7 |
| `data_messy` | 0.003 | +0.16% / 0.985 / 0.991 | −0.01% / 0.999 / 0.994 / 1.3 |
| `data_messy` | 0.01 | +0.60% / 0.950 / 0.968 | +0.02% / 0.988 / 0.984 / 2.6 |
| `data_n30` | 0.01 | +0.03% / 0.993 / 0.992 | −0.01% / 0.995 / 0.993 / 4.9 |
| `data_n30_messy` | 0.01 | +0.55% / 0.962 / 0.980 | −0.01% / 0.994 / 0.992 / 5.8 |

With rebuilds of radius 3, the built trees reach the score of the true tree within 0.05% in all 16 sets (the others are not shown). The outlier without rebuild (`data`, $\mu$ = 0.003, replicate 0) had only 4 and 6 differing splits but 248 extra substitutions, a misplaced large clade. Radius 5 gives the same results at 2 to 4 times the time.

### Real data

| Dataset | Reference tree | Without rebuild | Rebuild $r = 3$ (time, s) | Rebuild $r = 5$ (time, s) |
| --- | --- | --- | --- | --- |
| `flu/h3n2/500` | 1507 | 1515 (+0.5%) | 1507 (0.0%, 0.1) | |
| `rsv/a/2000` | 42,268 | 42,470 (+0.5%) | 42,279 (+0.03%, 10) | |
| `dengue/2000` | 55,018 | 55,673 (+1.2%) | 55,114 (+0.17%, 22) | |
| `sc2/4500` | 26,175 | 26,342 (+0.6%) | 26,268 (+0.36%, 26) | 26,247 (+0.28%, 39) |
| `mpox/clade-ii/2000` | 5562 | 5599 (+0.7%) | 5592 (+0.5%, 32) | 5578 (+0.29%, 37) |
| `lassa/L/500` | 75,263 | 78,478 (+4.3%) | 77,525 (+3.0%, 18) | 77,231 (+2.6%, 94) |

On `sc2/4500` and `mpox/clade-ii/2000` about a third of the balls of radius 5 exceed the size limit because of large polytomies. On `lassa/L/500`, five insertion orders end between 76,853 and 77,810, all 2 to 3.4% above the reference tree: every start reaches a different local optimum of similar quality.

The simulations reach the true tree even at 15.6% divergence, while `lassa/L/500` does not reach the reference tree. Neither the v0 nor the v1 simulator varies the substitution rate among sites. Rate variation concentrates homoplasy at fast sites and makes the parsimony landscape rough, so simulated data with gamma-distributed site rates is needed to test the builder on Lassa-like data.

## How local the topological differences are

`fn missed_splits()` takes every informative split of a reference tree that the built tree lacks and classifies it as unresolved (compatible with all informative splits of the built tree, so a polytomy and not an error) or conflicting. For a conflicting split it finds the smallest number of tips on the wrong side over all edges of the built tree, which does not depend on the root, and the distance of each of these tips from the split in the reference tree.

**Simulated data, reference = true tree.** Without rebuilds, a tree of 400 or 1200 tips has 0 to 10 conflicting splits; most need 1 to 5 tips moved, and the moved tips lie 1 to 10 branches from the split. After rebuilds of radius 3 no conflicting split is left at $\mu \le 0.003$. At $\mu = 0.01$ (15.6% divergence), 2 to 7 remain, 1 to 20 tips each, about 200 substitutions from the split; the built trees then score at or below the true tree, so these are cases where parsimony prefers another tree, not failures of the search.

**Real data, reference = maximum-likelihood tree, built with rebuilds of radius 3.** The reference tree is not the truth, so a conflict can be an error of either tree or a resolution that the data do not determine.

| Dataset | Missed splits: unresolved / conflicting | Tips to move per conflicting split: 1 / 2 / 3-5 / 6-20 / more | Misplaced tips within 5 branches of the split | Distance of misplaced tips in substitutions: median / 90th percentile |
| --- | --- | --- | --- | --- |
| `flu/h3n2/500` | 10 / 20 | 7 / 4 / 6 / 3 / 0 | 59% | 5 / 11 |
| `mpox/clade-ii/2000` | 115 / 74 | 33 / 20 / 10 / 11 / 0 | 71% | 3 / 34 |
| `rsv/a/2000` | 80 / 92 | 53 / 16 / 14 / 8 / 1 | 76% | 29 / 74 |
| `dengue/2000` | 49 / 103 | 59 / 15 / 18 / 10 / 1 | 70% | 81 / 1023 |
| `sc2/4500` | 194 / 284 | 122 / 44 / 53 / 38 / 27 | 33% | 11 / 33 |
| `lassa/L/500` | 2 / 122 | 48 / 17 / 25 / 26 / 6 | 55% | 366 / 665 |

- About half of the conflicting splits concern one tip and three quarters at most 5 tips; 0 to 10% move more than 20 tips
- In substitutions the distances scale with the divergence of the data. In branches most misplaced tips lie within 5 branches, except on `sc2/4500`, where 41% lie more than 10 branches away but only 11 substitutions (median) away: long chains of branches with one substitution each make a radius in branches a poor measure of a neighbourhood
- The rebuilds remove 38% of the conflicting splits on `rsv/a/2000` (149 to 92) and 52% on `dengue/2000` (216 to 103), but 18% on `lassa/L/500` (149 to 122) and 3 to 10% on `sc2/4500`, `mpox/clade-ii/2000`, and `flu/h3n2/500`, where the remaining score gap to the reference tree is 0 to 0.4%

## Batch re-placement

Single tips are already pruned and re-placed by the exhaustive SPR. Batch re-placement (`reprune=<batches>`) splits the samples into random batches; for each batch it removes the tips, refreshes the rest of the tree with Fitch, places the tips again, refreshes, and keeps the result if the score drops. Unlike SPR, the backbone's internal states are re-optimized without the removed tips before they are placed.

| Dataset (Fitch score vs reference tree) | Base | 10 batches × 3 rounds | Rebuild $r = 3$ | Rebuild + 10 × 3 batches |
| --- | --- | --- | --- | --- |
| `flu/h3n2/500` | +0.53% | +0.20% | 0.00% | 0.00% |
| `rsv/a/2000` | +0.48% | +0.43% | +0.03% | +0.03% (no batch accepted) |
| `dengue/2000` | +1.19% | +0.94% | +0.17% | +0.17% (no batch accepted) |
| `sc2/4500` | +0.64% | +0.64% | +0.36% | +0.36% (no batch accepted) |
| `lassa/L/500` | +4.27% | +2.43% | +3.00% | +2.02% |

On `lassa/L/500`, 10 rounds of 10 batches reach +2.06% and smaller batches do no better. Each batch needs two Fitch refreshes, which dominates the time (96 s on `sc2/4500`). On the simulated data without rate variation, batch re-placement does not repair the misplaced clade that rebuilds of radius 3 repair.

## Rate variation among sites

The v1 simulator now draws gamma-distributed site rates (`fn gamma_site_rates()`, `fn evolve_sequences_with_rates()`), and `examples/simulate_dataset.rs` writes datasets in the layout of `validation/`: 40 sampling times 0.1 years apart with 10 tips each, constant $T_c = 10$ years, random GTR, 10,000 sites, messy observation, 2 replicates per setting. Fitch score relative to the true tree, conflicting splits of the true tree, and recall of its informative splits:

| $\mu$ | Shape $\alpha$ | Base | Rebuild $r = 3$ | Batch re-placement | Both |
| --- | --- | --- | --- | --- | --- |
| 0.003 | equal rates | +0.32%, 1.0, 0.975 | 0.00%, 0.5, 0.994 | +0.21%, 1.0, 0.982 | 0.00%, 0.5, 0.994 |
| 0.003 | 1.0 | +0.32%, 1.5, 0.964 | −0.02%, 0.5, 0.997 | +0.24%, 1.5, 0.979 | −0.02%, 0.5, 0.997 |
| 0.003 | 0.3 | +1.34%, 15.5, 0.930 | −0.01%, 1.5, 0.992 | +0.92%, 13.5, 0.948 | −0.01%, 1.5, 0.992 |
| 0.01 | equal rates | +1.08%, 14.0, 0.927 | 0.00%, 1.0, 0.993 | +0.55%, 10.0, 0.962 | 0.00%, 1.0, 0.993 |
| 0.01 | 1.0 | +1.32%, 20.5, 0.916 | +0.07%, 4.5, 0.987 | +1.04%, 15.5, 0.942 | +0.07%, 4.5, 0.987 |
| 0.01 | 0.3 | +2.54%, 45.5, 0.861 | +0.23%, 12.0, 0.961 | +2.05%, 36.5, 0.889 | +0.22%, 11.5, 0.964 |

- Strong rate variation reproduces the behavior of `lassa/L/500`: without rebuilds the score ends 2.5% above the true tree with 45 conflicting splits
- At $\mu = 0.01$, $\alpha = 0.3$, rebuilds of radius 4 end +0.19% and +0.07% above the true tree; radius 5 and 6 end higher (+0.36%, +0.34%), because the rebuild of a large ball is itself a heuristic
- Refined from the true tree, the builder's moves lower its score by 5 substitutions while changing 2 to 3 splits: the parsimony optimum lies next to the true tree, and the remaining gap of the built trees is a search failure
- The remaining 8 to 9 conflicting splits move 1 or 2 tips each, about 180 substitutions away, and one clade of more than 20 tips per replicate. Moving a subtree, or a pseudo-leaf of a rebuild, keeps the state of its root fixed; moves that need a different root state are not found

## State sets, insertion order, and starting trees

**State sets.** Queries now hold a set of states per site. Ambiguity codes are scored as their sets, as Fitch does, instead of as missing data. SPR moves an internal subtree with the plurality set of its children at every site where they disagree, and rebuilds treat internal pseudo-leaves the same way. After a move, the subtree root takes a state from its set and the substitutions on its child edges are rewritten (`fn Mat.restate()`).

**Insertion order.** Comparing the two versions exposed a large effect of the insertion order. On the simulation with $\mu = 0.01$ and $\alpha = 0.3$ (true tree 23,321):

| Order | Score after placement | Final score (SPR and Fitch refresh) |
| --- | --- | --- |
| Fewest missing sites first | 28,142 | 24,427 |
| Input order (sampling time, most recent first) | 28,144 | 24,360 |
| Random, 3 seeds | 25,957 to 26,172 | 23,740 to 23,934 |
| Farthest-first scaffold of 20 to 200 anchors, then random | 27,261 to 28,027 | 23,936 to 24,235 |

Orders that group similar samples, by sampling time or by completeness, build a poor early backbone; farthest-first anchors start with outliers and do no better. The example program now uses random order by default.

**Starting trees.** Two alternatives to placing all samples one by one: a scaffold (`scaffold=<size>`), in which a random subset is built from several orders with full refinement and the rest is placed onto the best result; and neighbor joining (`initial=nj`) on Jukes-Cantor distances over the sites observed in both samples, converted to a MAT by the Fitch refresh.

On the rate-variation simulations, placement in random order with rebuilds of radius 3 and neighbor joining with rebuilds of radius 3 both end within 0.02% of the true tree's score in all settings ($\mu = 0.01$, $\alpha = 0.3$: +0.02% and 7 conflicting splits, and −0.02% and 8.5).

Real data, Fitch score relative to the reference tree, all with state sets and rebuilds of radius 3 (conflicting splits in parentheses):

| Dataset | Placement, random order | Neighbor joining | Neighbor joining and batch re-placement |
| --- | --- | --- | --- |
| `flu/h3n2/500` | −0.07% (9) | −0.07% (9) | |
| `rsv/a/2000` | −0.01% (80) | +0.14% (130) | +0.14% (130) |
| `dengue/2000` | +0.05% (92) | +0.12% (88) | +0.08% (91) |
| `mpox/clade-ii/2000` | +0.16% (82) | +0.31% (98) | |
| `sc2/4500` | +0.29% (291) | +0.16%, top hits (see below) | |
| `lassa/L/500` | +3.09% (128) | +0.96% (93) | +0.58% (84) |

- On `lassa/L/500`, scaffolds of 30 to 120 samples brought placement only to +2.2%. Neighbor joining sees all pairwise distances at once and gets the deep structure that sequential placement misses. The excess of the placement trees is spread evenly along the genome (100 to 200 extra substitutions per 500 sites), which argues against recombination as its cause
- On `rsv/a/2000`, `dengue/2000`, and `mpox/clade-ii/2000`, placement in random order is as good as or better than neighbor joining; the built trees of `flu/h3n2/500` and `rsv/a/2000` score slightly below the reference trees
- State sets made the rebuilds 2 to 4 times slower (the plurality sets of every pseudo-leaf), for example 85 s on `sc2/4500`

## Fast and local neighbor joining

**Top hits.** `fn neighbor_joining()` takes an optional number $m$ of top hits, as in FastTree:

- Each node keeps the $m$ partners that are best by the joining criterion when its list is made, and caches its best partner among them.
- A joined node inherits the lists of its two children and enters the lists of its own partners.
- All lists are rebuilt in parallel every $m$ joins.

The pair search then costs $O(n^3/m)$ instead of $O(n^3)$. The distance matrix still takes $O(n^2)$ time and memory. With $m \approx \sqrt n$ the trees are identical to exact neighbor joining on `lassa/L/500` ($m = 20$) and `dengue/2000` ($m = 44$): the same Fitch score before and after refinement. Lists ranked by raw distance instead of the joining criterion, or without the reciprocal entries and the rebuilds, gave starting trees 1.6 to 2.4% worse on `lassa/L/500`.

Neighbor-joining start, without refinement (`spr=0`):

| Dataset | $m$ | Time (s) | Exact neighbor joining |
| --- | --- | --- | --- |
| `lassa/L/500` | 20 | 1.2 | same tree |
| `dengue/2000` | 44 | 7.1 | same tree |
| `sc2/4500` | 72 | 16.6 | not run |

**Local neighbor joining.** `rebuild_nj` adds a neighbor-joining tree of the pseudo-leaves to the candidates of every local rebuild. The tree is fitted by Fitch on the variable columns of the ball and attached below the ball's top node. Fitch score relative to the reference tree, rebuilds of radius 3 unless stated:

| Dataset | Placement | Placement, local NJ | NJ start | NJ start, local NJ |
| --- | --- | --- | --- | --- |
| `flu/h3n2/500` | −0.07% | −0.07% | −0.07% | −0.07% |
| `rsv/a/2000` | −0.01% | −0.01% | +0.14% | +0.13% |
| `dengue/2000` | +0.05% | +0.05% | +0.12% | +0.09% |
| `mpox/clade-ii/2000` | +0.16% | +0.16% | +0.31% | +0.23% |
| `sc2/4500` | +0.29% | +0.31% | +0.16% | +0.13% |
| `lassa/L/500` | +3.09% | +2.19% | +0.96% | |
| `lassa/L/500`, radius 5, up to 64 pseudo-leaves | +3.20% | +1.28% (seeds 0 to 2: +1.28, +2.49, +1.53%) | | +0.55% |
| `lassa/L/500`, radius 8, up to 128 pseudo-leaves | | +2.77% | | |

Observations:

- On low-divergence data, local NJ changes nothing. Exhaustive placement in 4 orders already finds the best tree of a ball with at most 24 pseudo-leaves
- On `lassa/L/500`, local NJ removes a third to a half of the excess, but even the best placement build stays above the global NJ start. Balls of radius 8 improve less, because a rebuild is kept only if it beats the current ball, and NJ followed by SPR on 128 pseudo-leaves rarely does. The errors of the placement trees on Lassa are not local at the scale that a ball can hold
- The global NJ start with local NJ rebuilds of radius 5 is the best build of `lassa/L/500` so far: +0.55%, without batch re-placement
- With top hits, the NJ start is also the better start on `sc2/4500` (+0.16% against +0.29%). Which start wins differs between datasets, by at most 0.15% except on Lassa

**Fitch refresh on variable columns.** The refresh now encodes the leaves from their sparse queries and runs Fitch only on the columns where some sample differs from the reference. The other columns keep the reference state. Scores are identical to the refresh on full sequences, and the build of `mpox/clade-ii/2000` takes 13 s instead of 24 s. Its peak memory (0.9 GB) comes from reading the alignment, not from the refresh. The same function fits the local neighbor-joining trees.

## Site partitions and reweighting on divergent data

On `lassa/L/500` many sites are saturated. On the neighbor-joining start, 1994 of the 10,402 sites change 20 times or more. They carry 70% of the substitutions. The idea tested here: partition the sites by their number of changes on a rough tree, build the global clustering from the slow sites, then refine with all sites. The experiments use the `build_tree` example on modified alignments. Masked sites are written as `N`, and integer weights are emulated by duplicating columns. `start=<nwk>` then refines the result on the full alignment. The rough tree is neighbor joining without refinement. Both stages use `initial=nj spr=10 spr_slack=all refresh rebuild=3`. The table counts reference splits missing from the built tree (95 for the unmodified build).

| Sites | Weighted stage | After refining on all sites |
| --- | --- | --- |
| Only sites with at most 2, 5, 10, 20 changes | 78 to 194 of 128 to 385 splits; precision 42 to 52% | 142 to 165; Fitch 79,846 to 82,141 |
| Weight 4, 3, 2, 1 for at most 3, 8, 16, more changes | 122 | 117; Fitch 76,631 |
| Weight $\mathrm{round}(16/c)$ | 127 | 123; Fitch 77,799 |
| Drop sites with more than 25, 30, 40, 60 changes | 127, 98, 91, 92 | 109, 90, 88, 88 |
| Half weight for sites with more than 40 changes | 87 | 87; Fitch 75,839 |

- **Slow sites alone do not give the global clustering.** Few changes mostly means a rare allele: 583 of the 847 sites with at most 2 changes change once and are uninformative. The splits from slow sites are less precise (42 to 52%) than the rough tree on all sites (58%). The masked trees are poorly resolved, and their polytomies leave the full-alignment refinement at 80 to 82k against 76.0k
- **Up-weighting slow sites hurts.** Every scheme that gives slow sites more weight ends 20 to 30 splits further from the reference. The fast sites carry most of the topological signal
- **Down-weighting the fastest sites helps marginally.** Over three input orders, half weight for sites with more than 40 changes ends with 89, 92, 90 missing splits, against 93, 92, 90 without. The Fitch score is 75,898 against 75,969 on average. Dropping only 13 sites (more than 60 changes) moves the result as much as dropping 131, so most of the difference is a perturbation of the search
- **The gap is in the search, not the objective.** Refined by the builder's moves, the reference tree scores 75,258. That is 0.8% below every built tree, and it keeps 406 of its 411 splits. Of the 92 reference splits the best built tree misses, 66 sit on branches longer than 10 substitutions. The parsimony optimum lies next to the ML tree, and the built trees are stuck in other local optima
- **The parsimony ratchet makes slow progress.** Each iteration doubles the weight of a random 25% of sites, refines from the current tree, then refines on equal weights, and keeps the tree if its score does not rise. Eight iterations (6 min) lower the score from 75,877 to 75,790 and the missing splits from 89 to 82. The tree stays 0.7% above the reference

## Dating with built trees

The built Newick trees run through `treetime clock` and `treetime timetree` unchanged. On `sc2/4500` the clock estimates differ from those of the reference tree:

| Tree | Rooting | Clock rate | Root date | Correlation r |
| --- | --- | --- | --- | --- |
| Reference | least squares | 9.1e-4 | 2019.81 | 0.956 |
| Built, fewest missing first | least squares | 7.0e-4 | 2018.17 | 0.953 |
| Built, oldest first | least squares | 7.1e-4 | 2018.99 | 0.957 |
| Reference | outgroup Wuhan-Hu-1/2019 | 1.04e-3 | 2020.02 | 0.963 |
| Built, fewest missing first | outgroup Wuhan-Hu-1/2019 | 1.16e-3 | 2019.94 | 0.924 |
| Built, oldest first | outgroup Wuhan-Hu-1/2019 | 9.4e-4 | 2019.75 | 0.965 |

(`treetime clock` without the clock filter.) With the same outgroup root, the median root-to-tip distance of a tip differs from that in the reference tree by +24 substitutions when the samples with the fewest missing sites are placed first, and by −12 substitutions when the oldest samples are placed first. The parsimony scores of the two built trees differ by 9 substitutions. In the first order, most placed samples are from 2025 and 2026, and the early lineages join a backbone built around them; pairs of substitution and reversion on that backbone lengthen every root-to-tip path through it at almost no parsimony cost. Least-squares rerooting then picks a different root.

With the later builder (state sets, random insertion order, rebuilds of radius 3; Fitch score 26,250), the outgroup-rooted estimates agree with the reference tree: clock rate 9.95e-4, root date 2019.80, r = 0.962, and a median root-to-tip difference of −2 substitutions, although 1364 of 5100 tips differ by more than 5. Least-squares rerooting still gives 7.4e-4 and 2019.01. Collapsing the branches without mutations of the reference tree (`treetime prune --collapse-empty`) leaves its least-squares estimate unchanged (9.1e-4, 2019.81), so the polytomies of the built tree do not explain the difference; the least-squares root responds to the tips whose root-to-tip distances differ.

Consequences:

- Insertion order changes the dating more than it changes the parsimony score. Date order, which [feat-tree-infer.md](feat-tree-infer.md) proposes, is the better default, but it did not reproduce the reference clock estimates either
- The reference tree is not ground truth for dating. Simulations with known dates and rates ([Validation approach](feat-tree-infer.md#validation-approach)) are needed to decide which tree gives the better estimate
- `timetree` on the date-ordered built tree stopped with the internal error of [kb/issues/M-timetree-polytomy-resolution-fails-on-backward-branch.md](../issues/M-timetree-polytomy-resolution-fails-on-backward-branch.md). Built trees have many large polytomies and reach this error deterministically

## Changes to the plan

- **No separate initial-topology method.** Placement starts from the root and one sample, so the statement "Placement cannot build the first tree" of [feat-tree-infer.md](feat-tree-infer.md) does not hold for a MAT rooted at a reference sequence
- **Exact placement first, heuristics later.** Branch-and-bound makes exact placement affordable at 10^3 to 10^4 samples. Centroid starts belong to a later stage for 10^5 samples and more, as a seed of the bound and with a fallback, not as a replacement of the exact search
- **Fitch refresh as a pipeline step.** The existing Fitch pass improves the builder's internal states and should run between rearrangement rounds; its full-sequence memory cost makes the mutation-first representation a prerequisite for large genomes
- **Local rebuilds.** Rebuilding balls of radius 3 reaches the score of the true tree on simulated data and closes most of the gap to the reference trees on real data
- **Random insertion order, and neighbor joining for divergent data.** Sorted insertion orders build poor backbones. For data as divergent as Lassa, a neighbor-joining start followed by the same refinement beats every placement-based start, and local neighbor joining in the rebuilds does not replace it. For low-divergence data the two starts differ by at most 0.15%. Neighbor joining with top hits removes the cubic search, but its $O(n^2)$ distance matrix limits it to about $10^4$ samples
- **Dating in the acceptance criteria.** Parsimony and split agreement do not detect the backbone differences that change the clock estimates
- **Graph operations.** `fn Graph.remove_edge()` and `fn Graph.remove_node()` scanned all nodes or edges and now take time proportional to the degree. `fn Graph.add_edge()` still checks all children of the source for a duplicate edge, which is linear in the size of a polytomy

## Limitations of the prototype

- Nucleotides only; indels are not placed
- The SPR runs on one thread and moves subtrees with a fixed root state
- Insertion order by date is available only through a name list in the example program
- No branch support values, no command line interface, and no divergence test
- LBI parameters (time scale, pseudo-length) are not tuned

## Reproduction

```bash
./dev/docker/run cargo build --release -p treetime --example build_tree
target/release/examples/build_tree aln=data/sc2/4500/aln.fasta.xz tree=data/sc2/4500/tree.nwk \
  search=exhaustive spr=10 spr_slack=all refresh out=tmp/sc2_4500_built.nwk
```

`compare` adds the comparison of every placement with the full scan; `search=walk refs=lbi max_refs=32 slack=6` selects the greedy walk, and `exact` seeds the branch-and-bound with it; `search=references refs=lbi` or `search=references refs=branches cuts=64` starts the exact search at the nearest reference; `order_file=<file>` places the listed names first, in the listed order; `rebuild=3` adds local rebuilds of radius 3; `start_from_reference` refines the reference tree instead of building one; `sites` compares the built and the reference tree site by site; `locality` reports the missed splits in both directions; `reprune=10 reprune_rounds=3` adds batch re-placement; `order=random|scaffold`, `scaffold=60`, and `initial=nj` choose the insertion order and the starting tree; `sites_out=<tsv>` writes per-site Fitch counts; `start=<nwk>` refines a given tree on the alignment. `examples/simulate_dataset.rs out=<dir> mu=0.01 alpha=0.3 messy` writes a simulated dataset with its true tree.
