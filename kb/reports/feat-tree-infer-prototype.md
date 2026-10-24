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

## Reference regions

A second use of reference nodes separates their selection from the search, and does not depend on the root of the MAT:

- **Selection.** Either the LBI peaks with the time scale $\tau = r / \ln n$, where $r$ is the radius of the tree in substitutions and $n$ the number of samples, or the node of least eccentricity in each component left after cutting the $k$ longest internal branches
- **Cells.** Every node belongs to the nearest reference by substitutions along the tree (Voronoi cells)
- **Search.** The query is encoded relative to the nearest reference $c$, and the cell of $c$ is searched without regard to the root. The best placement in the cell, with cost $b$, is accepted when $b + d_{\mathrm{obs}}(q, c) + \lvert \mathrm{missing}(q) \cap P \rvert < m_c$, where $P$ is the set of polymorphic sites and $m_c$ the smallest distance from $c$ to a placement outside the cell. By the triangle inequality no outside placement can then match $b$. Otherwise the branch-and-bound runs with $b$ as its bound

The result equals the exhaustive search in all tests and on all datasets. Placement only, batches of 64:

| Dataset | Search | Time (s) | Node visits | References per batch | Certified without full search | Mean differences per query, global / local frame |
| --- | --- | --- | --- | --- | --- | --- |
| `lassa/L/500` | branch-and-bound | 0.32 | 21,132 | | | 825 |
| | LBI regions | 0.38 | 66,886 | 2 | 46% | 825 / 768 |
| | 64 cut branches | 0.67 | 40,492 | 56 | 16% | 825 / 468 |
| `dengue/2000` | branch-and-bound | 0.46 | 103,312 | | | 519 |
| | LBI regions | 0.82 | 642,584 | 5 | 28% | 519 / 329 |
| | 64 cut branches | 1.25 | 287,245 | 63 | 8% | 519 / 149 |
| `mpox/clade-ii/2000` | branch-and-bound | 4.95 | 157,908 | | | 15 |
| | LBI regions | 5.60 | 1,493,394 | 2 | 68% | 15 / 12 |
| `sc2/4500` | branch-and-bound | 1.43 | 781,762 | | | 58 |
| | LBI regions | 2.90 | 6,357,005 | 6 | 3% | 58 / 18 |
| | 64 cut branches | 2.73 | 1,898,303 | 64 | 1% | 58 / 19 |

Observations:

- The local frame shortens the difference lists by up to 3.5 times, but the cost of placement is set by the substitutions on the visited edges, not by these lists
- The search inside a cell has no bound and visits the whole cell, and uncertified queries search again from the root. At these sizes, regions are slower than the branch-and-bound from the root
- The certificate needs the query within about half the isolation of its reference. Few large cells (LBI) certify more queries but cost more per cell; many small cells (cut branches) are cheap to search but certify few
- With $\tau = r / \ln n$ the LBI has 2 to 6 peaks per batch, fewer than the lineages of Lassa or dengue
- The test of the idea at the scale it targets, 10^5 samples and more, and a bound inside the cells remain open

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

Consequences:

- Insertion order changes the dating more than it changes the parsimony score. Date order, which [feat-tree-infer.md](feat-tree-infer.md) proposes, is the better default, but it did not reproduce the reference clock estimates either
- The reference tree is not ground truth for dating. Simulations with known dates and rates ([Validation approach](feat-tree-infer.md#validation-approach)) are needed to decide which tree gives the better estimate
- `timetree` on the date-ordered built tree stopped with the internal error of [kb/issues/M-timetree-polytomy-resolution-fails-on-backward-branch.md](../issues/M-timetree-polytomy-resolution-fails-on-backward-branch.md). Built trees have many large polytomies and reach this error deterministically

## Changes to the plan

- **No separate initial-topology method.** Placement starts from the root and one sample, so the statement "Placement cannot build the first tree" of [feat-tree-infer.md](feat-tree-infer.md) does not hold for a MAT rooted at a reference sequence
- **Exact placement first, heuristics later.** Branch-and-bound makes exact placement affordable at 10^3 to 10^4 samples. Centroid starts belong to a later stage for 10^5 samples and more, as a seed of the bound and with a fallback, not as a replacement of the exact search
- **Fitch refresh as a pipeline step.** The existing Fitch pass improves the builder's internal states and should run between rearrangement rounds; its full-sequence memory cost makes the mutation-first representation a prerequisite for large genomes
- **Stronger SPR.** Subtree moves need state sets at the subtree root, moves of sibling groups, and a treatment of ambiguity codes that agrees with Fitch
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

`compare` adds the comparison of every placement with the full scan; `search=walk refs=lbi max_refs=32 slack=6` selects the greedy walk, and `exact` seeds the branch-and-bound with it; `search=regions refs=lbi` or `search=regions refs=branches cuts=64` selects the region search; `order_file=<file>` places the listed names first, in the listed order.
