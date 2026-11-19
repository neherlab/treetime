# TreeKnit output to named segment trees

TreeTime merges segment trees by node name ([kb/proposals/edge-scoped-inheritance.md](../proposals/edge-scoped-inheritance.md)). TreeKnit writes one resolved tree per segment with internal names local to each tree, and the maximum compatible clades (MCCs) of every pair of segments. This report records how far the TreeKnit output determines shared node identities, measured on `data/flu/treeknit/` (H3N2, HA, NA, PB1, PB2, 1564 leaves, all segments present for every leaf) with the prototype [test_scripts/treeknit_named_trees.py](../../test_scripts/treeknit_named_trees.py).

## Method

1. For each segment pair, assign each internal node of both trees an MCC with the rule of v0 (`assign_mccs()` in `packages/legacy/treetime/treetime/arg.py`): Fitch on the leaf MCCs, then the parent's MCC where the children allow it
2. For each MCC node of one tree, take the MCC leaves below it and match it with the MRCA of the same leaves in the other tree, when that MRCA has the same MCC and the match is mutual
3. Join the matches of all six pairs into node identities (union-find)
4. Split identities that hold two nodes of one tree, and identities that close a cycle in the merged graph (segments that order the same two identities in opposite ways), back into separate per-tree nodes

## Results

| Quantity                                | Value                                                          |
| --------------------------------------- | -------------------------------------------------------------- |
| Internal nodes per tree                 | 747 to 783 (resolved trees keep some polytomies)               |
| Matched nodes per pair                  | 558 to 596; 28 to 44 MCC nodes per pair have no mutual match   |
| Identities shared by 4 / 3 / 2 segments | 508 / 92 / 99                                                  |
| Internal nodes in one segment only      | 518                                                            |
| Identities with two nodes of one tree   | 10 (51 nodes)                                                  |
| Identities split to remove cycles       | 3                                                              |
| Merged graph                            | 2839 nodes, 4374 edges, 4 roots; 802 nodes with 2 to 4 parents |

`treetime ancestral --segment-trees --segment-alignments --method-anc parsimony` runs on the four named trees in 0.3 s and 60 MB, and gives each segment the same reconstructed sequences as a run on its own tree.

With pass-through insertion ([kb/algo/graph.md](../algo/graph.md#merging-trees-by-name)) the merge inserts 179 nodes into segments that pass through them, and leaves 9 edges whose other segments disagree on the path and 29 edges whose path holds a node of the segment elsewhere. The graph keeps 2839 nodes and has 4268 edges instead of 4374; 737 nodes keep 2 to 4 parents (289 with 2, 200 with 3, 248 with 4). The datasets of [data/flu/h3n2/segments/](../../data/flu/h3n2/segments/) hold these named trees, with subsets of 100 and 500 viruses.

## What the TreeKnit output does not determine

- **Identities across more than two segments**: MCCs are pairwise, and matches joined across pairs are not always consistent. Here 13 identities (10 that held two nodes of one tree, 3 that closed a cycle) had to be split by a heuristic. One set of node identities over all segments, from TreeKnit, would remove the heuristic
- **Where a reassortant lineage attaches, and which nodes are not reassortant**: an MCC root has a different parent in each segment tree, but TreeKnit does not say where on the branches of the other segments the lineage joins. Pass-through insertion recovers the nodes on a path that one segment shows and another skips; it cannot merge two lineages that both have intermediate nodes, because topology alone does not tell whether their child is a reassortant. An MCC member below its MCC root is not a reassortant, which TreeKnit knows
- **Backbone and roots**: nodes outside every MCC, including the four segment roots, stay per segment. Whether they mark reassortment or uncertain topology is not in the output
- **Leaves in singleton MCCs**: about half of the MCCs of each pair are single leaves, so many leaves have a different parent in each segment, which the graph represents as reassortment

## Possible upstream additions to TreeKnit

- Global node identities: one name per ancestor shared by several segment trees, consistent over all segments
- Attachment points of reassortment: for each MCC root and each other segment, the edge (and position) where its lineage joins, which gives the nodes that a segment passes through
