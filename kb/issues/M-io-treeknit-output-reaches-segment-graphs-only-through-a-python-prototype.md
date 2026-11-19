# TreeKnit output reaches segment graphs only through a Python prototype

`--segment-trees` merges segment trees by the names of their internal nodes ([kb/algo/graph.md](../algo/graph.md#merging-trees-by-name)). TreeKnit writes one tree per segment with names that are local to each tree, and the maximum compatible clades (MCCs) of every pair of segments. The only step from this output to named trees is the prototype [test_scripts/treeknit_named_trees.py](../../test_scripts/treeknit_named_trees.py): it matches nodes pairwise by MCC leaves, joins the matches of all pairs, and splits 13 identities on the H3N2 data with a heuristic (10 that held two nodes of one tree, 3 that closed a cycle). Method and results: [kb/reports/2026-11-15-treeknit-output-to-named-segment-trees.md](../reports/2026-11-15-treeknit-output-to-named-segment-trees.md).

So users of TreeKnit can run `--segment-trees` only with a test script whose identities depend on a heuristic. TreeTime also does not read an ARG as one file: the annotated network Newick with segment sets on the edges (`[&segments={0,1}]`, as CoalRe and treeknit-rs write it) is parsed by the Newick dialect tests, but no command takes it as segment input.

The heuristic belongs in TreeKnit: treeknit-rs has the same open problem for more than two trees, and its proposed fix (a union-find over the nodes of all trees, where a class with two nodes of one tree or a cycle is a conflict) is the algorithm of the prototype. The report lists the information that TreeKnit has and does not write: one set of node identities over all segments, and the attachment points of reassortant lineages, which would also give the pass-through nodes that the merge now infers from topology ([N-graph-merge-reroutes-segment-edges-past-reassortments-in-complete-trees.md](N-graph-merge-reroutes-segment-edges-past-reassortments-in-complete-trees.md)).

> [!IMPORTANT]
> **Decision required.** The exchange format between TreeKnit and TreeTime:
>
> - Named segment trees, one per segment, with node names that TreeKnit makes consistent over all segments (the current input of `--segment-trees`)
> - One ARG file, such as annotated network Newick with segment sets on the edges, that TreeTime reads into a graph with edge scopes
>
> Both need global node identities from TreeKnit. The choice is shared with treeknit-rs.

## Validation

On `data/flu/treeknit/`: the TreeKnit output goes into `timetree --segment-trees` without the prototype, and `ancestral --segment-trees --method-anc parsimony` gives each segment the reconstruction of a run on its own tree, as the report shows for the prototype's trees.
