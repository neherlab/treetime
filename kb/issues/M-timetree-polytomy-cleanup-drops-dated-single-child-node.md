# Timetree polytomy cleanup removes single-child nodes that carry a date

After polytomy resolution, `remove_single_child_nodes()` removes every node with one parent and one child, selected by degree alone ([packages/treetime/src/timetree/optimization/polytomy/resolve.rs#L212-L236](../../packages/treetime/src/timetree/optimization/polytomy/resolve.rs#L212-L236)). A single-child node with a date constraint carries data: the date constrains the time of that point on the lineage. Removing the node drops the constraint without a message.

## Scientific context

- An undated single-child node carries no data. Its removal is exact: substitution probabilities satisfy $P(t_1)P(t_2) = P(t_1 + t_2)$, and root-to-tip distances and their variances add along a path
- A dated single-child node is a sampled ancestor or a dated point on a branch. Its date cannot move to the child or to a leaf: the child is a later point, separated by a branch of unknown duration, which is the quantity the time inference estimates. A branch without mutations still has a positive duration, so such a transfer is never exact
- Converting the node into a zero-length leaf keeps the date but changes the model: the branch to the new leaf gets a positive duration, and the coalescent prior counts an extra lineage

## v0 comparison

v0 has the same defect: `TreeTime.resolve_polytomies()` removes every non-root single-child node, dated or not ([packages/legacy/treetime/treetime/treetime.py#L701-L705](../../packages/legacy/treetime/treetime/treetime.py#L701-L705)). On the tree from [neherlab/treetime#959](https://github.com/neherlab/treetime/issues/959) with an inserted single-child node `u1` dated 2020-03-01, v0 omits `u1` from `dates.tsv` and dates its child as if the date were absent. With `--keep-polytomies`, `u1` keeps its date and its child moves from 2020-05-20 to 2020-04-26.

## v1 status

Code evidence only. The v1 run that would show the dropped date fails earlier: [H-timetree-resolve-polytomies-fails-with-single-child-node.md](H-timetree-resolve-polytomies-fails-with-single-child-node.md). Without a polytomy, the cleanup does not run and v1 keeps a dated single-child node with its date (`U` at 2008.5 on `data/flu/h3n2/20` with the node inserted as in that issue).

## Proposed fix

Remove only single-child nodes without a date constraint. Keep dated ones in the graph, so their constraints stay in the time inference.
