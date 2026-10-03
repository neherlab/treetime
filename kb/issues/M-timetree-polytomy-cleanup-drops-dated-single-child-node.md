# Timetree polytomy cleanup removes single-child nodes that carry a date

After polytomy resolution, `remove_single_child_nodes()` removes every node with one parent and one child, selected by degree alone ([packages/treetime/src/timetree/optimization/polytomy/resolve.rs#L212-L236](../../packages/treetime/src/timetree/optimization/polytomy/resolve.rs#L212-L236)). A single-child node with a date constraint carries data: the date constrains the time of that point on the lineage. Removing the node drops the constraint without a message.

## Scientific context

- An undated single-child node carries no data. Its removal is exact: substitution probabilities satisfy $P(t_1)P(t_2) = P(t_1 + t_2)$, and root-to-tip distances and their variances add along a path
- A dated single-child node is a sampled ancestor or a dated point on a branch. Its date cannot move to the child or to a leaf: the child is a later point, separated by a branch of unknown duration, which is the quantity the time inference estimates. A branch without mutations still has a positive duration, so such a transfer is never exact
- Converting the node into a zero-length leaf keeps the date but changes the model: the branch to the new leaf gets a positive duration, and the coalescent prior counts an extra lineage

## v0 comparison

v0 has the same defect: `TreeTime.resolve_polytomies()` removes every non-root single-child node, dated or not ([packages/legacy/treetime/treetime/treetime.py#L701-L705](../../packages/legacy/treetime/treetime/treetime.py#L701-L705)). On the tree from [neherlab/treetime#959](https://github.com/neherlab/treetime/issues/959) with an inserted single-child node `u1` dated 2020-03-01, v0 omits `u1` from `dates.tsv` and dates its child as if the date were absent. With `--keep-polytomies`, `u1` keeps its date and its child moves from 2020-05-20 to 2020-04-26.

## v1 status

v1 has the same defect. `data/smoke/flu-h3n2-20-tree-single-child-polytomy.nwk` is `data/flu/h3n2/20/tree.nwk` with a single-child node `U` and a polytomy. With `U` dated 2008.5 in the metadata, `--resolve-polytomies` removes `U` from the tree and from the node data; `--keep-polytomies` keeps `U` at 2008.5:

```bash
cp data/flu/h3n2/20/metadata.tsv tmp/metadata.tsv && printf 'U\t2008.5\n' >> tmp/metadata.tsv
treetime timetree --resolve-polytomies --tree=data/smoke/flu-h3n2-20-tree-single-child-polytomy.nwk \
  --dates=tmp/metadata.tsv --aln=data/flu/h3n2/20/aln.fasta.xz --output-all=<dir>
```

## Proposed fix

Remove only single-child nodes without a date constraint. Keep dated ones in the graph, so their constraints stay in the time inference.
