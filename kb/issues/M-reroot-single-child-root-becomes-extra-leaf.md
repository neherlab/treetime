# Rerooting turns a single-child root into an extra leaf

When the input root has exactly one child and the reroot moves the root, the old root ends up as a leaf of the rerooted tree. The leaf has no date and no sequence, and it is written to the outputs as if it were a sample.

## Reproduction

```bash
d=data/flu/h3n2/20
mkdir -p tmp/single-child

# Add a single-child root STEM above the input root
sed -E 's/^(.*);$/(\1:0.001)STEM:0;/' $d/tree.nwk > tmp/single-child/tree_stem.nwk

./dev/docker/run just r treetime clock --tree=tmp/single-child/tree_stem.nwk --metadata=$d/metadata.tsv \
  --output-all=tmp/single-child/clock
```

`clock.nwk` has 20 leaves for the 19 samples, including `STEM`. `clock.clock.csv` has a row for `STEM` without a date. `--reroot-tips='A/Indiana/03/2012|KC892731|04/03/2012|USA|11_12|H3N2/1-1409'` gives the same result.

When the reroot keeps the old root position, nothing changes: on the tree from [neherlab/treetime#959](https://github.com/neherlab/treetime/issues/959), the least-squares reroot keeps `i373496` as a single-child root, and `--reroot-tips=255211` makes it a leaf.

## Analysis

The reroot inverts the edges on the path from the old root to the new root. A single-child old root loses its only outgoing edge and gains one incoming edge, so it becomes a leaf. `remove_node_if_trivial()` removes only a node with one incoming and one outgoing edge ([packages/treetime-graph/src/reroot.rs#L83-L93](../../packages/treetime-graph/src/reroot.rs#L83-L93)), so the leaf stays ([packages/treetime/src/clock/reroot.rs#L90-L93](../../packages/treetime/src/clock/reroot.rs#L90-L93)).

v0 has the same defect through Biopython's `Tree.root_with_outgroup()` and crashes on it in the root-to-tip plot: [kb/v0-errata/clock-bad-branch-numpy-bool-identity-check.md](../v0-errata/clock-bad-branch-numpy-bool-identity-check.md).

## Scientific context

A single-child root is valid input: it marks a stem above the most recent common ancestor of the samples, for example in a subtree cut from a larger tree. Rerooting replaces the root, so the old root has no role in the rerooted tree unless it carries data.

- An undated single-child root carries no data. Removing it before the reroot is exact: substitution probabilities satisfy $P(t_1)P(t_2) = P(t_1 + t_2)$, root-to-tip distances and their variances add along a path, and a branch to a node without data contributes a factor of 1 to the likelihood
- A dated single-child root carries a date. After the reroot it is a dated leaf, which the inference handles like any other sample
- With `--keep-root`, the stem is part of the requested tree and stays

## Proposed fix

Before rerooting, remove an undated single-child root by making its child the new root. Keep a dated one, and leave `--keep-root` runs unchanged. After the reroot, fail if a leaf without a date was created.
