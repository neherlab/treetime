# Polytomy resolution removes single-child nodes that carry a date

## v0 location

`TreeTime.resolve_polytomies()` (`#TreeTime`, `#resolve_polytomies`) [packages/legacy/treetime/treetime/treetime.py#L701-L705](../../packages/legacy/treetime/treetime/treetime.py#L701-L705)

## Erratum

After it resolves the polytomies of a round, v0 removes every non-root node with one child:

```python
obsolete_nodes = [n for n in self.tree.find_clades() if len(n.clades) == 1 and n.up is not None]
for node in obsolete_nodes:
    self.logger('TreeTime.resolve_polytomies: remove obsolete node ' + node.name, 4)
    if node.up is not None:
        self.tree.collapse(node)
```

The selection uses the number of children alone. A single-child node with a date constraint (a sampled ancestor, or a dated point on a lineage) is removed together with its date. The removal runs whenever polytomy resolution runs, also in rounds that resolve no polytomy, and polytomy resolution is on by default.

An undated single-child node carries no data, and its removal is exact: substitution probabilities satisfy $P(t_1)P(t_2) = P(t_1 + t_2)$, and root-to-tip distances and their variances add along a path. A dated single-child node constrains the time of its point on the lineage. Its date cannot move to the child, because the child is a later point, separated by a branch of unknown duration, which is the quantity the time inference estimates.

On the tree from [neherlab/treetime#959](https://github.com/neherlab/treetime/issues/959) with an inserted single-child node `u1` dated 2020-03-01, v0 omits `u1` from `dates.tsv` and dates its child as if the date were absent. With `--keep-polytomies`, `u1` keeps its date and its child moves from 2020-05-20 to 2020-04-26.

## Evidence

- The adjacent setup code reads a date for every named node, internal nodes included, and keeps a dated internal node out of the bad branches ([packages/legacy/treetime/treetime/clock_tree.py#L121-L152](../../packages/legacy/treetime/treetime/clock_tree.py#L121-L152)). The time inference then turns the date into a constraint of that node ([L378-L387](../../packages/legacy/treetime/treetime/clock_tree.py#L378-L387)). v0 thus treats the date of an internal node as data
- Where v0 ignores a date on purpose, it warns: "Date constraint will be ignored!" ([packages/legacy/treetime/treetime/clock_tree.py#L389-L396](../../packages/legacy/treetime/treetime/clock_tree.py#L389-L396)). The polytomy cleanup drops the date without a warning
- The cleanup calls the removed nodes obsolete. That holds only for a node without data, which a dated node is not

## v0 impact

- Only trees with a dated single-child node are affected, which requires a date for an internal node in the metadata
- The removed node's date does not constrain the times of its child or its ancestors, and the node is missing from the outputs

## v1 status

v1 removes only undated single-child nodes after polytomy resolution, in `fn remove_undated_single_child_nodes()` of [packages/treetime/src/timetree/optimization/polytomy/resolve.rs](../../packages/treetime/src/timetree/optimization/polytomy/resolve.rs). A dated single-child node stays in the tree with its date constraint.
