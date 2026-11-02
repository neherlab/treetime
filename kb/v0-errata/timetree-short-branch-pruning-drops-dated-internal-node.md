# Short-branch pruning removes internal nodes that carry a date

## v0 location

`TreeAnc.prune_short_branches()` (`#TreeAnc`, `#prune_short_branches`) [packages/legacy/treetime/treetime/treeanc.py#L1475-L1495](../../packages/legacy/treetime/treetime/treeanc.py#L1475-L1495), called by `TreeTime._run()` before the timetree loop [packages/legacy/treetime/treetime/treetime.py#L235-L243](../../packages/legacy/treetime/treetime/treetime.py#L235-L243)

## Erratum

Before the timetree loop, v0 removes every non-root internal node whose branch is shorter than a tenth of a mutation and whose reconstructed sequence equals the sequence of its parent:

```python
if (node.branch_length < 0.1 * self.one_mutation) and (
    self.gtr.prob_t(node.up._cseq, node._cseq, 0.0, pattern_multiplicity=...) > 0.1
):
    node.up.clades = [k for k in node.up.clades if k != node] + node.clades
```

The selection uses the branch length and the sequences alone. `ClockTree` has already read the dates at this point, so a node with a date constraint (a sampled ancestor, or a dated point on a lineage) is removed together with its date. Its children move to its parent, and no warning names the lost date.

A node without a date carries no data that the collapse loses: its sequence equals its parent's, and its branch is shorter than the resolution of the alignment. A dated node constrains the time of its point on the tree. The date cannot move to the parent or to the children, because they are other points, separated by branches of unknown duration, which is the quantity the time inference estimates.

## Evidence

- The setup code reads a date for every named node, internal nodes included, and keeps a dated internal node out of the bad branches ([packages/legacy/treetime/treetime/clock_tree.py#L121-L152](../../packages/legacy/treetime/treetime/clock_tree.py#L121-L152)). The time inference then turns the date into a constraint of that node ([L378-L387](../../packages/legacy/treetime/treetime/clock_tree.py#L378-L387)). v0 thus treats the date of an internal node as data
- Where v0 ignores a date on purpose, it warns: "Date constraint will be ignored!" ([packages/legacy/treetime/treetime/clock_tree.py#L389-L396](../../packages/legacy/treetime/treetime/clock_tree.py#L389-L396)). The pruning drops the date without a warning
- The same defect in the polytomy cleanup of v0: [timetree-polytomy-cleanup-drops-dated-single-child-node.md](timetree-polytomy-cleanup-drops-dated-single-child-node.md)

## v0 impact

- Only trees with a dated internal node are affected, which requires a date for an internal node in the metadata
- The removed node's date does not constrain the times of its parent and children, and the node is missing from the outputs

## v1 status

`fn collapse_short_internal_edges()` in [packages/treetime/src/optimize/topology/collapse.rs](../../packages/treetime/src/optimize/topology/collapse.rs) applies v0's rule, and the timetree pre-loop passes the nodes with a date constraint as nodes to keep ([packages/treetime/src/timetree/pre_loop.rs](../../packages/treetime/src/timetree/pre_loop.rs)). A dated internal node stays in the tree with its date constraint.
