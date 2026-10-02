# Tree preparation marks a dated internal node as a bad branch and drops its date

## v0 location

`TreeAnc._prepare_nodes()` (`#TreeAnc`, `#_prepare_nodes`) [packages/legacy/treetime/treetime/treeanc.py#L484-L488](../../packages/legacy/treetime/treetime/treeanc.py#L484-L488)

## Erratum

v0 sets the `bad_branch` flag of internal nodes in two places, with two different rules.

Rule A, at setup, in `ClockTree._assign_dates()` [packages/legacy/treetime/treetime/clock_tree.py#L121-L152](../../packages/legacy/treetime/treetime/clock_tree.py#L121-L152): a node with a date constraint is never bad. An internal node without a date is bad when all its children are bad.

Rule B, on every tree preparation, in `TreeAnc._prepare_nodes()`:

```python
for clade in self.tree.find_clades(order='postorder'):  # children first
    if clade.is_terminal():
        clade.bad_branch = clade.bad_branch if hasattr(clade, 'bad_branch') else False
    else:
        clade.bad_branch = all([c.bad_branch for c in clade])
```

An internal node is bad when all its children are bad. Its own date is not consulted.

`prepare_tree()` runs after the clock filter ([packages/legacy/treetime/treetime/treetime.py#L509](../../packages/legacy/treetime/treetime/treetime.py#L509)), inside `reroot()` ([L642](../../packages/legacy/treetime/treetime/treetime.py#L642)) and after polytomy resolution ([L330](../../packages/legacy/treetime/treetime/treetime.py#L330)). In a default run, rule B therefore overwrites rule A before the first time inference.

The backward pass then drops all evidence of a bad node, its own date included ([packages/legacy/treetime/treetime/clock_tree.py#L685-L687](../../packages/legacy/treetime/treetime/clock_tree.py#L685-L687)):

```python
if node.bad_branch:
    # no information
    node.marginal_pos_Lx = None
```

An internal node that the user dated, but whose leaves are all undated or all clock-filter outliers, is thus treated as if it had no date. The date does not constrain the node's own time and does not reach its parent.

## Evidence

- The adjacent setup code `_assign_dates()` applies rule A and states it as the intent: "If all branches dowstream are 'bad', and there is no date constraint for this node, the branch is marked as 'bad'" ([packages/legacy/treetime/treetime/clock_tree.py#L150-L152](../../packages/legacy/treetime/treetime/clock_tree.py#L150-L152))
- A bad branch means "no information" in the backward pass ([L686](../../packages/legacy/treetime/treetime/clock_tree.py#L686)), but a node with its own date carries information whatever its children do
- Rule B discards a date the user supplied, without a warning

## v0 impact

- Only internal nodes with a date constraint whose subtree has no usable leaf date are affected. Such nodes exist only when the user supplies a date for an internal node
- The affected node is dated from its parent alone, and its date does not contribute to the times of its ancestors

## v1 status

v1 applies rule A for every time inference. [`derive_bad_branches()`](../../packages/treetime/src/timetree/inference/bad_branches.rs) derives the flags for the current topology from the leaf flags: a node with its own date constraint is never bad. See [kb/decisions/timetree-dated-internal-node-never-bad-branch.md](../decisions/timetree-dated-internal-node-never-bad-branch.md).
