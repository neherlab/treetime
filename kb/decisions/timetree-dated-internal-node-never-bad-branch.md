# A dated internal node is never a bad branch

v1 derives the bad-branch flags anew for every time inference, and a node with its own date constraint is never bad. v0 applies this rule only at setup and then overwrites it: after the clock filter, reroot and polytomy resolution, v0 marks an internal node bad whenever all its children are bad, even when the node has its own date.

**Type**: Bug fix (v0 erratum correction).

**v0 location**: `TreeAnc._prepare_nodes()` at `packages/legacy/treetime/treetime/treeanc.py#L484-L488`, which overwrites the rule of `ClockTree._assign_dates()` at `packages/legacy/treetime/treetime/clock_tree.py#L121-L152`.

**v1 location**: [`derive_bad_branches()`](../../packages/treetime/src/timetree/inference/bad_branches.rs), called by [`run_timetree()`](../../packages/treetime/src/timetree/inference/runner.rs) at the start of each time inference.

**v0 erratum**: [kb/v0-errata/timetree-bad-branch-ignores-internal-date.md](../v0-errata/timetree-bad-branch-ignores-internal-date.md).

## Rule

The inputs are the leaf flags: a leaf is bad when it has no date, or when the clock filter marked it as an outlier. From these, each time inference derives the flags of the current tree in postorder:

- A leaf is bad per its leaf flag
- An internal node with its own date constraint is not bad
- Any other internal node is bad when all its children are bad

A bad node sends no message to its parent. It still receives a posterior: when it has no subtree evidence, its posterior is the message from its parent alone, as in v0 (`clock_tree.py#L890-L891`).

## Why

A bad branch means that the subtree carries no information about time. A node with a date carries information whatever its children do. v0 states this intent in `_assign_dates()`, and rule B silently discards a date the user supplied.

Deriving the flags per inference, instead of carrying them from round to round, also keeps them correct for the current topology after reroot and polytomy resolution. The leaf flags are derived once from the tree after the last reroot, because polytomy resolution adds only internal nodes and the refinement rounds do not reroot.

## Impact

- Differs from v0 only for internal nodes with a date constraint whose subtree has no usable leaf date
- Such a node keeps its given date, and the date contributes to the times of its ancestors
