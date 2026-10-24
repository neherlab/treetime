# Rooting on a named node roots at the middle of its unrooted branch

`--outgroup` and `--reroot=oldest` root the tree at the middle of the branch above the named node in the unrooted tree. The named node is the most recent common ancestor of the outgroup tips, or the oldest dated tip. A second reroot of the result changes nothing, also when a user reroots the output of an earlier run.

**Type**: Bug fix (v0 erratum correction). Approved by the user on 2026-10-24.

**v0 erratum**: [kb/v0-errata/reroot-named-node-halves-branch-on-each-call.md](../v0-errata/reroot-named-node-halves-branch-on-each-call.md).

**v1 location**: `fn named_node_point()` in `packages/treetime/src/reroot/named_node.rs`, and the placement modes of `fn RootTarget::on_edge()` in `packages/treetime/src/reroot/placement.rs`.

## Decision

Let $M$ be the named node and $R$ its parent:

- When $R$ is not the root, or $R$ does not have exactly two children, the root splits the edge `R -> M` in the middle. The branch of $M$ in the unrooted tree is then the edge `R -> M`
- When $R$ is the root with exactly the two children $M$ and $S$, the branch of $M$ in the unrooted tree is `R -> M` plus `R -> S`, of length $a + b$, where $a$ and $b$ are the lengths of the two edges. Its middle lies $(a+b)/2$ from $M$:
  - $a > b$: on the edge `R -> M` at the split $x = (a-b)/(2a)$ from $R$
  - $b > a$: on the edge `R -> S` at the split $x = (b-a)/(2b)$ from $R$
  - $a = b$, or a length is missing or invalid: the root stays at $R$

After the move, the old root $R$ has one parent and one child, and the reroot merges its two edges.

A root position that is computed from branch lengths is placed exactly. It moves to an edge end only when its distance from that end, $x \ell$ or $(1-x)\ell$ for an edge of length $\ell$, is zero within 4 ulps of $\ell$. The least-squares and `min-dev` root searches keep their edge-end band of $10^{-5}$ of the edge length, which keeps the root when a later search finds the same root after rounding ([kb/algo/reroot.md](../algo/reroot.md)). In both modes a leaf never becomes the root. Biopython, and so augur, also place an outgroup root exactly.

## v0 divergence

TreeTime v0 places the root at half of the named node's current branch on every call of `reroot()`. When the node is a child of a two-child root, this halves its branch again and moves the root toward the node. v1 roots at the middle of the unrooted branch, which is the root v0 gives on a tree whose root is not next to the node.

## Outputs this changes

Every reroot on a named node where the node is a child of a two-child root at that moment:

- an input tree that is already rooted next to the node, as augur and IQ-TREE trees with a two-child root often are: one reroot is enough
- a run that reroots more than once: `--reroot=oldest` or `--outgroup` in `clock` with the outlier filter, and in `timetree`
- a reroot of the output of an earlier run

A named node whose branch has length zero, or no length, now roots at its parent instead of at a new node with two branches of length zero.

## Validation contract

- Two reroots on the same outgroup give the same topology and branch lengths as one reroot
- An outgroup next to a two-child root with branches $a$ and $b$ gets $(a+b)/2$ on both sides of the root, also when the split falls inside the edge-end band
- `clock --reroot=oldest` gives the same root branch lengths with the outlier filter on as with it off
