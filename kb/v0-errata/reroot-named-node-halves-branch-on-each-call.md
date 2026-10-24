# Rooting on a named node halves its branch on each reroot

## Summary

TreeTime v0 roots on a named node (a tip name, a list of tips, or `oldest`) by placing the root on the branch above that node, at half of the node's current branch length. When the node is already a child of a two-child root, as after any earlier reroot on it, its current branch is only one half of the branch it has in the unrooted tree. Each further reroot then halves the node's branch again and moves the root toward the node. The root therefore depends on how many times a run reroots, and on whether the input tree is already rooted next to the node.

TreeTime v1 corrects this defect: rooting on a named node places the root at the middle of the node's branch in the unrooted tree, so a repeated reroot changes nothing.

## v0 code

`def TreeTime.reroot()` resolves the named node and calls Biopython's `root_with_outgroup()` with half of the node's current branch length [packages/legacy/treetime/treetime/treetime.py#L604-L620](../../packages/legacy/treetime/treetime/treetime.py#L604-L620):

```python
self.tree.root_with_outgroup(new_root, outgroup_branch_length=new_root.branch_length / 2)
```

`new_root.branch_length` is the length of the branch between the node and its current parent. When that parent is a two-child root, the node's branch in the unrooted tree also contains the branch from the root to its sibling, of length $b$. With $a$ the length of the node's own branch, the root then lies $a/2$ from the node instead of $(a+b)/2$.

A run of `treetime clock` or `treetime timetree` calls `reroot()` more than once: `def TreeTime._run()` and `def TreeTime.clock_filter()` reroot before and after the outlier filter [packages/legacy/treetime/treetime/treetime.py#L260](../../packages/legacy/treetime/treetime/treetime.py#L260), [#L279](../../packages/legacy/treetime/treetime/treetime.py#L279), [#L487](../../packages/legacy/treetime/treetime/treetime.py#L487) and [#L513](../../packages/legacy/treetime/treetime/treetime.py#L513). After the first call the node is a child of the two-child root, so every later call halves its branch again.

## Evidence

- **Runtime**: on `data/zika/20` with outgroup `Bahia12|KX101067|2015-05-15|brazil`, whose input branch is $L = 2.76761 \cdot 10^{-4}$, `treetime clock --reroot <tip> --name-column '#name' --date-column date --sequence-length 10000` gives a tip branch of $1.3838 \cdot 10^{-4}$ ($L/2$) with `--clock-filter 0`, and $3.4595 \cdot 10^{-5}$ ($L/8$) with `--clock-filter 3`. The input and the requested root are the same; only the number of reroot calls differs
- **Docstring**: `reroot()` says that it reroots "to the node" with the given name or "to the MRCA" of the given nodes [packages/legacy/treetime/treetime/treetime.py#L558-L560](../../packages/legacy/treetime/treetime/treetime.py#L558-L560). It says nothing about earlier calls, and a root that moves on each call is not a root "to the node"
- **Code comment**: the comment above the call gives the intent, a bifurcating root on the node's branch [packages/legacy/treetime/treetime/treetime.py#L617-L618](../../packages/legacy/treetime/treetime/treetime.py#L617-L618). Halving the node's own half of an already bifurcating root does not serve this intent

## v1 status

Corrected. For a named node $M$ with parent $R$:

- when $R$ is not the root, or $R$ does not have exactly two children, the root splits the edge `R -> M` in the middle, as v0 does on the first call
- when $R$ is the root with the two children $M$ (branch length $a$) and $S$ (branch length $b$), the root lies $(a+b)/2$ from $M$: on `R -> M` at the split $(a-b)/(2a)$ when $a > b$, on `R -> S` at the split $(b-a)/(2b)$ when $b > a$, and at $R$ when $a = b$ or a length is missing

The root position is placed exactly, without the edge-end band of the least-squares root search. The approved v1 behavior is recorded in [kb/decisions/reroot-named-node-at-unrooted-branch-middle.md](../decisions/reroot-named-node-at-unrooted-branch-middle.md).
