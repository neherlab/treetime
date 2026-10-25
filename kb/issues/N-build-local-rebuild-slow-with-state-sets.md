# Local rebuilds became 2 to 4 times slower with state-set queries

Since queries carry state sets, every local rebuild computes the plurality sets of each internal pseudo-leaf and encodes a new query for it. This happens for every ball, although neighbouring balls share most of their pseudo-leaves. Rebuilds of radius 3 take 85 s on `sc2/4500`, against 26 s before state sets ([kb/reports/feat-tree-infer-prototype.md](../reports/feat-tree-infer-prototype.md)).

## Location

`fn best_local_tree()`: [`packages/treetime/src/build/rebuild.rs#L130-L194`](../../packages/treetime/src/build/rebuild.rs#L130-L194), which calls `fn Mat.plurality_sets()` per pseudo-leaf.

## Fix

- Cache the plurality sets per node, and invalidate them only on splice
- In later rounds, rebuild only the balls that contain a node changed since the previous round
