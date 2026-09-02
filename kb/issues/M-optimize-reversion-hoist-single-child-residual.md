# Reversion hoist leaves a residual reversion when every child reverts the same position

## Symptom and reproduction

In the optimize loop, `resolve_polytomies` resolves a reversion polytomy by merging siblings that share substitutions, then hoisting the reverting group under a new node. When _every_ child of a polytomy node reverts the same parent-edge substitution, the shared-mutation merge groups all of them under one helper node, leaving the polytomy node with a single child. The hoist then declines to fire, because moving that single reverting child under a new node would leave the polytomy node childless (a spurious leaf).

The result keeps one residual reversion on the merged child edge. The total mutation count still drops relative to the input (the merge removes the shared reversion from every duplicate), but it does not reach the parsimony optimum, which would remove the reversion entirely.

Reproduction: a node `V` under parent edge `U -> V = {A0C}` whose children all carry `C0A`. After the merge, `V` has a single child (the helper) carrying `C0A`; the optimum is zero mutations (drop `A0C` and the reversion), but the routine stops at two (`A0C` on `U -> V`, `C0A` on the helper edge).

## Impact and scope

Narrow. It affects only polytomies where all children revert the same position, and only leaves a reversion that the merge already reduced from many to one. The output is never worse than the input and never invalid. Correctness (no added mutations, distance preservation, single-rooted tree) is unaffected.

## Root cause

`fn try_hoist_reverting_child()` requires the node to keep at least two children so the hoisted child leaves a sibling behind. A single-child node drops below that threshold. The optimal resolution here is to collapse the now-degree-two node into its parent (composing the two edges cancels the reversion), but collapsing a pre-existing node contradicts the routine's helper-retirement scope, which deliberately never dissolves nodes the input asserted.

## Fix approach

Recognize a degree-two node whose single child reverts its parent edge and collapse that node's parent edge (compose `M_v` with `M_c`, cancelling the reversion). This removes a topologically vacuous node (a single-child node carries no phylogenetic split), so it preserves all bipartitions. Requires distinguishing this case from the helper-retirement scope, and a decision on whether removing pre-existing degree-two nodes during `optimize` is acceptable.

This is a greedy limitation in the same family as the incompatible-splits case the routine already accepts (disagreeing reverting children leave irreducible homoplasy); exact resolution is subtree parsimony re-optimization, which is out of scope.

## Locations

- `fn try_hoist_reverting_child()` two-children guard: `packages/treetime/src/optimize/topology/resolve_polytomy.rs`
- Hoist move: `packages/treetime/src/optimize/topology/hoist_reversions.rs`
