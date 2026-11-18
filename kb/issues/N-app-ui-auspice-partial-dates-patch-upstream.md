# Auspice offers a time axis for trees where only some nodes have dates

## Problem

Auspice 3.2.0 offers the time branch length, the clock layout and the date controls when the root node of the tree has a date (`num_date`). In a tree where some nodes have no date, Auspice offers the time view and leaves these nodes out of the drawing, because their positions on the time axis are not numbers. Augur trees and trees of other tools can have such nodes, and the compare page draws a second tree that may lack dates where the first has them.

TreeTime fixes this in the Bun patch of Auspice, `patches/auspice@3.2.0.patch`: `findAvailableDistanceMeasures()` offers `num_date` only when every node of the tree has a date. Auspice offers the branch lengths that both trees have when a second tree is shown, so a second tree with undated nodes also removes the time view. Trees with undated nodes are shown by divergence or mutation counts, with no time branch length, clock layout or date controls.

## Code

- Time branch length availability: [util/computeDistanceMetric.ts#L101-L108](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/util/computeDistanceMetric.ts#L101-L108)

## Fix

Contribute the change to [nextstrain/auspice](https://github.com/nextstrain/auspice), then remove it from the patch after upgrading to an Auspice release that contains it. [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md) describes how to edit the patch.

## Proposed implementation plan

1. Check upstream state: on the current default branch of nextstrain/auspice (clone in `.repos/auspice`), check whether `findAvailableDistanceMeasures()` in `src/util/computeDistanceMetric.ts` still offers `num_date` from the root node alone (`numDateAtRoot`)
2. Take the `src/util/computeDistanceMetric.ts` hunk of `patches/auspice@3.2.0.patch` (`allNodesHaveNumDate = nodes.every(...)`) as the upstream change
3. Open a pull request at nextstrain/auspice with a unit test in its Jest suite: a tree whose root and some tips have `num_date` and one tip lacks it yields no `num_date` among the available measures, and a fully dated tree still yields it. Oracle: the time axis needs a numeric date for every drawn node, so a node without `num_date` has no position on it
4. After an Auspice release contains the change: bump `auspice` in `packages/app-ui/package.json` and the `patchedDependencies` key in the root `package.json`, regenerate the patch without this hunk (procedure in [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md)), and delete this issue

Validation: `just typecheck`, `just test-ts`, and a manual check in the web2 dev server that a tree with an undated tip opens by divergence and offers no time controls.
