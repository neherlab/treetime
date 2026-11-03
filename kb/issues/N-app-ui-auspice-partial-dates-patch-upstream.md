# Auspice offers a time axis for trees where only some nodes have dates

## Problem

Auspice 3.2.0 offers the time branch length, the clock layout and the date controls when the root node of the tree has a date (`num_date`). In a tree where some nodes have no date, Auspice offers the time view and leaves these nodes out of the drawing, because their positions on the time axis are not numbers. Augur trees and trees of other tools can have such nodes, and the compare page draws a second tree that may lack dates where the first has them.

TreeTime fixes this in the Bun patch of Auspice, `patches/auspice@3.2.0.patch`: `findAvailableDistanceMeasures()` offers `num_date` only when every node of the tree has a date. Auspice offers the branch lengths that both trees have when a second tree is shown, so a second tree with undated nodes also removes the time view. Trees with undated nodes are shown by divergence or mutation counts, with no time branch length, clock layout or date controls.

## Code

- Time branch length availability: [util/computeDistanceMetric.ts#L101-L108](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/util/computeDistanceMetric.ts#L101-L108)

## Fix

Contribute the change to [nextstrain/auspice](https://github.com/nextstrain/auspice), then remove it from the patch after upgrading to an Auspice release that contains it. [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md) describes how to edit the patch.
