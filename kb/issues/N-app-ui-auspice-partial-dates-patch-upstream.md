# Auspice offers a time axis for trees where only some nodes have dates

## Problem

Auspice 3.1.0 offers the time metric, the clock layout and the date controls when the root node of the tree has a date (`num_date`). In a tree where some nodes have no date, Auspice offers the time view and leaves these nodes out of the drawing, because their positions on the time axis are not numbers. Augur trees and trees of other tools can have such nodes, and the compare page draws a second tree that may lack dates where the first has them.

TreeTime fixes this in the Bun patch of Auspice, `patches/auspice@3.1.0.patch`: Auspice offers the time axis only when every node of the main tree and of the second tree has a date. Trees with undated nodes are shown by divergence, with no time metric, clock layout or date controls.

## Code

- Time axis availability: [actions/recomputeReduxState.js#L476-L480](https://github.com/nextstrain/auspice/blob/v3.1.0/src/actions/recomputeReduxState.js#L476-L480)

## Fix

Contribute the change to [nextstrain/auspice](https://github.com/nextstrain/auspice), then remove it from the patch after upgrading to an Auspice release that contains it. [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md) describes how to edit the patch.
