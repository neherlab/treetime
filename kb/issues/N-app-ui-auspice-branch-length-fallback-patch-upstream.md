# Auspice keeps a branch length that the shown trees do not have

## Problem

Auspice 3.2.0 computes the branch lengths that the trees offer (`availableDistanceMeasures`), but keeps the selected branch length (`distanceMeasure`, by default the time) when the trees do not offer it:

- Each load of a tree without dates starts with the time, and only the final error check switches it to the first branch length the tree has, with the console error "Requested distance measure 'num_date' is not valid"
- A second tree without dates, opened next to a tree in the time view, gets no positions: its nodes have no value for the selected branch length, and the Branch Length select shows no selection. The compare page of TreeTime shows two runs this way

TreeTime fixes this in the Bun patch of Auspice, `patches/auspice@3.2.0.patch`: `modifyControlsStateViaTree()` switches to the first branch length the trees offer when they do not offer the selected one, and `createTreeTooState()` computes the branch lengths of the main tree again, because the second tree can change the selection.

## Code

- Branch lengths of the trees: [actions/recomputeReduxState.js#L482-L490](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/actions/recomputeReduxState.js#L482-L490)
- Branch lengths after a second tree: [actions/recomputeReduxState.js#L1158](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/actions/recomputeReduxState.js#L1158)

## Fix

Contribute the change to [nextstrain/auspice](https://github.com/nextstrain/auspice), then remove it from the patch after upgrading to an Auspice release that contains it. [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md) describes how to edit the patch.
