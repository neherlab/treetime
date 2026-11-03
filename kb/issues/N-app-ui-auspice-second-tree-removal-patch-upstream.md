# Auspice keeps the branch lengths of a removed second tree

## Problem

Auspice 3.2.0 offers the branch lengths that both trees have while a second tree is shown. When the user removes the second tree, the reducer of `REMOVE_TREE_TOO` keeps that list, so the main tree offers no time view or mutation count that only the second tree lacked, until the page loads again. TreeTime builds a new store for each pair of runs and never removes a second tree, but the patch keeps the same change as the other apps that embed Auspice.

TreeTime fixes this in the Bun patch of Auspice, `patches/auspice@3.2.0.patch`: the thunk `removeSecondTree()` dispatches `REMOVE_TREE_TOO` with the branch lengths of the main tree, the reducer takes them, and the second-tree select dispatches the thunk.

## Code

- Second tree removal: [reducers/controls.ts#L509](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/reducers/controls.ts#L509)
- Second-tree select: [controls/choose-second-tree.js#L63](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/components/controls/choose-second-tree.js#L63)

## Fix

Contribute the change to [nextstrain/auspice](https://github.com/nextstrain/auspice), then remove it from the patch after upgrading to an Auspice release that contains it. [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md) describes how to edit the patch.
