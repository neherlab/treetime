# Auspice animates the tree when its size changes

## Problem

When the SVG of an Auspice 3.2.0 tree changes size, `change()` maps the tree to the new size with the usual transition. The SVG takes its new size at once, so the transition starts from positions computed for the old size and draws the tree out of place until it ends. The tree view of the TreeTime apps changes the tree size whenever the window or the sidebar changes width, the panel layout changes, or the tree header moves to its own row.

The Bun patch of Auspice, `patches/auspice@3.2.0.patch`, sets the transition time to 0 when the SVG size has changed.

## Code

- Transition time of an update: [tree/phyloTree/change.ts#L342-L345](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/components/tree/phyloTree/change.ts#L342-L345)

## Fix

Contribute the change to [nextstrain/auspice](https://github.com/nextstrain/auspice), then remove it from the patch after upgrading to an Auspice release that contains it. [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md) describes how to edit the patch.
