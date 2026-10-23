# Auspice tree clip areas share one element ID

## Problem

In Auspice 3.1.0, every tree object creates a `<clipPath>` with the fixed ID `treeClip`, and all branch, tip, label and grid groups clip to `url(#treeClip)`. Element IDs are unique per document, and the browser resolves a duplicated ID to the first element. The compare page draws the trees of two runs side by side, so both trees would clip to the clip area of the tree drawn first. The SVG download would also write both clip areas under one ID.

The legend has the same problem: every legend clips its first column to the clip area `legendFirstColumnClip`. The tree and the map each show a legend, so the second legend would clip to the area of the first.

TreeTime fixes this in the Bun patch of Auspice, `patches/auspice@3.1.0.patch`, registered in `patchedDependencies` of the root `package.json`: the clip ID of a tree comes from its tree ID (`treeClip-LEFT`, `treeClip-RIGHT`), and each legend instance gets its own clip ID.

## Code

- Clip area creation: [tree/phyloTree/renderers.ts#L524](https://github.com/nextstrain/auspice/blob/v3.1.0/src/components/tree/phyloTree/renderers.ts#L524)
- Legend clip area: [legend/legend.tsx#L220](https://github.com/nextstrain/auspice/blob/v3.1.0/src/components/legend/legend.tsx#L220), [legend/legend.tsx#L234](https://github.com/nextstrain/auspice/blob/v3.1.0/src/components/legend/legend.tsx#L234)
- Clip references: [tree/phyloTree/renderers.ts#L224](https://github.com/nextstrain/auspice/blob/v3.1.0/src/components/tree/phyloTree/renderers.ts#L224), [tree/phyloTree/labels.js#L138](https://github.com/nextstrain/auspice/blob/v3.1.0/src/components/tree/phyloTree/labels.js#L138), [tree/phyloTree/grid.js#L36-L39](https://github.com/nextstrain/auspice/blob/v3.1.0/src/components/tree/phyloTree/grid.js#L36-L39)

## Fix

Contribute the change to [nextstrain/auspice](https://github.com/nextstrain/auspice), then remove it from the patch after upgrading to an Auspice release that contains it. [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md) describes how to edit the patch.
