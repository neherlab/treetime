# Auspice tree clip areas share one element ID

## Problem

In Auspice 3.2.0, every tree object creates a `<clipPath>` with the fixed ID `treeClip`, and all branch, tip, label and grid groups clip to `url(#treeClip)`. Element IDs are unique per document, and the browser resolves a duplicated ID to the first element. The compare page draws the trees of two runs side by side, so both trees would clip to the clip area of the tree drawn first. The SVG download would also write both clip areas under one ID.

The legend has the same problem: every legend clips its first column to the clip area `legendFirstColumnClip`. The tree and the map each show a legend, so the second legend would clip to the area of the first.

TreeTime fixes this in the Bun patch of Auspice, `patches/auspice@3.2.0.patch`, registered in `patchedDependencies` of the root `package.json`: the clip ID of a tree comes from its tree ID (`treeClip-LEFT`, `treeClip-RIGHT`), and each legend instance gets its own clip ID.

## Code

- Clip area creation: [tree/phyloTree/renderers.ts#L524](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/components/tree/phyloTree/renderers.ts#L524)
- Legend clip area: [legend/legend.tsx#L220](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/components/legend/legend.tsx#L220), [legend/legend.tsx#L234](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/components/legend/legend.tsx#L234)
- Clip references: [tree/phyloTree/renderers.ts#L224](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/components/tree/phyloTree/renderers.ts#L224), [tree/phyloTree/labels.js#L138](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/components/tree/phyloTree/labels.js#L138), [tree/phyloTree/grid.js#L37-L40](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/components/tree/phyloTree/grid.js#L37-L40)

## Fix

Contribute the change to [nextstrain/auspice](https://github.com/nextstrain/auspice), then remove it from the patch after upgrading to an Auspice release that contains it. [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md) describes how to edit the patch.

## Proposed implementation plan

1. Open a pull request on `nextstrain/auspice` with the clip ID hunks of `patches/auspice@3.2.0.patch`: `renderers.ts` (clip area creation and the `tips`, `branchTee`, `branchStem`, and `regression` groups), `labels.js` (`branchLabels`), `grid.js` (`majorGrid`, `minorGrid`) use `treeClip-${id}`, and `legend.tsx` gives each legend instance its own first-column clip ID from an instance counter. Auspice 3.2.0 (tag `v3.2.0`, commit `db211cc6`), the release that the patch targets, still uses the fixed IDs `treeClip` and `legendFirstColumnClip`. Check first whether upstream `master` changed them after that release.
2. Upstream test: render two `PhyloTree` objects with the IDs `LEFT` and `RIGHT` into one document, and assert that the document holds two `<clipPath>` elements with distinct IDs and that each tree's groups reference the clip path of their own tree. Oracle: the SVG rule that an element ID is unique in a document.
3. After an Auspice release contains the fix, upgrade `auspice` in `package.json`, regenerate the Bun patch without these hunks (procedure in [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md)), and delete this issue.
4. Validate: `just typecheck`, `just test-ts`, and the compare page with two runs: each tree clips to its own area, and the downloaded SVG holds unique clip IDs.
