# Auspice entropy chart crowds its genome axes in narrow viewers

## Problem

Auspice 3.2.0 asks d3 for 16 ticks on the main genome axis of the entropy chart and 20 on the navigation axis below it, whatever the chart width. In a narrow viewer, the labels (such as "28,000") run into each other and cannot be read.

TreeTime fixes this in the Bun patch of Auspice, `patches/auspice@3.2.0.patch`: each axis asks for one tick per 60 px of its width, at most the numbers above.

## Code

- Axis ticks: [entropy/entropyD3.js#L627-L630](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/components/entropy/entropyD3.js#L627-L630)

## Fix

Contribute the change to [nextstrain/auspice](https://github.com/nextstrain/auspice), then remove it from the patch after upgrading to an Auspice release that contains it. [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md) describes how to edit the patch.

## Proposed implementation plan

1. Check upstream state: on the current default branch of nextstrain/auspice (clone in `.repos/auspice`), check whether `src/components/entropy/entropyD3.js` still calls `.ticks(16)` on `axes.xMain` and `.ticks(20, ",f")` on `axes.xNav`; the local clone is at `v3.1.0`, so fetch the default branch first
2. Take the hunk of `patches/auspice@3.2.0.patch` that adds `minPxPerAxisTick` and `tickCount(scale, max)` and use it as the upstream change, against the current upstream default branch
3. Open a pull request at nextstrain/auspice with a before and after screenshot of a narrow entropy panel (for example 400 px wide), where the labels such as "28,000" overlap before the change and are separate after it
4. After an Auspice release contains the change: bump `auspice` in `packages/app-ui/package.json` and the `patchedDependencies` key in the root `package.json`, regenerate the patch without this hunk (procedure in [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md)), and delete this issue

Tests: none in TreeTime, because the change lives in Auspice. The oracle is the visual check of step 3: at the chart width W, each axis shows at most floor(W / 60) ticks and at least 2.

Validation: `just typecheck`, `just test-ts`, and a manual check of the entropy panel at a narrow and a wide viewer width in the web app.
