# Auspice lacks hooks for apps that embed its components

## Problem

The TreeTime apps render Auspice components inside their own pages and build the Auspice sources with their own Vite configuration (`fn auspice()` in `packages/app-ui/build/auspice-vite.ts`). Auspice 3.1.0 assumes its own app and build in several places, which the Bun patch of Auspice, `patches/auspice@3.1.0.patch`, changes:

- **Extensions registry**: Auspice reads its extensions from a build-time `EXTENSION_DATA` variable and loads them through webpack's `require.context("@extensions")`, which Vite does not provide. The patch imports the registry from the `@extensions` module, which `fn auspice()` aliases to `packages/app-ui/src/auspice/extensions.ts`. `getExtension()` also returns registered values that are falsy, such as `false`. The registry holds plain data only, because Vite's dependency optimizer bundles the alias target into the pre-bundled Auspice modules, so a component or a React context in it would exist twice in the dev server
- **Figure caption**: the SVG download writes the text lines of a Nextstrain dataset below the figure: the page address, the data use terms and the CC-BY license. The patch adds a `caption` parameter to the asynchronous `SVG()` that replaces these lines, which TreeTime uses to name the run, the command and the TreeTime version
- **Download modal**: the "d" key opens Auspice's download modal when it is rendered. The patch skips the key when the `downloadModal` extension is `false`. TreeTime renders no download modal
- **Kept view**: once the tree is drawn, Auspice keeps the zoom on the drawn nodes only. The patch lets a new drawing of the same nodes take the zoom from the previous drawing

## Code

- Extensions registry: [util/extensions.ts#L8-L35](https://github.com/nextstrain/auspice/blob/v3.1.0/src/util/extensions.ts#L8-L35)
- SVG download entry: [download/helperFunctions.js#L789](https://github.com/nextstrain/auspice/blob/v3.1.0/src/components/download/helperFunctions.js#L789)
- Download modal key: [modal/Modal.jsx#L29](https://github.com/nextstrain/auspice/blob/v3.1.0/src/components/modal/Modal.jsx#L29)
- Zoom of a new drawing: [tree/phyloTree/phyloTree.ts#L39](https://github.com/nextstrain/auspice/blob/v3.1.0/src/components/tree/phyloTree/phyloTree.ts#L39)

## Fix

Contribute the hooks to [nextstrain/auspice](https://github.com/nextstrain/auspice), then remove them from the patch after upgrading to an Auspice release that contains them. [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md) describes how to edit the patch.
