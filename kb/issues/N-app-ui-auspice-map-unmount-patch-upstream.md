# Auspice map fails when it is mounted again

## Problem

The tree view mounts the Auspice 3.2.0 map when the user turns it on in the sidebar, and unmounts it when the user turns it off or leaves the result. The map component assumes that it is mounted once per page load:

- **Stale export**: the map sets `window.L.getMapSvg`, the function that the SVG download uses, only once, and the function keeps the first map instance. After the map is mounted again, the download reads the detached first map and writes an empty map
- **Crash on early unmount**: the Leaflet map is created on the first update, but `componentWillUnmount()` reads it without a check, so an unmount before that update (for example, React StrictMode's immediate remount, or a failed WebGL setup) throws `Cannot read properties of null (reading 'off')`
- **Leaked resources**: `componentWillUnmount()` removes only two event handlers. It does not remove the Leaflet map, which holds a WebGL context, nor its timers

TreeTime fixes this in the Bun patch of Auspice, `patches/auspice@3.2.0.patch`: each mount registers the export function of its own instance and removes it on unmount when it still belongs to that instance; the cleanup after an unmount removes the Leaflet map and clears the timers, works when the map was never created, and runs after the current task, so that React StrictMode's immediate remount keeps the map.

## Code

- Export function: [map/map.js#L105-L139](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/components/map/map.js#L105-L139)
- Unmount: [map/map.js#L736-L739](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/components/map/map.js#L736-L739)

## Fix

Contribute the changes to [nextstrain/auspice](https://github.com/nextstrain/auspice), then remove them from the patch after upgrading to an Auspice release that contains them. [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md) describes how to edit the patch.

## Proposed implementation plan

1. Check upstream state: on the current default branch of nextstrain/auspice (clone in `.repos/auspice`), check whether `src/components/map/map.js` still sets `window.L.getMapSvg` only when it is unset, and `componentWillUnmount()` still calls `this.state.map.off(...)` without a check and removes neither the map nor its timers; the local clone is at `v3.1.0`, so fetch the default branch first
2. Take the `src/components/map/map.js` hunks of `patches/auspice@3.2.0.patch` as the upstream change, split in three commits that match the three problems: per-instance `this.getMapSvg` registered in `componentDidMount()` and deleted on unmount when `window.L.getMapSvg === this.getMapSvg`; a null check of `this.state.map` in `componentWillUnmount()`; the deferred cleanup (`this.pendingCleanup = window.setTimeout(...)`) that removes the Leaflet map and clears the timers, cancelled by a remount
3. Open a pull request at nextstrain/auspice with reproductions: toggle the map panel off and on, then download the SVG (empty map before, full map after); render the map under `<StrictMode>` in a development build (TypeError before, none after); count WebGL contexts after ten toggles (growing before, constant after)
4. After an Auspice release contains the change: bump `auspice` in `packages/app-ui/package.json` and the `patchedDependencies` key in the root `package.json`, regenerate the patch without these hunks (procedure in [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md)), and delete this issue

Tests: none in TreeTime, because the change lives in Auspice. The oracles are the observable results of step 3.

Validation: `just typecheck`, `just test-ts`, and a manual check in the web2 dev server: toggle the map off and on, then download the SVG and check that it contains the map.
