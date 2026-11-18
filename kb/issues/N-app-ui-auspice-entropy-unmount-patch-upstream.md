# Auspice entropy chart is drawn twice under React StrictMode

## Problem

In Auspice 3.2.0, `componentWillUnmount()` of the entropy chart removes the chart only when the component state holds it (`if (this.state.chart)`). `setUp()` draws the chart first and stores it later with `this.setState({chart})`. All three TreeTime app shells render `StrictMode` (`packages/app-web/src/main.tsx`, `packages/app-web2/src/main.tsx`, `packages/app-desktop/renderer/main.tsx`), which mounts, unmounts and mounts each component again in development builds. At the simulated unmount the state does not hold the chart yet, so the cleanup is skipped, and the second mount draws a second chart over the first in the same SVG. A probe with the React 19.3 development build and a class component that draws in `componentDidMount()` and then calls `setState()` showed this order.

TreeTime fixes this in the Bun patch of Auspice, `patches/auspice@3.2.0.patch`: `componentWillUnmount()` stops the visibility observer and removes everything under the chart's SVG element, whatever the state holds.

## Code

- Unmount: [entropy/index.js#L337-L341](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/components/entropy/index.js#L337-L341)
- Drawing before the state update: [entropy/index.js#L230-L236](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/components/entropy/index.js#L230-L236)

## Fix

Contribute the change to [nextstrain/auspice](https://github.com/nextstrain/auspice), then remove it from the patch after upgrading to an Auspice release that contains it. [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md) describes how to edit the patch.

## Proposed implementation plan

1. Check upstream state: on the current default branch of nextstrain/auspice (clone in `.repos/auspice`), check whether `componentWillUnmount()` of `src/components/entropy/index.js` still calls `this.state.chart.destroy()` only under `if (this.state.chart)`; the local clone is at `v3.1.0`, so fetch the default branch first
2. Take the `src/components/entropy/index.js` hunk of `patches/auspice@3.2.0.patch` (`select(this.d3entropy).selectAll("*").remove()` in place of the guarded `destroy()`) as the upstream change. Alternative for upstream review: store the chart on the instance (`this.chart`) before the first draw, so that `destroy()` always finds it
3. Open a pull request at nextstrain/auspice with a reproduction: render `<StrictMode>` around the entropy panel in a development build and count the bar groups under the panel's SVG element (two sets before the change, one after it)
4. After an Auspice release contains the change: bump `auspice` in `packages/app-ui/package.json` and the `patchedDependencies` key in the root `package.json`, regenerate the patch without this hunk (procedure in [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md)), and delete this issue

Tests: none in TreeTime, because the change lives in Auspice. The oracle is React's documented StrictMode behavior (mount, unmount, mount in development builds), which the reproduction of step 3 shows.

Validation: `just typecheck`, `just test-ts`, and a manual check in the web2 dev server that the entropy panel shows one set of bars after a view switch.
