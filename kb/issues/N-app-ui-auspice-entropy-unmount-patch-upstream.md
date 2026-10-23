# Auspice entropy chart is drawn twice under React StrictMode

## Problem

In Auspice 3.1.0, `componentWillUnmount()` of the entropy chart removes the chart only when the component state holds it (`if (this.state.chart)`). `setUp()` draws the chart first and stores it later with `this.setState({chart})`. All three TreeTime app shells render `StrictMode` (`packages/app-web/src/main.tsx`, `packages/app-web2/src/main.tsx`, `packages/app-desktop/renderer/main.tsx`), which mounts, unmounts and mounts each component again in development builds. At the simulated unmount the state does not hold the chart yet, so the cleanup is skipped, and the second mount draws a second chart over the first in the same SVG. A probe with the React 19.3 development build and a class component that draws in `componentDidMount()` and then calls `setState()` showed this order.

TreeTime fixes this in the Bun patch of Auspice, `patches/auspice@3.1.0.patch`: `componentWillUnmount()` stops the visibility observer and removes everything under the chart's SVG element, whatever the state holds.

## Code

- Unmount: [entropy/index.js#L337-L341](https://github.com/nextstrain/auspice/blob/v3.1.0/src/components/entropy/index.js#L337-L341)
- Drawing before the state update: [entropy/index.js#L230-L236](https://github.com/nextstrain/auspice/blob/v3.1.0/src/components/entropy/index.js#L230-L236)

## Fix

Contribute the change to [nextstrain/auspice](https://github.com/nextstrain/auspice), then remove it from the patch after upgrading to an Auspice release that contains it. [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md) describes how to edit the patch.
