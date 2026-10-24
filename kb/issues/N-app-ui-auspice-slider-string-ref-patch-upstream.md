# Auspice date range slider fails under React 19

## Problem

The slider of Auspice 3.1.0's "Date Range" control renders its track, its handles and its bars with string refs (`ref="slider"`, `ref={'handle' + i}`, `ref={'bar' + i}`). React 19 removed string refs and throws "Expected ref to be a function, an object returned by React.createRef(), or undefined/null" when it renders such an element, so the whole tree view fails as soon as the sidebar shows the date range of a dated tree. The TreeTime apps run React 19, because Auspice runs in-process with the app; Auspice itself declares React 16. Nothing reads the refs: the slider computes its positions from fixed lengths.

TreeTime fixes this in the Bun patch of Auspice, `patches/auspice@3.1.0.patch`, which removes the string refs.

## Code

- String refs: [controls/slider.js#L677](https://github.com/nextstrain/auspice/blob/v3.1.0/src/components/controls/slider.js#L677), [controls/slider.js#L735](https://github.com/nextstrain/auspice/blob/v3.1.0/src/components/controls/slider.js#L735), [controls/slider.js#L804](https://github.com/nextstrain/auspice/blob/v3.1.0/src/components/controls/slider.js#L804)

## Fix

Contribute the change to [nextstrain/auspice](https://github.com/nextstrain/auspice), then remove it from the patch after upgrading to an Auspice release that contains it. [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md) describes how to edit the patch.
