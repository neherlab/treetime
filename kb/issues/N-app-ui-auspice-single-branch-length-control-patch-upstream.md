# Auspice shows a Branch Length select with one option

## Problem

Auspice 3.2.0 shows the Branch Length select of the sidebar for every tree in the rectangular, radial and unrooted layouts. A tree with one branch length, for example a TreeTime run without dates and without sequences, gets a select whose only option is "divergence". Auspice 3.1.0 hid the control unless the tree had both divergence and dates.

TreeTime fixes this in the Bun patch of Auspice, `patches/auspice@3.2.0.patch`: `ChooseMetric` renders nothing when the trees offer fewer than two branch lengths.

## Code

- Branch Length control: [controls/choose-metric.tsx#L64-L66](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/components/controls/choose-metric.tsx#L64-L66)

## Fix

Contribute the change to [nextstrain/auspice](https://github.com/nextstrain/auspice), then remove it from the patch after upgrading to an Auspice release that contains it. [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md) describes how to edit the patch.
