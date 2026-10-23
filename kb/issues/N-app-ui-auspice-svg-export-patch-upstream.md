# Auspice SVG download cuts tree labels and misses off-screen charts

## Problem

The SVG download of Auspice 3.1.0, which the "Download" menu of the tree view offers as "Figure (SVG)", writes each panel as a nested `<svg>` of the panel's size, with `overflow="hidden"` and a `<clipPath>` of the panel's box:

- Tree labels that extend past the tree panel, such as long tip labels or axis labels at its edge, are cut off in the downloaded image, although the page shows them
- The entropy chart computes its bars only while it is on screen. Below the fold, the downloaded image shows the chart without bars

TreeTime fixes this in the Bun patch of Auspice, `patches/auspice@3.1.0.patch`: the download measures the labels of each tree panel and records the box of the panel and its labels, drops labels that are anchored outside the panel and never visible on screen, gives such a panel `overflow="visible"` and a clip path of that box (the nested viewport clips at the panel edge even when the clip path is larger), widens the image with a `viewBox` that holds the boxes, and computes the entropy bars before writing the image when the chart has none.

## Code

- Panel serialization: [download/helperFunctions.js#L427](https://github.com/nextstrain/auspice/blob/v3.1.0/src/components/download/helperFunctions.js#L427)
- Panel clipping: [download/helperFunctions.js#L568-L580](https://github.com/nextstrain/auspice/blob/v3.1.0/src/components/download/helperFunctions.js#L568-L580)
- SVG download entry: [download/helperFunctions.js#L789](https://github.com/nextstrain/auspice/blob/v3.1.0/src/components/download/helperFunctions.js#L789)

## Fix

Contribute the changes to [nextstrain/auspice](https://github.com/nextstrain/auspice), then remove them from the patch after upgrading to an Auspice release that contains them. [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md) describes how to edit the patch.
