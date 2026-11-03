# Auspice SVG download cuts tree labels, misses off-screen charts and shifts text after links

## Problem

The SVG download of Auspice 3.2.0, which the "Download" menu of the tree view offers as "Figure (SVG)", writes each panel as a nested `<svg>` of the panel's size, with `overflow="hidden"` and a `<clipPath>` of the panel's box:

- Tree labels that extend past the tree panel, such as long tip labels or axis labels at its edge, are cut off in the downloaded image, although the page shows them
- The entropy chart computes its bars only while it is on screen. Below the fold, the downloaded image shows the chart without bars
- Text that follows an inline element on the same line, such as the words after the link of the map attribution, starts one space too far left, because the line keeps the space before it and an SVG `<text>` drops leading spaces

TreeTime fixes this in the Bun patch of Auspice, `patches/auspice@3.2.0.patch`: the download measures the labels of each tree panel and records the box of the panel and its labels, drops labels that are anchored outside the panel and never visible on screen, gives such a panel `overflow="visible"` and a clip path of that box (the nested viewport clips at the panel edge even when the clip path is larger), widens the image with a `viewBox` that holds the boxes, computes the entropy bars before writing the image when the chart has none, and starts each line of HTML text at its first visible character.

## Code

- Panel serialization: [download/helperFunctions.js#L429](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/components/download/helperFunctions.js#L429)
- Panel clipping: [download/helperFunctions.js#L570-L582](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/components/download/helperFunctions.js#L570-L582)
- Text lines of HTML blocks: [download/domToSvg.js#L106](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/components/download/domToSvg.js#L106)
- SVG download entry: [download/helperFunctions.js#L791](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/components/download/helperFunctions.js#L791)

## Fix

Contribute the changes to [nextstrain/auspice](https://github.com/nextstrain/auspice), then remove them from the patch after upgrading to an Auspice release that contains them. [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md) describes how to edit the patch.
