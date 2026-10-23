# Auspice reads the coordinates of dropped metadata wrongly

## Problem

Auspice 3.1.0 turns a dropped CSV or TSV file with `latitude` and `longitude` (or `__latitude` and `__longitude`) columns into a map of the samples. The TreeTime apps do not offer dropping metadata onto the tree, so these defects do not show in TreeTime. The patch carries the fix because the Auspice patches of Nextclade, treeknit-rs and auspice-vscode carry it with the same code. Its parser has these defects:

- **Coordinates-only files**: a file whose only columns besides the names are the coordinates is refused with "No valid columns", because the check for other columns runs before the coordinates are read
- **Merged locations**: samples are grouped by the joined strings of their coordinates, so `(1, 23)` and `(12, 3)` are one location
- **Empty cells**: an empty latitude or longitude is read as 0, so the sample is placed at `(0, 0)`
- **Missing legend placement**: the placements of the legends are computed from the panels of the dataset when it loads. When a dataset without a map gets coordinates from a dropped file, the map has no legend placement, and its legend fails with "Cannot read properties of undefined (reading 'horizontal')"

The Bun patch of Auspice, `patches/auspice@3.1.0.patch`, fixes these: coordinates alone make a valid file, samples are grouped by their parsed coordinates, rows with an empty coordinate and rows of names that are not in the tree are skipped, and the metadata update adds the legend placement of a map that it makes available.

## Code

- Column check: [actions/filesDropped/parseCsv.ts#L81](https://github.com/nextstrain/auspice/blob/v3.1.0/src/actions/filesDropped/parseCsv.ts#L81)
- Coordinate grouping: [actions/filesDropped/parseCsv.ts#L313-L318](https://github.com/nextstrain/auspice/blob/v3.1.0/src/actions/filesDropped/parseCsv.ts#L313-L318)
- Legend placements: [util/metadataJsonParsing.ts#L113](https://github.com/nextstrain/auspice/blob/v3.1.0/src/util/metadataJsonParsing.ts#L113)
- Geographic metadata update: [actions/updateMetadata/updateMetadata.ts#L182-L189](https://github.com/nextstrain/auspice/blob/v3.1.0/src/actions/updateMetadata/updateMetadata.ts#L182-L189)

## Fix

Contribute the changes to [nextstrain/auspice](https://github.com/nextstrain/auspice), then remove them from the patch after upgrading to an Auspice release that contains them. [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md) describes how to edit the patch.
