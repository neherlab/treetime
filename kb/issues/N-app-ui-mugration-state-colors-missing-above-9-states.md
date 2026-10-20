# Mugration tables show no state colors for traits with more than 9 states

`StateLabel` in `packages/app-ui/src/runs/MugrationResultsView.tsx` draws a colored dot before each state in the "State changes" and "Uncertain ancestors" tables. `MugrationResultsView` [packages/app-ui/src/runs/MugrationResultsView.tsx#L78-L85](../../packages/app-ui/src/runs/MugrationResultsView.tsx#L78-L85) takes the colors from the `scale` of the coloring whose key is the mugration attribute, in the colorings that the back end returns with the tree.

The back end owns these colors. `fn coloring_scale()` [packages/app-backend/src/results/tree.rs#L137-L159](../../packages/app-backend/src/results/tree.rs#L137-L159), called by `fn result_colorings()` [packages/app-backend/src/results/tree.rs#L80](../../packages/app-backend/src/results/tree.rs#L80), pairs the sorted states of a categorical coloring with the 9 colors of `MUTED_PALETTE` when the coloring has 1 to 9 states. `fn display_auspice()` [packages/app-backend/src/results/auspice.rs#L22-L28](../../packages/app-backend/src/results/auspice.rs#L22-L28) writes the same scale into the Auspice JSON that the app displays, so the tree and the tables agree.

With more than 9 states, `coloring_scale()` falls back to the `scale` of the Auspice file. The Auspice file that `treetime mugration` writes has no `scale` for the attribute, so the scale is empty: the tables draw no dot, and Auspice assigns the tree colors in the browser.

## Evidence

Mugration on `data/dengue/1000` with `--attribute country` (56 distinct countries in `metadata.tsv`): `meta.colorings` of `mugration.auspice.json` holds `{"key": "country", "type": "categorical"}` without `scale`, the tree is colored, and neither table shows a dot.

## Fix direction

Give the back end a color for every state of a categorical coloring, and keep the tree and the tables on the same scale.

> [!IMPORTANT]
> **Decision required.** The back end colors at most 9 states. Options for more states: extend the palette (for example a longer qualitative palette, or interpolated colors in the back end), or reproduce the color assignment that Auspice applies to categorical colorings without a scale, so the tree keeps its current colors. Extending the palette changes the tree colors that users see for traits above 9 states; reproducing the Auspice assignment couples the back end to Auspice's coloring rules.
