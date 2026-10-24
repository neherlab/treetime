# The web and desktop apps show no Auspice frequencies panel

The results page embeds the Auspice tree, its sidebar controls, the map of a `mugration` run with coordinates (`mugration --lat-longs`), and the diversity (entropy) panel. nextstrain.org also shows the frequencies panel when a dataset provides its data.

## Evidence

- `function AuspicePanels()` in `packages/app-ui/src/auspice/AuspicePanels.tsx` renders the tree, the map and the entropy chart, and no other Auspice panel
- The Auspice JSON that TreeTime writes declares at most `meta.panels: ["tree", "map", "entropy"]` (`fn auspice_data()` in `packages/app-output/src/auspice.rs`)
- Auspice shows the frequencies panel only when a tip-frequencies sidecar file is available (`auspice/src/actions/recomputeReduxState.js`)

## Impact

- Users cannot see how the frequencies of clades or trait values change over time

## Open question

> [!IMPORTANT]
> **Decision required.** Should the apps show the frequencies panel, and what must the CLI write for it? The panel needs a tip-frequencies file, which TreeTime does not compute, so the option adds a new computation to the CLI, because the apps show only what the CLI writes.
