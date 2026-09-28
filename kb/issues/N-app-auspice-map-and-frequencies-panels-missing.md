# The web and desktop apps show no Auspice map or frequencies panel

The results page embeds the Auspice tree, its sidebar controls, and the diversity (entropy) panel. nextstrain.org also shows the map and frequencies panels when a dataset provides their data.

## Evidence

- `fn AuspiceTree()` in `packages/app-ui/src/auspice/AuspiceTree.tsx` renders the controls, `auspice/src/components/tree`, and `auspice/src/components/entropy`, and no other Auspice panel
- The Auspice JSON that TreeTime writes declares at most `meta.panels: ["tree", "entropy"]` and has no `meta.geo_resolutions` (`fn auspice_data()` in `packages/app-output/src/tree_output.rs`)
- Auspice removes the map panel when `meta.geo_resolutions` is empty, and shows the frequencies panel only when a tip-frequencies sidecar file is available (`auspice/src/actions/recomputeReduxState.js`)

## Impact

- Users of a `mugration` run on a geographic attribute cannot see the inferred locations on a map

## Open question

Which of these panels should the apps show, and what must the CLI write for them?

- Map: needs latitude and longitude for each value of a metadata column in `meta.geo_resolutions`, which TreeTime does not read today
- Frequencies: needs a tip-frequencies file, which TreeTime does not compute
