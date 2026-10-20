# The web and desktop apps show no Auspice map or frequencies panel

The results page embeds the Auspice tree, its sidebar controls, and the diversity (entropy) panel. nextstrain.org also shows the map and frequencies panels when a dataset provides their data.

## Evidence

- `function AuspiceTree()` in `packages/app-ui/src/auspice/AuspiceTree.tsx` renders the controls, `auspice/src/components/tree`, and `auspice/src/components/entropy`, and no other Auspice panel [packages/app-ui/src/auspice/AuspiceTree.tsx#L121-L122](../../packages/app-ui/src/auspice/AuspiceTree.tsx#L121-L122)
- The Auspice JSON that TreeTime writes declares at most `meta.panels: ["tree", "entropy"]` and has no `meta.geo_resolutions` (`fn auspice_data()` [packages/app-output/src/auspice.rs#L102-L105](../../packages/app-output/src/auspice.rs#L102-L105))
- Auspice removes the map panel when `meta.geo_resolutions` is empty, and shows the frequencies panel only when a tip-frequencies sidecar file is available (`auspice/src/actions/recomputeReduxState.js`)

## Impact

- Users of a `mugration` run on a geographic attribute cannot see the inferred locations on a map

## Open question

> [!IMPORTANT]
> **Decision required.** Which of these panels should the apps show, and what must the CLI write for them? Options, one or both:
>
> - Map: needs latitude and longitude for each value of a metadata column in `meta.geo_resolutions`, which TreeTime does not read today
> - Frequencies: needs a tip-frequencies file, which TreeTime does not compute
>
> Each option adds a new input or a new computation to the CLI, because the apps show only what the CLI writes.
