# The web and desktop apps show no Auspice panel besides the tree

The results page embeds only the Auspice tree and its sidebar controls. Nextstrain apps that embed Auspice also show the diversity (entropy) panel below the tree, and nextstrain.org shows the map and frequencies panels when a dataset provides their data.

## Evidence

- `fn AuspiceTree()` in `packages/app-ui/src/auspice/AuspiceTree.tsx` renders the controls and `auspice/src/components/tree`, and no other Auspice panel
- The Auspice JSON that TreeTime writes declares `meta.panels: ["tree"]` and has no `meta.genome_annotations` and no `meta.geo_resolutions` (checked on a `timetree` run of `data/mpox/clade-ii/1000`)
- The branches of that JSON carry nucleotide mutations (`branch_attrs.mutations.nuc`), which is the input of the Auspice entropy panel
- Nextclade embeds `auspice/src/components/entropy` below the tree and imports Auspice's `entropy.css`, `select.css` and `notifications.css` (`packages/nextclade-web/src/components/Tree/Tree.tsx` and `packages/nextclade-web/src/styles/auspice.scss` in the Nextclade repository)

## Impact

- Users do not see which alignment positions vary along the tree, which Auspice can show from data TreeTime already computes

## Open question

Which Auspice panels should the apps show, and what must the CLI write for them?

- Entropy panel: needs `meta.genome_annotations.nuc` (start and end of the alignment) and `"entropy"` in `meta.panels` in the Auspice JSON the CLI writes, which changes the CLI output, and the panel in `AuspiceTree`
- Map: needs geographic coordinates for a metadata column, which TreeTime does not read today
- Frequencies: needs a tip-frequencies file, which TreeTime does not compute
