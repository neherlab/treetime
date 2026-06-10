# Newick writer places all annotations after branch length

`nwk_write_with` at `packages/treetime-io/src/nwk.rs#L242-L258` emits: `name:length[&annotations]`.

BEAST/FigTree convention distinguishes node vs branch annotations by position:

- Node annotations before colon: `name[&node_stuff]:length`
- Branch annotations after colon: `name:length[&branch_stuff]`

Current code conflates both into one post-length position, losing the node-vs-branch distinction. Downstream tools (FigTree, BEAST) interpret post-length annotations as branch properties, so node-level annotations (confidence, dates, states) would be misclassified.

## Impact

Node-level annotations are written in the branch-annotation position. Tools that distinguish positions will misinterpret them.

## Locations

- `packages/treetime-io/src/nwk.rs#L242-L258`

## Related issues

- [kb/tickets/io-nwk-writer-style-dispatch.md](../tickets/io-nwk-writer-style-dispatch.md)
