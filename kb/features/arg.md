# ARG Inference

- [/] **Partial in v1**: `ancestral --segment-trees --segment-alignments` reconstructs segments by parsimony on one graph of named segment trees

## v0 Features (All Missing)

- [ ] Two-tree input (`--trees`, `--alignments`)
- [ ] MCC (most compatible clades) file input
- [ ] Per-tree timetree inference
- [ ] Combined ARG output

## v1 direction

One graph for all segments, with edge scopes that record which segments (or sites) each edge carries: [kb/proposals/edge-scoped-inheritance.md](../proposals/edge-scoped-inheritance.md).
