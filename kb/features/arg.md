# ARG Inference

- [/] **Partial in v1**: `ancestral` and `optimize` with `--segment-trees --segment-alignments` reconstruct segments and fit shared branch lengths on one graph of named segment trees, and `timetree` dates the graph at a fixed clock rate ([kb/issues/M-timetree-segment-graphs-fixed-clock-only.md](../issues/M-timetree-segment-graphs-fixed-clock-only.md))

## v0 Features (All Missing)

- [ ] Two-tree input (`--trees`, `--alignments`)
- [ ] MCC (most compatible clades) file input
- [ ] Per-tree timetree inference
- [ ] Combined ARG output

## v1 direction

One graph for all segments, with edge scopes that record which segments (or sites) each edge carries: [kb/proposals/edge-scoped-inheritance.md](../proposals/edge-scoped-inheritance.md).
