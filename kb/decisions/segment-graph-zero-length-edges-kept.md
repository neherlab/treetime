# Zero-length edges of segment graphs are kept

A segment graph keeps every zero-length edge of its input trees and does not contract it. On the H3N2 graph made from TreeKnit output, 391 of the 801 nodes with several parents have a zero-length parent edge ([kb/issues/M-timetree-segment-graphs-fixed-clock-only.md](../issues/M-timetree-segment-graphs-fixed-clock-only.md)).

**Type**: Scientific default of a new v1 input. Approved.

## Why

TreeKnit expects collapsed input trees, in which every branch is supported, and resolves their polytomies with splits from compatible other segments. A zero-length edge in its output is therefore a split that another segment supports, not an arbitrary resolution of a polytomy. Contracting it would remove a node that another segment shares and would give that segment's split to no node of this segment.

## Consequences

- **Inputs to TreeKnit should be collapsed trees**: in a tree that is not collapsed, a zero-length edge can be an arbitrary resolution, which the graph would then keep as shared structure
- **No short-branch collapse on segment graphs**: `timetree --segment-trees` keeps the polytomies and the short branches of its input ([kb/issues/M-timetree-segment-graphs-fixed-clock-only.md](../issues/M-timetree-segment-graphs-fixed-clock-only.md))
