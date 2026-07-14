# Time-based branch lengths not written to Newick/Nexus output

## Summary

The timetree command infers node dates and branch lengths in time units, but the Newick and Nexus output writers only emit divergence-based branch lengths (substitutions per site). v0 writes time-based branch lengths as an alternative output.

## Details

v0 `write_json` and `BranchLengthMode` support emitting branch lengths in calendar time units (`node.branch_length` vs `node.mutation_length`). v1 shared graph output always writes divergence-based lengths from `edge.branch_length()`.

After timetree inference, each node has both a divergence-based branch length and inferred dates. The time-based branch length is the date difference between parent and child. This is available in the graph but not exposed as an output option.

## v0 reference

`packages/legacy/treetime/treetime/treetime.py` `BranchLengthMode` enum, used in `write_json()` and `save_nexus()`.

## Related

- [kb/proposals/output-format-selection.md](../proposals/output-format-selection.md) -- timetree `divergence-tree` identity addresses this

## Related tickets

- [kb/tickets/io-timetree-divergence-tree-output.md](../tickets/io-timetree-divergence-tree-output.md)
