# Add timetree divergence-tree output

Timetree should output two trees: one with time-based branch lengths and one with divergence-based branch lengths (substitutions/site). v0 does this via `BranchLengthMode` in `CLI_io.py:210-233`.

## Scope

- Add `divergence-tree` tree identity to timetree command alongside existing `timetree` identity
- `timetree` output: branch lengths in time units (date difference between parent and child)
- `divergence-tree` output: branch lengths in substitutions/site (current behavior, from `edge.branch_length()`)
- Both tree identities get the full set of selected output formats from ticket #4
- Per-file flags use tree identity prefix: `--output-timetree-nwk`, `--output-divergence-tree-nwk`

## v0 reference

`packages/legacy/treetime/treetime/CLI_io.py:210-233`: `export_sequences_and_tree` writes `timetree.nexus` with time branch lengths, then swaps `n.branch_length = n.mutation_length` and writes `divergence_tree.nexus`.

## Tests

### Happy paths

- `timetree.nwk` branch lengths are time differences between parent and child dates (years)
- `divergence_tree.nwk` branch lengths are substitutions/site (current behavior)
- Both trees have identical topology (same node count, same parent-child relationships, same leaf names)
- Both trees produced in all selected formats (nwk, nexus, etc.)
- Per-file flags: `--output-timetree-nwk=time.nwk`, `--output-divergence-tree-nwk=div.nwk`
- `--output-selection=nwk` produces both `timetree.nwk` and `divergence_tree.nwk`

### Edge cases

- Zero time interval: parent and child have same inferred date -> time branch length is 0.0
- Root node: time branch length is undefined or 0.0 (no parent). Document convention
- Polytomy: multiple children at same date -> all time branch lengths to parent are equal
- Node with very small time difference (1e-10 years) -> precision preserved

### Golden master

- Compare `timetree.nwk` branch lengths against v0's `timetree.nexus` for each test dataset (`data/flu/h3n2/20`, `data/ebola`, etc.)
- Compare `divergence_tree.nwk` branch lengths against v0's `divergence_tree.nexus`
- Tolerance: 1e-6 (v0 comparison standard)

### Regression

- Existing timetree output content unchanged (divergence-tree is additive)
- Existing timetree tests still pass

## Related issues

- Source: [kb/issues/N-io-time-based-branch-lengths-not-implemented.md](../issues/N-io-time-based-branch-lengths-not-implemented.md)
