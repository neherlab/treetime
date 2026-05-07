# Auspice JSON output missing mutations, branch confidence, and genome annotations

The timetree `auspice_tree.json` output omits three features present in v0.

## Missing features

### Branch mutations (`branch_attrs.mutations`)

v0 populates `branch_attrs.mutations.nuc` with nucleotide substitutions for each branch (e.g. `["A100G", "C200T"]`). v1 emits an empty mutations object. The partition-layer mutation data exists but is not accessible from the graph payload layer. Blocked by the unified mutation API gap (`M-core-branch-mutations-no-unified-api`).

### Branch support confidence (`node_attrs.confidence`)

v0 computes a pseudo-bootstrap confidence value per internal node: `1 - exp(-n_mutations)` where `n_mutations` is the count of ACGT substitutions on the branch. Leaves get confidence 1.0. This requires branch mutation counts, so it is blocked by the same mutation API gap. The corresponding `meta.colorings` entry `{'title': 'Branch Support', 'type': 'continuous', 'key': 'confidence'}` is also missing.

### Genome annotations (`meta.genome_annotations`)

v0 writes `genome_annotations.nuc` with `start`, `end`, `type`, and `strand` derived from `tt.data.full_length` (the alignment length). v1's timetree graph data type is `()` and does not carry sequence length. The alignment length is available in the partition layer but not passed to the output writer.

## Current state

`packages/treetime/src/commands/timetree/output/auspice.rs` writes `auspice_tree.json` with `num_date` (including confidence intervals), `div`, and `bad_branch`. The `panels` metadata includes only `["tree"]`; v0 also includes `"entropy"` which depends on genome annotations.

## Impact

Auspice visualization works for time tree display and coloring. The entropy panel and mutation tooltips are unavailable. Branch confidence bars (distinct from date confidence) are absent.

## Related issues

- Source: [N-timetree-auspice-json-incomplete.md](../issues/N-timetree-auspice-json-incomplete.md) -- delete after full resolution
