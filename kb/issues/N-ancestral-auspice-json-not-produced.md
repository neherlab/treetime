# Ancestral command does not produce auspice JSON

v0's `ancestral` command produces `auspice_tree.json` via `export_sequences_and_tree()` (`CLI_io.py:220`). v1's ancestral command produces Newick, Nexus, Graphviz, FASTA, and GTR JSON but no auspice JSON.

v0's ancestral auspice JSON contains `node_attrs.div` (cumulative mutation_length), `branch_attrs.mutations.nuc` (per-branch mutations), `node_attrs.confidence` (pseudo-bootstrap), `meta.genome_annotations.nuc`, and `node_attrs.bad_branch`. No dates (timetree=False).

This is a standalone visualization convenience for viewing ancestral reconstruction results in auspice without going through `augur export v2`. The augur pipeline path is covered by node data JSON output.

## Locations

- v0: `packages/legacy/treetime/treetime/CLI_io.py:277-341` `create_auspice_json()`, called at line 220
- v0: `packages/legacy/treetime/treetime/wrappers.py:641-648` ancestral command output
- v1: `packages/treetime/src/commands/ancestral/run.rs` - no auspice output
- v1 timetree auspice: `packages/treetime/src/commands/timetree/output/auspice.rs` - exists but ancestral has no equivalent

## Related

- [N-timetree-auspice-json-incomplete.md](N-timetree-auspice-json-incomplete.md) - same three features missing in timetree's auspice output
- [N-io-write-graph-files-missing-formats.md](N-io-write-graph-files-missing-formats.md) - auspice as missing generic output format
- [../reports/augur-node-data-json.md](../reports/augur-node-data-json.md) - node data JSON shares the same data sources (partition mutations, annotations)
