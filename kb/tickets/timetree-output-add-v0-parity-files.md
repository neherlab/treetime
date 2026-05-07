# Add missing output files for v0 parity

v1 produces fewer output files than v0. Several output types are not implemented.

## v0 outputs not present in v1

- `dates.tsv` - node date estimates in tabular format
  (see [Node dates output unimplemented](timetree-output-node-dates-unimplemented.md))
- `ancestral_sequences.fasta` - reconstructed ancestral sequences
- `branch_mutations.txt` - mutations mapped to branches in tabular form (Nexus `[&mutations=...]` annotations are emitted; only the separate tabular file is missing)
- `molecular_clock.txt` - clock model summary (v1 has `timetree.json` with more detail)
- `divergence_tree.nexus` - tree with divergence (not time) branch lengths
- `sequence_evolution_model.txt` - GTR model details (v1 has `gtr.json`)
- `root_to_tip_regression.pdf` - RTT regression plot
  (see [Plot commands unimplemented](timetree-output-plot-rtt-and-plot-tree-unimplemented.md))
- `timetree.pdf` - time tree visualization
  (see [Plot commands unimplemented](timetree-output-plot-rtt-and-plot-tree-unimplemented.md))
- `skyline.tsv` - effective population size over time
  (see [Missing skyline output files](timetree-output-missing-skyline-files.md))
- `skyline.pdf` - skyline plot

## Produced but incomplete

- `auspice_tree.json` - v1 produces this file, but it is missing mutations, branch confidence, and genome annotations compared to v0 (see [Auspice JSON output missing mutations, branch confidence, and genome annotations](../issues/N-timetree-auspice-json-incomplete.md))

## v1-only outputs

- `timetree.json` - clock model in JSON format
- `timetree.nwk` - tree in Newick format (v0 only writes Nexus)
- `gtr.json` - GTR model in JSON format

## Related issues

- Source: [N-timetree-missing-output-files.md](../issues/N-timetree-missing-output-files.md) -- delete after full resolution
