# Missing output files compared to v0

v1 produces fewer output files than v0. Several output types are not implemented.

## v0 outputs not present in v1

- `dates.tsv` - node date estimates in tabular format
  (see [Node dates output unimplemented](N-timetree-node-dates-output-unimplemented.md))
- `ancestral_sequences.fasta` - reconstructed ancestral sequences
- `branch_mutations.txt` - mutations mapped to branches in tabular form (Nexus `[&mutations=...]` annotations are emitted; only the separate tabular file is missing)
- `molecular_clock.txt` - clock model summary (v1 has `timetree.json` with more detail)
- `divergence_tree.nexus` - tree with divergence (not time) branch lengths
- `sequence_evolution_model.txt` - GTR model details (v1 has `gtr.json`)
- `root_to_tip_regression.pdf` - RTT regression plot
  (see [Plot commands unimplemented](N-timetree-plot-unimplemented.md))
- `timetree.pdf` - time tree visualization
  (see [Plot commands unimplemented](N-timetree-plot-unimplemented.md))
- `skyline.tsv` - effective population size over time
  (see [Missing skyline output files](N-timetree-missing-skyline-output.md))
- `skyline.pdf` - skyline plot

## v1-only outputs

- `timetree.json` - clock model in JSON format
- `timetree.nwk` - tree in Newick format (v0 only writes Nexus)
- `gtr.json` - GTR model in JSON format
- `timetree.augur-node-data.json` - augur-compatible node data JSON (treetime equivalent of `augur refine` output): per-node dates, branch lengths, divergence metrics, confidence intervals, and the clock model. Carries the node date estimates that the unimplemented `dates.tsv` would provide (see [Node dates output unimplemented](N-timetree-node-dates-output-unimplemented.md))

## Potential solutions

- O1. Implement a v0 output when its distinct file contract remains required.
- O2. Treat an existing structured v1 output as the approved replacement after verifying that it preserves the same information and workflow interoperability.
- O3. Omit an obsolete visualization or duplicate representation through an explicit parity decision.

## Recommendation

Decide each output independently. Use the existing focused issues for node dates, plots, skyline files, and TreeTime metadata; create separate issues for remaining files whose semantics are still required. Do not use one ticket for unrelated tabular, model, tree, and visualization artifacts.

## Ticket readiness

This inventory has no executable aggregate ticket. Focused source issues determine readiness for each output.
