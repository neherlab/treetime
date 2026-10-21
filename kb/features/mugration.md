# Mugration (Discrete Trait Reconstruction)

- [x] **Implemented** (marginal reconstruction complete, ~80%)

## Input Processing

- [x] Tree loading from `--tree`
- [x] Discrete state table from `--states`
- [x] Attribute column selection with `--attribute`
- [x] Custom taxon-name column with `--name-column`
- [x] Auto-detect taxon-name column from `name`, `strain`, or `accession`
- [x] Weights file parsing and validation
- [x] Merge categories observed in states and weights files
- [x] Warn on categories missing from weights file
- [x] Enforce `--missing-weights-threshold`
- [x] Fill missing weights with mean observed weight
- [x] Normalize weights to sum to one
- [x] Build discrete state alphabet of names (sorted, no state limit) excluding missing-data token
- [x] Reject datasets with fewer than two non-missing states
- [x] Treat tree leaves without a metadata row as missing data, with a warning per leaf and a count warning

## Reconstruction

- [x] GTR model construction (uniform rates, optional weight-based equilibrium frequencies)
- [x] Discrete marginal reconstruction (backward pass, forward pass)
- [x] Missing data handling (all-ones profile for tips with the missing-data value or without a metadata row)
- [x] Trait assignment from argmax of posterior profiles
- [x] Confidence profile and entropy per node (`StatePosterior`)

## Output

- [x] `traits.csv` (per-node trait assignments, all nodes)
- [x] `annotated_tree.nexus` (Newick with trait annotations)
- [x] `annotated_tree.nwk` (Newick with NHX-style annotations)
- [x] `gtr.json` (model document with one `trait` entry: attribute, ordered state labels, model parameters)

## Additional Features

- [x] Iterative GTR inference in the v0 order (golden master parity tracked in [kb/issues/H-mugration-iterative-gtr.md](../issues/H-mugration-iterative-gtr.md))
- [x] Sampling bias correction (`--sampling-bias-correction`)
- [x] Confidence CSV output (`--confidence`)
- [x] `--pc` pseudo-counts
