# Mugration (Discrete Trait Reconstruction)

- [x] **Implemented** (marginal reconstruction complete, ~80%)

## Input Processing

- [x] Tree loading from `--tree`
- [x] Discrete state table from `--states`
- [x] Attribute column selection with `--attribute`
- [x] Custom taxon-name column with `--metadata-id-columns`
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
- [x] Trait assignment from the most likely state of the posterior profiles, with the tie rule of [kb/decisions/posterior-most-likely-state-tie-rule.md](../decisions/posterior-most-likely-state-tie-rule.md)
- [x] Confidence profile and entropy per node (`StatePosterior`)

## Output

- [x] `traits.csv` (per-node trait assignments, all nodes)
- [x] `annotated_tree.nexus` (`mugration.nexus`: Nexus with the trait in BEAST comments)
- [x] `annotated_tree.nwk` (`mugration.nwk`: Newick with the trait in comments with `--newick-annotations beast` or `nhx`, plain otherwise)
- [x] `gtr.json` (model document with one `trait` entry: attribute, ordered state labels, model parameters)

## Additional Features

- [x] Iterative GTR inference in the v0 order (golden masters: traits and confidences match v0 on nine cases)
- [x] Sampling bias correction (`--sampling-bias-correction`)
- [x] Confidence CSV output (`--confidence`)
- [x] `--pc` pseudo-counts
- [x] Optional ambiguity report in the traits CSV, augur node data, and Auspice JSON (`--ambiguity-threshold`, `--ambiguity-outputs`; new in v1, [kb/decisions/ambiguity-report-of-plausible-states.md](../decisions/ambiguity-report-of-plausible-states.md))
