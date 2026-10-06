# Homoplasy Analysis

- [x] **Homoplasy command** (`treetime homoplasy`), with the divergences of [kb/decisions/homoplasy-mutation-mapping-and-counting.md](../decisions/homoplasy-mutation-mapping-and-counting.md)
- [x] Tree input from `--tree` (required; v0 tree building is not ported, [kb/issues/H-timetree-tree-inference-unimplemented.md](../issues/H-timetree-tree-inference-unimplemented.md))
- [x] Alignment input from `--alignment`
- [ ] VCF input with `--vcf-reference` ([kb/issues/M-io-vcf-input-output-unimplemented.md](../issues/M-io-vcf-input-output-unimplemented.md))

## Mutation Mapping

- [x] Ancestral reconstruction with `--method-anc` (marginal default, parsimony); v0 uses joint reconstruction
- [x] Substitution model flags of `ancestral` (`--model`, `--gtr-iterations`, `--dense`, `--gap-fill`)
- [x] Unknown-state bridging of substitutions

## v0 Features

- [x] Mutation multiplicity distribution
- [x] Position hit count distribution
- [x] Poisson comparison (expected vs observed, log-likelihood difference)
- [x] Top-N homoplasic mutations display
- [x] `--detailed` (terminal branch mutations, strains with homoplasies)
- [x] `--drms` (drug resistance mutation annotation)
- [x] `--const` (constant sites correction)
- [x] `--rescale` (branch length rescaling)
- [x] `-n` (number of mutations/nodes to display)
- [x] `--zero-based` (positions counted from 0)

## v1 Additions

- [x] Changes involving ambiguous characters in their own lists
- [x] Insertions and deletions in their own lists
- [x] Statistics JSON (`--output-homoplasy-stats`) and text report (`--output-homoplasy-report`)
- [x] Tree outputs with branch mutations (Newick, Nexus, Auspice, MAT, graph JSON, Graphviz)
- [x] Web and desktop apps: results page with the recurrent mutations linked to the tree, the sites hit more than once along the genome, and the observed site counts against the Poisson expectation
