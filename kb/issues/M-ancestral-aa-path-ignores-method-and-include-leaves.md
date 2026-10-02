# Amino-acid reconstruction ignores `--method-anc`, `--gtr-iterations` and `--include-leaves`

The per-CDS amino-acid path of `treetime ancestral --translations` (`fn reconstruct_aa` in [packages/treetime/src/ancestral/aa.rs](../../packages/treetime/src/ancestral/aa.rs), which calls the shared `fn reconstruct_partition` in [packages/treetime/src/ancestral/plan.rs](../../packages/treetime/src/ancestral/plan.rs) once per CDS) does not follow three flags that the nucleotide path follows. None of them produces a warning.

- **`--method-anc`**: the AA path always builds a `ReconPlan::Marginal`. With `--method-anc parsimony`, nucleotides use Fitch while amino acids use marginal reconstruction
- **`--gtr-iterations`**: `AaParams` has no iteration field and each CDS plan sets `gtr_refinement: None`, so the AA path never refines the GTR
- **`--include-leaves`**: the AA FASTA sink writes every graph node, tips included, and ignores the emitted-node list that the nucleotide FASTA uses

## Compatibility target

Augur, the target for the AA path, states that amino-acid sequences "are inferred with the same method as the nucleotide sequences" (`augur/ancestral.py`, module docstring). Augur writes all nodes, tips included, in both its nucleotide and amino-acid FASTA outputs. The Fitch arm needed for parsimony already works on the AA alphabet: `ReconPlan::Fitch` in [packages/treetime/src/ancestral/plan.rs](../../packages/treetime/src/ancestral/plan.rs) takes any alphabet.

## Open question

For each flag: honor it on the AA path, reject the combination with an error, or document it as nucleotide-only.

## Related issues

- [M-ancestral-gtr-iterations-refit-against-frozen-messages.md](M-ancestral-gtr-iterations-refit-against-frozen-messages.md)
