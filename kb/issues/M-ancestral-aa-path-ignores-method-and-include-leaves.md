# Amino-acid reconstruction ignores `--method-anc` and `--include-leaves`

The per-CDS amino-acid path of `treetime ancestral --translations` (`fn reconstruct_aa()` in [packages/treetime/src/ancestral/aa.rs](../../packages/treetime/src/ancestral/aa.rs), which calls the shared `fn reconstruct_partition()` in [packages/treetime/src/ancestral/plan.rs](../../packages/treetime/src/ancestral/plan.rs) once per CDS) does not follow two flags that the nucleotide path follows. None of them produces a warning.

- **`--method-anc`**: the AA path always builds a `ReconstructionPlan::Marginal` ([packages/treetime/src/ancestral/aa.rs#L64-L67](../../packages/treetime/src/ancestral/aa.rs#L64-L67)). With `--method-anc parsimony`, nucleotides use Fitch while amino acids use marginal reconstruction
- **`--include-leaves`**: `AaParams.include_leaves` reaches `stream_sequences()`, which sets the `emitted` field of each `SeqItem` ([packages/treetime/src/seq/sink.rs#L15-L20](../../packages/treetime/src/seq/sink.rs#L15-L20)). The AA FASTA sink `AaFastaSink` ignores this field and writes every node it receives, tips included ([packages/app-commands/src/commands/ancestral/run.rs#L413-L428](../../packages/app-commands/src/commands/ancestral/run.rs#L413-L428)). The nucleotide sink writes only items with `emitted` set ([packages/app-commands/src/commands/ancestral/run.rs#L242](../../packages/app-commands/src/commands/ancestral/run.rs#L242))

## Compatibility target

Augur, the target for the AA path, states that amino-acid sequences "are inferred with the same method as the nucleotide sequences" (`augur/ancestral.py`, module docstring). Augur writes all nodes, tips included, in both its nucleotide and amino-acid FASTA outputs. The Fitch arm needed for parsimony already works on the AA alphabet: `ReconstructionPlan::Fitch` in [packages/treetime/src/ancestral/plan.rs](../../packages/treetime/src/ancestral/plan.rs) takes any alphabet.

## Open question

> [!IMPORTANT]
> **Decision required.** The AA path silently diverges from the nucleotide path for two flags. For each flag the options are: honor it on the AA path, reject the combination with an error, or document it as nucleotide-only. Augur infers amino acids with the nucleotide method (see above), which supports honoring `--method-anc`.
