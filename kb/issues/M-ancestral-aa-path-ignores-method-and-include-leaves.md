# Amino-acid reconstruction ignores `--method-anc` and `--include-leaves`

The per-CDS amino-acid path of `treetime ancestral --translations` (`fn reconstruct_aa()` in [packages/treetime/src/ancestral/aa.rs](../../packages/treetime/src/ancestral/aa.rs), which calls the shared `fn reconstruct_partition()` in [packages/treetime/src/ancestral/plan.rs](../../packages/treetime/src/ancestral/plan.rs) once per CDS) does not follow two flags that the nucleotide path follows. None of them produces a warning.

- **`--method-anc`**: the AA path always passes `ReconstructionMethod::Marginal(representation)` to `fn reconstruct_partition()` ([packages/treetime/src/ancestral/aa.rs#L66](../../packages/treetime/src/ancestral/aa.rs#L66)). With `--method-anc parsimony`, nucleotides use Fitch while amino acids use marginal reconstruction
- **`--include-leaves`**: `AaParams.include_leaves` reaches `stream_sequences()`, which sets the `emitted` field of each `SeqItem` ([packages/treetime/src/seq/sink.rs#L20-L27](../../packages/treetime/src/seq/sink.rs#L20-L27)). The AA FASTA sink `AaFastaSink` ignores this field and writes every node it receives, tips included ([packages/app-commands/src/commands/ancestral/run.rs#L456-L475](../../packages/app-commands/src/commands/ancestral/run.rs#L456-L475)). The nucleotide sink writes only items with `emitted` set ([packages/app-commands/src/commands/ancestral/run.rs#L274](../../packages/app-commands/src/commands/ancestral/run.rs#L274))

## Compatibility target

Augur, the target for the AA path, states that amino-acid sequences "are inferred with the same method as the nucleotide sequences" (`augur/ancestral.py`, module docstring). Augur writes all nodes, tips included, in both its nucleotide and amino-acid FASTA outputs. The Fitch arm needed for parsimony already works on the AA alphabet: `ReconstructionMethod::Fitch` in [packages/treetime/src/ancestral/plan.rs#L259](../../packages/treetime/src/ancestral/plan.rs#L259) takes any alphabet.

## Open question

> [!IMPORTANT]
> **Decision required.** The AA path silently diverges from the nucleotide path for two flags. For each flag the options are: honor it on the AA path, reject the combination with an error, or document it as nucleotide-only. Augur infers amino acids with the nucleotide method (see above), which supports honoring `--method-anc`.

## Related

- [kb/proposals/codon-aware-ancestral-reconstruction.md](../proposals/codon-aware-ancestral-reconstruction.md): axis 1 recommends that TreeTime translates its own reconstruction, which limits this path to protein-only input

## Proposed implementation plan

One short plan per option of the Decision required block; the steps apply per flag.

- **Honor the flag**:
  - `--method-anc`: add the method to `AaParams` and map it to `ReconstructionMethod` with the same function that the nucleotide path uses ([packages/treetime/src/ancestral/plan.rs#L41](../../packages/treetime/src/ancestral/plan.rs#L41)), in place of the fixed `ReconstructionMethod::Marginal` in `fn reconstruct_aa()`. With Fitch, the fitted GTR model and rate of each CDS are unused, as on the nucleotide path
  - `--include-leaves`: in `impl SeqSink for AaFastaSink` ([packages/app-commands/src/commands/ancestral/run.rs#L456](../../packages/app-commands/src/commands/ancestral/run.rs#L456)) skip items whose `emitted` is false, as the nucleotide sink does at line 274
- **Reject the combination**: return an error from argument validation of `treetime ancestral` when `--translations` is given with `--method-anc parsimony` or without `--include-leaves`, with a message that names the flag
- **Nucleotide-only**: state the limit in the `--help` doc comment of each flag and in the docs of `treetime ancestral`

Tests (oracle: augur, which infers amino acids with the nucleotide method and writes tips in both FASTA outputs; `augur/ancestral.py`):

- `packages/treetime/src/ancestral/__tests__/test_aa.rs`: with the parsimony method on a small tree with one CDS, the amino acids of internal nodes equal the Fitch result of `ReconstructionMethod::Fitch` on the same peptides
- A sink test in `packages/app-commands`: an `AaFastaSink` given one item with `emitted` false and one with `emitted` true writes only the second
- For the reject option: one `assert_err!` case per flag

Validation: `just check`, `just test-rs`, and `just r treetime ancestral` with `--translations` and `--annotation` on `rsv/a/20` with and without `--include-leaves` and with `--method-anc parsimony`.
