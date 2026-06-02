# Ancestral AA: empirical model maps out-of-alphabet characters to unknown

When an empirical amino-acid model (e.g. JTT92) is selected via `--aa-model`, input translations may contain characters outside the model's 20-state alphabet (`AaNoStop`). The most common case is the stop codon `*`, which appears in terminal positions of CDS translations.

## Behavior

Out-of-alphabet characters (notably `*`) are mapped to the unknown state `X` before reconstruction, and a warning is emitted with the count of mapped characters per CDS. The unknown state is marginalized over all 20 canonical amino acids during the probabilistic inference, producing the posterior-most-likely residue at each such position.

## Alternatives considered

Reject-on-stop: error if any input sequence contains `*`. Rejected because stop codons in terminal positions are routine in CDS translations and would force users to pre-strip stops before running the pipeline. augur does not reject stops.

## Default model

The default AA model is `infer`, which operates on the 22-state `aa` alphabet (20 canonical + gap + stop). Stop codons are modeled as a real state under `infer`, so the mapping to unknown only applies when an empirical model is selected explicitly.

## Implementation

- `commands/ancestral/aa_model.rs`: `AaModelName::Jtt92` resolves to `AlphabetName::AaNoStop` (20-state)
- `ancestral/attach.rs`: `sanitize_to_alphabet()` maps characters absent from the target alphabet to its unknown state
- `commands/ancestral/run.rs`: `run_aa_reconstructions()` reads translations with the stop-inclusive `Aa` alphabet, then sanitizes to the reconstruction alphabet, logging the count
