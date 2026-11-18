# Simulation substitution models support only nucleotides

`treetime-simulation` builds substitution models over the 4 nucleotides only ([packages/treetime-simulation/src/substitution_model.rs#L28](../../packages/treetime-simulation/src/substitution_model.rs#L28)). It cannot simulate codon evolution, so no test can measure how well a codon model or codon decoding recovers known ancestral states.

The one published benchmark of codon against nucleotide ancestral states is a tie (98.4% against 98.3% per-site accuracy, section "Accuracy of ancestral states" of [kb/reports/codon-aware-ancestral-reconstruction.md](../reports/codon-aware-ancestral-reconstruction.md#accuracy-of-ancestral-states)), so a codon-model default needs a benchmark at viral divergence.

## Fix direction

Add an n-state substitution model, and codon states mapped to triplets in the simulated sequences. Optional: rate categories among codons and multi-nucleotide events. The benchmark ([kb/proposals/codon-aware-ancestral-reconstruction.md](../proposals/codon-aware-ancestral-reconstruction.md), axis 15) compares codon decoding of the nucleotide reconstruction with a codon model on simulated trees at the divergence of `rsv/a` and `sc2`: amino-acid and nucleotide accuracy at internal nodes, stop codons, and calibration of posteriors.

## Validation

Simulated codon frequencies at stationarity match the model's stationary distribution, and simulated substitution counts per branch match the expected counts.

## Proposed implementation plan

1. Generalize `struct SubstitutionModel` in `packages/treetime-simulation/src/substitution_model.rs` from the fixed `NUCLEOTIDES` alphabet to an alphabet given at construction (a `Vec<AsciiChar>` or a states type), with `pi` of length n, `w` of shape n×n, and the zero-diagonal check against `Array1::zeros(n)` instead of the literal `[0.0; 4]`.
2. Add codon states: a 61-sense-codon alphabet (standard genetic code) and an expansion step in `evolve.rs` or `observe.rs` that writes each codon state as its triplet into the simulated nucleotide sequence. Keep the nucleotide path as the n = 4 case.
3. Optional, behind separate parameters: gamma rate categories among codon sites, and multi-nucleotide events (nonzero exchangeabilities between codons that differ at 2 or 3 positions).
4. Update callers (`observe.rs`, `evolve.rs`) and the existing tests in `packages/treetime-simulation/src/__tests__/` (`test_substitution_model.rs`, `test_prop_substitution_model_normalized.rs`, `test_stat_evolve.rs`) to the new constructor.

Tests:

- Property test (`test_prop_substitution_model_normalized.rs`): for random n and random valid `pi`, `w`, the rate matrix rows sum to zero and the expected rate under `pi` is `mu`. Oracle: the GTR definition.
- Statistical test (`test_stat_evolve.rs`, slow tier): codon frequencies at long branches match `pi`, and substitution counts per branch match `mu * t` within a stated false-fail rate. Oracle: stationarity and the Poisson count of a CTMC.
- Unit test of the triplet expansion against the standard genetic code table.

Validation: `just test-rs -E 'package(treetime-simulation)'`, `just lint-rs`, and the slow simulation tests on user request. The codon benchmark itself (proposal axis 15) is a separate step that uses this model.
