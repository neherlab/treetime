# Simulation substitution models support only nucleotides

`treetime-simulation` builds substitution models over the 4 nucleotides only ([packages/treetime-simulation/src/substitution_model.rs#L28](../../packages/treetime-simulation/src/substitution_model.rs#L28)). It cannot simulate codon evolution, so no test can measure how well a codon model or codon decoding recovers known ancestral states.

The one published benchmark of codon against nucleotide ancestral states is a tie (98.4% against 98.3% per-site accuracy, section "Accuracy of ancestral states" of [kb/reports/codon-aware-ancestral-reconstruction.md](../reports/codon-aware-ancestral-reconstruction.md#accuracy-of-ancestral-states)), so a codon-model default needs a benchmark at viral divergence.

## Fix direction

Add an n-state substitution model, and codon states mapped to triplets in the simulated sequences. Optional: rate categories among codons and multi-nucleotide events. The benchmark ([kb/proposals/codon-aware-ancestral-reconstruction.md](../proposals/codon-aware-ancestral-reconstruction.md), axis 15) compares codon decoding of the nucleotide reconstruction with a codon model on simulated trees at the divergence of `rsv/a` and `sc2`: amino-acid and nucleotide accuracy at internal nodes, stop codons, and calibration of posteriors.

## Validation

Simulated codon frequencies at stationarity match the model's stationary distribution, and simulated substitution counts per branch match the expected counts.
