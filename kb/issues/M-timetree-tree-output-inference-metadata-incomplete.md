# TreeTime tree-output inference metadata is incomplete

TreeTime now produces TreeIR-backed Auspice and PhyloXML output, but its projection remains incomplete as one output contract.

## Problem

- Auspice substitution projection exists only where the command constructs TreeIR; insertion and deletion handling is governed by [M-core-mutation-representation-and-format-projection-inconsistent.md](M-core-mutation-representation-and-format-projection-inconsistent.md).
- `branch_attrs` omits the v0 pseudo-bootstrap branch support confidence.
- `meta.colorings` omits the corresponding continuous `confidence` coloring.
- `meta.genome_annotations.nuc` omits alignment coordinates, type, and strand.
- The entropy panel metadata cannot be complete without genome annotations.
- TreeIR defines `date_is_inferred`, `date_raw_value`, and relaxed-clock `gamma`, but the TreeTime projection leaves them at defaults even though command state contains the source data.
- PhyloXML can serialize inferred-date and `gamma` properties, so command-generated relaxed-clock output remains incomplete.

The Auspice fields share one command projection and must be validated together. The already-defined non-mutation PhyloXML date and gamma properties use the same command state. PhyloXML mutation vocabulary and indel representation remain excluded pending [N-io-phyloxml-mutation-property-contract-undecided.md](N-io-phyloxml-mutation-property-contract-undecided.md).

## Branch support formula

For $m$ canonical nucleotide substitutions on a branch, v0 reports confidence

$$c = 1 - e^{-m}$$

where $c$ is branch support confidence. Leaves receive $c=1$.

## Potential solutions

- O1. Extend the shared TreeIR projection with all mutation and inference metadata needed by Auspice and PhyloXML.
- O2. Build format-specific TreeTime projections. This can preserve format-specific data but duplicates the shared boundary.

## Recommendation

Extend the TreeTime projection once, supplying typed Auspice mutations, confidence, alignment length, annotations, raw/inferred dates, and relaxed-clock rate multipliers. Add an Auspice fixture for the complete visualization payload and a PhyloXML fixture limited to the already-defined non-mutation inference properties.

## Related issues

- [M-core-mutation-representation-and-format-projection-inconsistent.md](M-core-mutation-representation-and-format-projection-inconsistent.md)
- [N-ancestral-auspice-json-not-produced.md](N-ancestral-auspice-json-not-produced.md)
- [N-io-phyloxml-mutation-property-contract-undecided.md](N-io-phyloxml-mutation-property-contract-undecided.md)
- [kb/reports/augur-node-data-json.md](../reports/augur-node-data-json.md)
