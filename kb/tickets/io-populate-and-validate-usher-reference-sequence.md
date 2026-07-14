# Populate and validate the UShER MAT reference sequence

Mutation-bearing command projections leave `TreeIrData.root_sequence` empty, and `TreeIrUsherWriter` then substitutes each mutation's branch-local parent allele for the global `ref_nuc` field.

## Implementation

- Obtain the nucleotide root sequence from the projected ancestral partition and store it under the `nuc` TreeIR track.
- Make MAT serialization return an actionable error when the nucleotide reference is absent, a mutation coordinate is outside it, or the reference nucleotide is not A/C/G/T.
- Remove the branch-parent fallback in [`packages/treetime-io/src/tree_ir/usher.rs#L75-L79`](../../packages/treetime-io/src/tree_ir/usher.rs#L75-L79).
- Cover ancestral, optimize, and timetree projection paths.
- Add a recurrent-mutation regression in which two mutations at the same coordinate have different parent alleles and both export the same global `ref_nuc`.
- Add parameterized writer rejection tests for a missing reference, an out-of-range coordinate, and a non-ACGT reference nucleotide, asserting actionable error context.
- Assert the exact `nuc` root sequence projected by ancestral, optimize, and timetree.

## Related issues

Source: [kb/issues/H-io-usher-ref-nuc-uses-parent-allele.md](../issues/H-io-usher-ref-nuc-uses-parent-allele.md)
