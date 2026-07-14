# UShER MAT export fabricates global reference nucleotides from branch-local alleles

UShER `message mut` [`packages/util-usher-mat/schemas/parsimony.proto#L4-L11`](../../packages/util-usher-mat/schemas/parsimony.proto#L4-L11) defines `ref_nuc` as the reference nucleotide at a mutation's position. TreeTime's mutation-bearing command projections leave `TreeIrData.root_sequence` empty:

- ancestral constructs default `TreeIrData` with only `has_mutations` set in [`packages/treetime/src/commands/ancestral/run.rs#L243-L263`](../../packages/treetime/src/commands/ancestral/run.rs#L243-L263);
- optimize does the same in [`packages/treetime/src/commands/optimize/run.rs#L84-L102`](../../packages/treetime/src/commands/optimize/run.rs#L84-L102);
- timetree constructs default root-sequence data in [`packages/treetime/src/commands/timetree/output/ir.rs#L22-L28`](../../packages/treetime/src/commands/timetree/output/ir.rs#L22-L28).

`TreeIrUsherWriter` looks up the reference nucleotide and falls back to the mutation's branch-local parent allele when the reference is absent or cannot encode the position. [`packages/treetime-io/src/tree_ir/usher.rs#L75-L84`](../../packages/treetime-io/src/tree_ir/usher.rs#L75-L84)

For a recurrent mutation, parent alleles at one coordinate can differ across branches. The exported MAT can therefore assign different `ref_nuc` values to mutations at the same coordinate even though the field denotes one global reference.

## Required behavior

- Populate the nucleotide root reference when projecting a mutation-bearing domain graph to TreeIR.
- Reject MAT export when the reference is absent, the coordinate is outside the reference, or the reference nucleotide is not representable as A/C/G/T.
- Verify a recurrent-mutation case where a branch parent allele differs from the global reference.

## Related issues

- [N-io-tree-ir-architecture-unapproved.md](N-io-tree-ir-architecture-unapproved.md)

## Related tickets

- [kb/tickets/io-populate-and-validate-usher-reference-sequence.md](../tickets/io-populate-and-validate-usher-reference-sequence.md)
