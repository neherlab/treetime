# Amino-acid outputs disagree with the translated nucleotide outputs

`treetime ancestral --translations <peptides> --annotation <gff3>` reconstructs amino acids in one independent partition per CDS, from tip peptides, with its own GTR model and fitted rate per CDS (`fn reconstruct_aa()` [packages/treetime/src/ancestral/aa.rs#L25-L92](../../packages/treetime/src/ancestral/aa.rs#L25-L92)). The amino-acid output of a node can therefore differ from the translation of that node's nucleotide output:

- A node can carry an amino acid that its own codon does not code
- A branch can carry an amino-acid mutation without any nucleotide change in that codon

Auspice derives ancestral states from the "from" letters of the mutations on the path from the root and only logs `console.error` on a mismatch, so such chains give wrong colors without a visible error. In Nextstrain pipelines the peptides are Nextclade translations, which Nextclade aligns to the reference peptide, so peptide columns and CDS codons also diverge near indels and frameshifts.

TreeTime has no translator of its own: every amino-acid output depends on peptides from an external tool.

## Evidence

Release build, default settings, the bundled peptide files (translations of the same alignment, so peptide columns match the CDS codons). Internal nodes only, terminal stop codon excluded:

| Dataset              | Internal-node residues | Both determined, different | One side `X` or `-`, different |
| -------------------- | ---------------------- | -------------------------- | ------------------------------ |
| `rsv/a/20`           | 82,026                 | 11                         | 3                              |
| `rsv/a/2000`         | 9,100,329              | 76                         | 496                            |
| `mpox/clade-ii/2000` | 109,640,673            | 17                         | 4,915                          |

Model differences at confident codons cause most both-determined differences; the terminal-gap fill of each peptide causes most of the others ([M-ancestral-terminal-gap-fill-acts-on-each-peptide.md](M-ancestral-terminal-gap-fill-acts-on-each-peptide.md)). Method, commands and analysis: section "Measured disagreement of the two routes" of [kb/reports/codon-aware-ancestral-reconstruction.md](../reports/codon-aware-ancestral-reconstruction.md#measured-disagreement-of-the-two-routes).

## Fix direction

[kb/proposals/codon-aware-ancestral-reconstruction.md](../proposals/codon-aware-ancestral-reconstruction.md) recommends that TreeTime translates its own reconstruction:

- **Coding view**: each CDS maps to its ordered alignment columns (strand, segments, overlaps, frameshift junctions)
- **Translator** with the Nextclade rules: an ambiguous codon whose resolutions all code one amino acid gives that amino acid, any other ambiguous or partly gapped codon gives `X`, `---` gives `-`. Tips with frameshifting indels need detection and masking, as Nextclade does, or a documented limit
- **Decoder**: the codon posterior is the product of the three site posteriors, exact under the current site-independent models. Choose the most probable amino acid with `*` as a 21st symbol, then its codon with the largest sum of site posteriors. Decode each overlap region jointly over all its frames by dynamic programming. A codon with one or two gaps at an internal node keeps its gap pattern, decodes its other positions per site, and translates as `X`
- **One decoding step before every writer**, so FASTA, node data, Nexus, Auspice and MAT outputs read the same codons; imputed tips, ambiguity codes and sampled sequences use the same rule
- **`--translations` only for protein-only input**: reject a run that gives the same CDS as nucleotides and as peptides. [M-ancestral-aa-path-ignores-method-and-include-leaves.md](M-ancestral-aa-path-ignores-method-and-include-leaves.md) then applies only to protein-only input

> [!IMPORTANT]
> **Decision required.** The fix changes scientific outputs, and each part needs approval:
>
> - **Source of amino acids**: translate TreeTime's own decoded codons, keep the independent reconstruction, or both with a consistency check
> - **Decoding rule**: amino acid first (lowest expected amino-acid error; changes nucleotides at diffuse coding sites against v0 and augur, which needs an entry in `kb/decisions/`), nucleotide first (smallest change of nucleotide output), or most probable codon
> - **Stop codons at internal nodes**: `*` as a symbol (keeps real ancestral stop codons such as SARS-CoV-2 ORF8 Q27*), sense codons only, or exact conditioning on sense codons
> - **Overlapping frames**: joint decoding over all frames, per-site argmax in overlaps, or a primary frame per overlap
> - **Translator**: `X` for codons with several meanings (Auspice ignores only `X`), or Biopython's `B`, `Z` and `J`
>
> The proposal gives options, trade-offs and evidence for each (axes 1 to 5 and 7).

> [!IMPORTANT]
> **Investigation required.** A spliced codon inside an overlap region (influenza M2, NEP) breaks the chain structure of the overlap decoder; the decoder needs state carried across the junction or a small junction tree, and no test covers it. The share of diffuse codons, where the decoding rules differ, is not measured beyond two shallow datasets.

## Validation

- Decoder tests defined by the specification: the worked examples of the report, and properties over random posteriors (the translation of the decoded nucleotides equals the decoded amino acid; the decoder picks the most probable amino acid; it keeps the per-site argmax whenever that codes the most probable amino acid; a stop codon is decoded only when the stop mass is the largest)
- Brute-force oracles: all ancestral triplet assignments on a small tree against the product of site posteriors; all nucleotide strings of short overlap regions against the overlap decoder
- The route comparison of the report repeated after the change: zero differences between translated nucleotides and amino acids at every node

## Related

- [M-ancestral-node-data-annotation-drops-cds-phase.md](M-ancestral-node-data-annotation-drops-cds-phase.md), [M-io-gff3-cds-segment-order-wrong-across-circular-origin.md](M-io-gff3-cds-segment-order-wrong-across-circular-origin.md), [M-ancestral-translation-annotation-not-checked-against-alignment.md](M-ancestral-translation-annotation-not-checked-against-alignment.md): the annotation reader that the coding view needs
- [M-ancestral-sample-from-profile-draws-nodes-independently.md](M-ancestral-sample-from-profile-draws-nodes-independently.md): sampled sequences
- [M-ancestral-codon-substitution-models-unimplemented.md](M-ancestral-codon-substitution-models-unimplemented.md): the optional codon model, which this fix does not need

## Proposed implementation plan

The steps below need no answer of the Decision required block above; the decoder and writer steps wait for it.

1. **Coding view** (no decision needed): a type in `packages/treetime` that maps each CDS of the GFF3 annotation (read by `packages/treetime-io/src/gff.rs`) to its ordered alignment columns, with strand, segments, phase, overlaps and frameshift junctions. Depends on the annotation-reader issues under "Related"
2. **Translator** (needs the "Translator" answer): a pure function from a codon of three `AsciiChar` to an amino acid, with the standard genetic code as a static table. No translator exists in v1 yet (only `packages/util-newick/src/nexus/read.rs` has an unrelated Nexus `translate`)
3. **Decoder** (needs the "Decoding rule", "Stop codons" and "Overlapping frames" answers): a pure function from the per-site posteriors of the nucleotide partition and the coding view to decoded codons and amino acids, one decoding step before every writer of `treetime ancestral`
4. **Route change** (needs the "Source of amino acids" answer): replace or check `fn reconstruct_aa()` [packages/treetime/src/ancestral/aa.rs#L25-L92](../../packages/treetime/src/ancestral/aa.rs#L25-L92) for CDS that also have nucleotides, and reject `--translations` for such CDS when the answer is "translate own codons"
5. **Decision file**: record the approved rules in `kb/decisions/`, because the decoder changes scientific outputs

Tests: the tests listed under "Validation" above, in `packages/treetime/src/ancestral/__tests__/` (translator unit table with the standard genetic code as oracle; decoder property tests and brute-force oracles). Validate with `just check`, `just test-rs`, and the route comparison of [kb/reports/codon-aware-ancestral-reconstruction.md](../reports/codon-aware-ancestral-reconstruction.md#measured-disagreement-of-the-two-routes) on `rsv/a/20`, `rsv/a/2000` and `mpox/clade-ii/2000` (target: zero differences).
