# Nexus tree and PHYLIP alignment input are not supported

The commands read trees only as Newick and alignments only as FASTA. A Nexus tree fails in the Newick parser, and a PHYLIP alignment fails in the FASTA reader:

```
Error:
   0: Failed to load tree from file
   1: When reading file 'data/flu/h3n2/20/tree.nex'
   2: Failed to parse Newick string
      2 | Begin Taxa;

Error:
   0: FASTA input is incorrectly formatted: expected at least one FASTA record starting with character '>', but none found
```

v0 reads both formats through Biopython (`Bio.Phylo` for trees, `Bio.AlignIO` for alignments), and the bundled dataset `data/flu/h3n2/20` ships `tree.nex` and `aln.phylip` for them.

## Scope

- `packages/treetime-io/src/nex.rs` has Nexus writers only; no reader exists for Nexus trees or PHYLIP alignments
- `kb/features/io.md` listed "Nexus tree" and "Phylip tree" as implemented input formats; both are now marked not done

## Reproduction

```bash
treetime ancestral --tree=data/flu/h3n2/20/tree.nex --aln=data/flu/h3n2/20/aln.fasta.xz --output-all=<dir>
treetime ancestral --tree=data/flu/h3n2/20/tree.nwk --aln=data/flu/h3n2/20/aln.phylip --output-all=<dir>
```

These are the `ancestral/flu/h3n2/20/input-nexus-tree`, `ancestral/flu/h3n2/20/input-phylip-alignment` and `timetree/flu/h3n2/20/input-nexus-tree` cases of `dev/smoke`, declared expected failures in `dev/smoke.toml`.

## Decision needed

Implement both readers for parity with v0, or record the reduced input set as a decision in `kb/decisions/` and remove the bundled `tree.nex` and `aln.phylip`.
