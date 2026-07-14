# Unified input format support for analysis commands

## Motivation

Analysis commands (`ancestral`, `timetree`, `optimize`, `clock`, `mugration`) currently accept only Newick trees and FASTA alignments as separate files. This two-file model requires name-based matching to reconcile sequences with tree nodes, which is unreliable and inefficient (see [kb/issues/M-io-sequence-name-matching-unreliable.md](../issues/M-io-sequence-name-matching-unreliable.md) and [kb/issues/M-io-sequence-attachment-quadratic.md](../issues/M-io-sequence-attachment-quadratic.md)).

The codebase already includes format adapters for unified formats where sequence data is bound to tree structure:

- Auspice JSON: edge mutations plus root sequence
- UShER MAT (protobuf): edge mutations in preorder traversal
- PhyloXML: per-clade sequence elements

Schema adapters exist in `treetime-io` and are used by the `convert` command. They do not yet define a lossless shared mutation vocabulary: UShER MAT preservation policy, PhyloXML mutation properties, and shared substitution/insertion/deletion projection remain unresolved in [kb/issues/M-io-usher-mat-mutation-loss-is-implicit.md](../issues/M-io-usher-mat-mutation-loss-is-implicit.md), [kb/issues/N-io-phyloxml-mutation-property-contract-undecided.md](../issues/N-io-phyloxml-mutation-property-contract-undecided.md), and [kb/issues/M-core-mutation-representation-and-format-projection-inconsistent.md](../issues/M-core-mutation-representation-and-format-projection-inconsistent.md). Integrating a reader therefore requires an explicit preservation contract in addition to schema parsing.

## Proposal

Allow analysis commands to accept unified input formats directly. When a unified format is provided, populate partitions from the embedded sequence/mutation data without name matching.

### Input paths

Two input paths converge to the same internal representation:

1. Newick + FASTA (existing): parse tree, parse sequences, match by name, populate partitions
2. Unified format (proposed): parse tree with embedded data, populate partitions directly

The unified path eliminates the name-matching phase entirely.

### Format capabilities

| Format         | Tree | Branch lengths | Mutations | Full sequences | Root sequence |
| -------------- | ---- | -------------- | --------- | -------------- | ------------- |
| Newick + FASTA | yes  | yes            | no        | yes (separate) | no            |
| Auspice JSON   | yes  | via div        | partial   | no             | yes           |
| UShER MAT      | yes  | no             | partial   | no             | implicit      |
| PhyloXML       | yes  | yes            | unresolved | yes (embedded) | no            |

For mutation-based formats (Auspice, MAT), full sequences are derived from root sequence plus accumulated mutations along the path from root to each node.

## Implementation approach

### Input format detection

Add a shared `--input-format` flag to analysis commands. Auto-detect from file extension when not specified.

### Unified format readers

Create adapter functions that read unified formats and return populated partition structures:

```rust
fn partitions_from_auspice(path: &Path, alphabet: &Alphabet)
    -> Result<(Graph, Vec<PartitionMarginalSparse>), Report>;

fn partitions_from_usher_mat(path: &Path, alphabet: &Alphabet)
    -> Result<(Graph, Vec<PartitionMarginalSparse>), Report>;
```

These reuse existing format reading logic from `treetime-io` but produce partition structures instead of `ConverterGraph`.

### Sequence derivation for mutation-based formats

For Auspice and MAT inputs, derive full sequences by traversing from root:

1. Start with root sequence (explicit in Auspice, implicit reference in MAT)
2. For each node, apply accumulated mutations along path from root
3. Store derived sequence in partition node data

This allows marginal reconstruction algorithms to operate on the derived sequences.

### Sparse-native mode

For workflows that only need mutations (not full probability distributions), operate directly on the mutation representation without deriving full sequences. This matches how UShER processes pandemic-scale datasets.

## Interaction with existing proposals

- [config-file-multi-partition.md](config-file-multi-partition.md): the config file could specify input format per partition, mixing Newick+FASTA with unified formats
- [kb/decisions/multi-format-tree-io.md](../decisions/multi-format-tree-io.md): documents the existing format adapters that this proposal would integrate

## Open questions

1. Should unified format input be a separate flag (`--input-format auspice`) or inferred from file extension?

2. For mutation-based formats without explicit root sequence, what is the fallback? Require `--root-sequence` flag?

3. Should output format match input format by default? If input is Auspice, should output also be Auspice with updated annotations?

4. How to handle format-specific metadata (Auspice colorings, MAT clade annotations) through the analysis pipeline?

## Recommended complete design

The production design should combine all of the following:

- Index Newick-plus-FASTA input by validated unique names and report duplicate or missing mappings precisely.
- Parse each unified format into a typed representation whose preservation and unsupported-state contracts are explicit before constructing partitions.
- Configure input format per partition through [config-file-multi-partition.md](config-file-multi-partition.md), including multi-segment datasets and explicit name mappings.
- Support indexed or memory-mapped sequence access so the same input contract remains usable for datasets that do not fit in memory.

These capabilities share the same typed partition-construction boundary. None substitutes for the others.

## Orthogonal concerns

The question of whether to remove payloads from graph nodes (making graph topology-only) is **orthogonal** to input format. If graph becomes topology-only:

- Input still needs to populate data stores
- Name-based matching still needed for FASTA input
- The input format problem remains unchanged

The graph payload question affects internal organization. This proposal affects input paths. They can be addressed independently.

## Related issues

- [kb/issues/M-io-sequence-name-matching-unreliable.md](../issues/M-io-sequence-name-matching-unreliable.md)
- [kb/issues/M-io-sequence-attachment-quadratic.md](../issues/M-io-sequence-attachment-quadratic.md)
- [kb/issues/N-io-large-dataset-memory-constraint.md](../issues/N-io-large-dataset-memory-constraint.md)
- [kb/issues/M-core-mutation-representation-and-format-projection-inconsistent.md](../issues/M-core-mutation-representation-and-format-projection-inconsistent.md)
- [kb/issues/M-timetree-tree-output-inference-metadata-incomplete.md](../issues/M-timetree-tree-output-inference-metadata-incomplete.md)
- [kb/issues/N-io-multi-segment-genome-input.md](../issues/N-io-multi-segment-genome-input.md)

## Related documentation

- [kb/decisions/multi-format-tree-io.md](../decisions/multi-format-tree-io.md) - existing format adapter implementations
- [kb/decisions/partition-system-architecture.md](../decisions/partition-system-architecture.md) - target internal representation
