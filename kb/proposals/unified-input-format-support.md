# Unified input format support for analysis commands

## Motivation

Analysis commands (`ancestral`, `timetree`, `optimize`, `clock`, `mugration`) currently accept only Newick trees and FASTA alignments as separate files. This two-file model requires name-based matching to reconcile sequences with tree nodes, which is unreliable and inefficient (see [M-io-sequence-name-matching-unreliable](../issues/M-io-sequence-name-matching-unreliable.md) and [M-io-sequence-attachment-quadratic](../issues/M-io-sequence-attachment-quadratic.md)).

The codebase already includes format adapters for unified formats where sequence data is bound to tree structure:

- Auspice JSON: edge mutations plus root sequence
- UShER MAT (protobuf): edge mutations in preorder traversal
- PhyloXML: per-clade sequence elements

These adapters are implemented in `treetime-io` and used by the `convert` command, but not integrated into analysis commands.

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
| Auspice JSON   | yes  | via div        | yes       | no             | yes           |
| UShER MAT      | yes  | no             | yes       | no             | implicit      |
| PhyloXML       | yes  | yes            | no        | yes (embedded) | no            |

For mutation-based formats (Auspice, MAT), full sequences are derived from root sequence plus accumulated mutations along the path from root to each node.

## Implementation approach

### Phase 1: Input format detection

Add `--input-format` flag to analysis commands, mirroring the `convert` command at [packages/treetime-cli/src/convert/args.rs#L7-L18](../../packages/treetime-cli/src/convert/args.rs#L7-L18). Auto-detect from file extension when not specified.

### Phase 2: Unified format readers

Create adapter functions that read unified formats and return populated partition structures:

```rust
fn partitions_from_auspice(path: &Path, alphabet: &Alphabet)
    -> Result<(Graph, Vec<PartitionMarginalSparse>), Report>;

fn partitions_from_usher_mat(path: &Path, alphabet: &Alphabet)
    -> Result<(Graph, Vec<PartitionMarginalSparse>), Report>;
```

These reuse existing format reading logic from `treetime-io` but produce partition structures instead of `ConverterGraph`.

### Phase 3: Sequence derivation for mutation-based formats

For Auspice and MAT inputs, derive full sequences by traversing from root:

1. Start with root sequence (explicit in Auspice, implicit reference in MAT)
2. For each node, apply accumulated mutations along path from root
3. Store derived sequence in partition node data

This allows marginal reconstruction algorithms to operate on the derived sequences.

### Phase 4: Sparse-native mode

For workflows that only need mutations (not full probability distributions), operate directly on the mutation representation without deriving full sequences. This matches how UShER processes pandemic-scale datasets.

## Interaction with existing proposals

- [Configuration file format for multi-partition analysis](config-file-multi-partition.md): the config file could specify input format per partition, mixing Newick+FASTA with unified formats
- [Multi-format tree I/O](../decisions/multi-format-tree-io.md): documents the existing format adapters that this proposal would integrate

## Open questions

1. Should unified format input be a separate flag (`--input-format auspice`) or inferred from file extension?

2. For mutation-based formats without explicit root sequence, what is the fallback? Require `--root-sequence` flag?

3. Should output format match input format by default? If input is Auspice, should output also be Auspice with updated annotations?

4. How to handle format-specific metadata (Auspice colorings, MAT clade annotations) through the analysis pipeline?

## Implementation priority

Recommended order for addressing input architecture:

1. Immediate: Improve Newick+FASTA attachment (A0 in design). Build name index, detect duplicates, clear errors. This is the primary input path.

2. Short-term: Prototype Auspice JSON input for one analysis command. Validates the unified-format approach with minimal risk.

3. Medium-term: Config file format per [config-file-multi-partition](config-file-multi-partition.md). Enables multi-segment and explicit name mappings.

4. Long-term: Out-of-core processing for large datasets (indexed FASTA, memory-mapped access). Requires architectural changes.

## Orthogonal concerns

The question of whether to remove payloads from graph nodes (making graph topology-only) is **orthogonal** to input format. If graph becomes topology-only:

- Input still needs to populate data stores
- Name-based matching still needed for FASTA input
- The input format problem remains unchanged

The graph payload question affects internal organization. This proposal affects input paths. They can be addressed independently.

## Related issues

- [Sequence-to-node name matching is unreliable](../issues/M-io-sequence-name-matching-unreliable.md)
- [Sequence attachment has O(n squared) complexity](../issues/M-io-sequence-attachment-quadratic.md)
- [Large datasets require all sequences in memory](../issues/N-io-large-dataset-memory-constraint.md)
- [Auspice JSON output incomplete](../issues/N-timetree-auspice-json-incomplete.md)
- [Multi-segment genome input not wired](../issues/N-io-multi-segment-genome-input.md)

## Related documentation

- [Multi-format tree I/O](../decisions/multi-format-tree-io.md) - existing format adapter implementations
- [Partition system architecture](../decisions/partition-system-architecture.md) - target internal representation
