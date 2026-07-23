# PartitionCompressed exposes the concrete Fitch storage layout

`PartitionCompressed` returns concrete `BTreeMap<GraphNodeKey, SparseNodePartition>` and `BTreeMap<GraphEdgeKey, SparseEdgePartition>` values in read-only, mutable, and paired mutable forms [`packages/treetime/src/partition/traits.rs#L227`](../../packages/treetime/src/partition/traits.rs#L227). Its only implementation is `PartitionFitch` [`packages/treetime/src/partition/fitch.rs#L54`](../../packages/treetime/src/partition/fitch.rs#L54).

Clients perform map insertion, replacement, and in-place mutation themselves. Any alternate compressed representation would have to imitate the same maps, and the partition cannot enforce invariants across node and edge storage.

## Evidence of speculative generality

- The implementation count is one: `PartitionFitch`.
- No test double or second compressed representation implements the trait.
- Principal consumers are explicitly Fitch-specific: sequence compression, Fitch GTR inference, and Fitch mutation counting.
- Tests construct and inspect `PartitionFitch` directly rather than substituting another implementation.

The trait does not hide storage variation because it exposes the exact map keys, values, and mutability. Replacing the maps would break every caller despite the nominal abstraction.

## Decision constraint

The current trait is documented in [kb/decisions/partition-system-architecture.md](../decisions/partition-system-architecture.md). Removing or replacing it requires explicit user approval. The open question is whether Fitch-specific inherent operations are sufficient or whether a second representation demonstrates a narrower behavioral abstraction.

## Options

- O1. Move storage access to inherent `PartitionFitch` operations and accept `PartitionFitch` in Fitch-specific algorithms. Introduce no abstraction until a second representation supplies concrete shared behavior.
- O2. Replace map access with behavioral messages for compression, node/edge initialization, mutation extraction, and counting. Retain a trait only if more than one representation can honor those messages.

O1 is supported by the current implementation set, but it conflicts with the recorded architecture and cannot proceed without approval.

## Validation

- Inventory every current map mutation and identify the invariant it serves.
- Fitch compression, GTR inference, reconstruction, and mutation-count tests retain whole-value behavior.
- Alternate partition types are not forced to mimic Fitch maps.
- No compatibility trait or re-export remains because v1 has not shipped.

## Related issues

- [M-core-partition-init-orchestration-duplication.md](M-core-partition-init-orchestration-duplication.md)
- [N-representation-dense-sparse-partition-asymmetry.md](N-representation-dense-sparse-partition-asymmetry.md)
