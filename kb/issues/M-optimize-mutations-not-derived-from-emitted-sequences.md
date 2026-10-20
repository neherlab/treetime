# Optimize writes engine substitutions instead of mutations derived from node sequences

`treetime optimize` writes the substitutions of the marginal reconstruction engine: at each position, the most likely parent state against the most likely child state. The `ancestral` and `timetree` commands derive their mutations from the emitted node sequences instead, so a leaf-branch mutation reports the observed tip letter ([kb/decisions/ancestral-marginal-tip-reconstruction-and-imputation.md](../decisions/ancestral-marginal-tip-reconstruction-and-imputation.md)). The two sibling commands therefore write different mutations for the same leaf.

## Example

A parent `C` and a tip that observes `Y` (`C` or `T`) whose posterior favors `T`:

| Command                 | Written mutation | `--divergence-units=mutations` count |
| ----------------------- | ---------------- | ------------------------------------ |
| `ancestral`, `timetree` | `C1Y`            | 0 (`C` and `Y` share `C`)            |
| `optimize`              | `C1T`            | 1                                    |

## Mechanism

- `fn optimize_result()` builds the partition result of each partition from `ReconstructionRef::edge_mutations()` for every edge ([packages/treetime/src/optimize/pipeline.rs#L199-L233](../../packages/treetime/src/optimize/pipeline.rs#L199-L233))
- `fn optimize_output_maps()` passes these mutations through the unknown-mutation filter to the outputs ([packages/app-commands/src/commands/optimize/run.rs#L170-L185](../../packages/app-commands/src/commands/optimize/run.rs#L170-L185))
- `ReconstructionRef::edge_mutations()` combines the engine substitutions with the edge indels ([packages/treetime/src/partition/marginal/reconstruction.rs#L245-L252](../../packages/treetime/src/partition/marginal/reconstruction.rs#L245-L252))
- The dense engine compares the argmax states of the parent and child profiles and skips non-canonical states (`fn edge_subs()` in [packages/treetime/src/partition/marginal/dense/partition.rs#L80-L111](../../packages/treetime/src/partition/marginal/dense/partition.rs#L80-L111)), so the unknown-mutation filter that optimize applies never removes anything

## Fix direction

Derive the optimize mutations in `fn optimize_result()` through `stream_sequences()` in `packages/treetime/src/partition/marginal/reconstruction.rs`, as the `ancestral` and `timetree` commands do, so the root, the mutations and the mutation counts of all three commands come from one source and one rule. That method compares the emitted node sequences (`fn stream_sequence_mutations()` in `packages/treetime/src/seq/mutation.rs`) for dense reconstructions and derives the same result from the sparse node states without materializing the sequences (`fn sparse_edge_mutations()` in `packages/treetime/src/partition/marginal/sparse/mutations.rs`) for sparse ones.
