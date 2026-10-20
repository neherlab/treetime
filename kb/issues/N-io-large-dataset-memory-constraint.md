# Large datasets require all sequences in memory simultaneously

The current input model loads all sequences into memory before attachment, then copies them into partition data structures. For very large datasets, this memory requirement may exceed available RAM.

## Current flow

Sequence loading in `fn read_sequence_source()` [packages/app-commands/src/commands/shared/sequence_inputs.rs#L84-L124](../../packages/app-commands/src/commands/shared/sequence_inputs.rs#L84-L124), called once per `--alignment` file:

```rust
let mut records = fasta_read_file(path, alphabet)?;
```

All records of a file must fit in memory, and the paired sequences of every file stay in memory until the partitions are built. After attachment, sequences are copied/transformed into partition data structures, then the original vector is dropped.

## Memory characteristics

> [!IMPORTANT]
> **Investigation required.** The peak of twice the sequence data (FASTA records plus partition copies) is an estimate from the code structure, not a measurement. Measure the peak resident memory of a run on a large dataset (for example `data/sc2/4500` or `data/mpox/clade-ii/2000`) against the size of the decompressed alignment, for dense and sparse partitions.

1. Peak memory: estimated at 2x sequence data (FASTA records + partition copies)
2. Post-attachment: original `Vec<FastaRecord>` is dropped, partitions remain
3. Algorithm phase: partitions must remain in memory for random access during tree traversal

## Why streaming is blocked

True streaming (process each sequence once as it arrives) is blocked by:

1. Length validation: `LeafSequences::common_length()` [packages/treetime/src/seq/alignment.rs#L91-L93](../../packages/treetime/src/seq/alignment.rs#L91-L93), through the free function `common_length()` [packages/treetime/src/seq/alignment.rs#L128](../../packages/treetime/src/seq/alignment.rs#L128), verifies that all sequences of a file have uniform length before processing begins

2. Tree-order access: algorithms traverse the tree in topological order (postorder, then preorder). Leaf data is accessed in tree-traversal order, not FASTA file order

3. Random access during traversal: backward and forward passes access arbitrary leaves based on tree structure

## Scale considerations

| Sequences | Genome size | Raw data | Peak memory (2x estimate) |
| --------- | ----------- | -------- | -------------------- |
| 1,000     | 30 kb       | 30 MB    | 60 MB                |
| 10,000    | 30 kb       | 300 MB   | 600 MB               |
| 100,000   | 30 kb       | 3 GB     | 6 GB                 |
| 1,000,000 | 30 kb       | 30 GB    | 60 GB                |
| 10,000    | 4.4 Mb (TB) | 44 GB    | 88 GB                |

Pandemic-scale datasets (millions of SARS-CoV-2 genomes) or large bacterial genomes (tuberculosis at 4.4 Mb) can exceed typical workstation memory.

Public reports ask about TreeTime operation with approximately 10 million tips [[issue](https://github.com/neherlab/treetime/issues/207)] and 12.68 million tips [[issue](https://github.com/neherlab/treetime/issues/343)]. A maintainer stated that the Python implementation would not work with more than a few thousand tips [[comment](https://github.com/neherlab/treetime/issues/207#issuecomment-1240512184)]. These reports establish the user-facing scale requirement; the Rust peak-memory mechanism described here is implementation-specific.

## Possible approaches

### Indexed FASTA with lazy loading

First pass: scan FASTA, build index mapping sequence name to file byte offset. Store only the index in memory. During attachment: seek to offset, read single sequence on demand.

Requires a FASTA reader that supports seeking; `fasta_read()` reads the whole input into records.

### Memory-mapped FASTA

Use `memmap2` crate to memory-map the FASTA file. The operating system handles paging, and sequences are accessed as slices without explicit loading.

Requires parsing directly from memory-mapped slices.

### Sparse representation throughout

For mutation-annotated inputs (UShER MAT, Auspice JSON), never materialize full sequences. Store only root sequence plus per-edge mutations. Derive sequence at any node by traversing the path from root.

This is how UShER handles millions of SARS-CoV-2 genomes in compact form.

### External storage

Store sequences in a key-value database (SQLite, RocksDB) indexed by name. Random access without full dataset in memory.

Adds external dependency and I/O overhead.

## Related issues

- [Multi-segment genome input not wired](N-io-multi-segment-genome-input.md) - related input architecture

## Related documentation

- [kb/decisions/multi-format-tree-io.md](../decisions/multi-format-tree-io.md) - UShER MAT format stores mutations, not full sequences
