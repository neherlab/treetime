# Large datasets require all sequences in memory simultaneously

The current input model loads all sequences into memory before attachment, then copies them into partition data structures. For very large datasets, this memory requirement may exceed available RAM.

## Current flow

Sequence loading at [packages/treetime/src/commands/ancestral/run.rs#L61-L69](../../packages/treetime/src/commands/ancestral/run.rs#L61-L69):

```rust
let aln = if input_fastas.is_empty() {
    // read from stdin
} else {
    read_many_fasta(input_fastas, &alphabet)?  // loads ALL sequences
};
```

The entire `Vec<FastaRecord>` must fit in memory. After attachment, sequences are copied/transformed into partition data structures, then the original vector is dropped.

## Memory characteristics

1. **Peak memory**: approximately 2x sequence data (FASTA records + partition copies)
2. **Post-attachment**: original `Vec<FastaRecord>` is dropped, partitions remain
3. **Algorithm phase**: partitions must remain in memory for random access during tree traversal

## Why streaming is blocked

True streaming (process each sequence once as it arrives) is blocked by:

1. **Length validation**: `get_common_length()` at [packages/treetime/src/commands/ancestral/fitch.rs#L520](../../packages/treetime/src/commands/ancestral/fitch.rs#L520) verifies all sequences have uniform length before processing begins

2. **Tree-order access**: algorithms traverse the tree in topological order (postorder, then preorder). Leaf data is accessed in tree-traversal order, not FASTA file order

3. **Random access during traversal**: backward and forward passes access arbitrary leaves based on tree structure

## Scale considerations

| Sequences | Genome size | Raw data | Peak memory (approx) |
| --------- | ----------- | -------- | -------------------- |
| 1,000     | 30 kb       | 30 MB    | 60 MB                |
| 10,000    | 30 kb       | 300 MB   | 600 MB               |
| 100,000   | 30 kb       | 3 GB     | 6 GB                 |
| 1,000,000 | 30 kb       | 30 GB    | 60 GB                |
| 10,000    | 4.4 Mb (TB) | 44 GB    | 88 GB                |

Pandemic-scale datasets (millions of SARS-CoV-2 genomes) or large bacterial genomes (tuberculosis at 4.4 Mb) can exceed typical workstation memory.

## Possible approaches

### Indexed FASTA with lazy loading

First pass: scan FASTA, build index mapping sequence name to file byte offset. Store only the index in memory. During attachment: seek to offset, read single sequence on demand.

Requires changes to `FastaReader` API to support seeking.

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

- [Sequence attachment has O(n squared) complexity](M-io-sequence-attachment-quadratic.md) - attachment performance
- [Multi-segment genome input not wired](N-io-multi-segment-genome-input.md) - related input architecture

## Related documentation

- [Multi-format tree I/O](../decisions/multi-format-tree-io.md) - UShER MAT format stores mutations, not full sequences
