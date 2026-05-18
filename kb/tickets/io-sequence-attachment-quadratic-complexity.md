# Sequence attachment has O(n squared) complexity

The sequence attachment loop performs a linear search through all sequences for each leaf node, resulting in O(n squared) complexity where n is the number of leaves/sequences.

## Root cause

The attachment logic at [packages/treetime/src/commands/ancestral/fitch.rs#L45-L71](../../packages/treetime/src/commands/ancestral/fitch.rs#L45-L71) iterates over leaves and searches sequences:

```rust
for leaf in graph.get_leaves() {           // O(n) leaves
    let leaf_fasta = aln
        .iter()
        .find(|fasta| fasta.seq_name == leaf_name)  // O(n) sequences per leaf
}
```

The same pattern appears in dense marginal initialization at [packages/treetime/src/partition/marginal_dense.rs#L191-L195](../../packages/treetime/src/partition/marginal_dense.rs#L191-L195).

## Impact

For typical phylogenetic datasets where the number of leaves approximately equals the number of sequences:

| Sequences | Comparisons | Time (relative) |
| --------- | ----------- | --------------- |
| 100       | 10,000      | 1x              |
| 1,000     | 1,000,000   | 100x            |
| 10,000    | 100,000,000 | 10,000x         |

For large datasets (10,000+ sequences), attachment becomes a noticeable bottleneck before the actual phylogenetic algorithms run.

## Fix

Build a name index before the attachment loop:

```rust
let seq_by_name: HashMap<&str, &FastaRecord> = aln
    .iter()
    .map(|f| (f.seq_name.as_str(), f))
    .collect();

for leaf in graph.get_leaves() {
    let leaf_fasta = seq_by_name.get(leaf_name.as_str())  // O(1) lookup
        .ok_or_else(|| ...)?;
}
```

This reduces complexity to O(n) for index construction plus O(n) for attachment, giving O(n) overall.

## Additional improvement

The index construction can also detect duplicate names upfront:

```rust
let mut seq_by_name: HashMap<&str, &FastaRecord> = HashMap::new();
for fasta in aln {
    if let Some(existing) = seq_by_name.insert(&fasta.seq_name, fasta) {
        warn!("Duplicate sequence name: '{}', using last occurrence", fasta.seq_name);
    }
}
```

## Related issues

- Source: [M-io-sequence-attachment-quadratic.md](../issues/M-io-sequence-attachment-quadratic.md) -- delete after full resolution
- [Sequence-to-node name matching is unreliable](io-sequence-name-matching-unreliable.md) - the matching itself has reliability issues beyond performance
- [Dense optimize iteration is slow](../issues/N-optimize-dense-iteration-slow.md) - other performance concerns
