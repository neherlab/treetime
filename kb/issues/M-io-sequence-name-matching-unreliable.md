# Sequence-to-node name matching is unreliable

The current input model reads tree topology from Newick and sequences from FASTA as separate files, then reconciles them by matching sequence names to leaf node names. This name-based matching is inherently fragile.

## Failure modes

### Duplicate leaf names

When a tree contains multiple nodes with the same name (e.g., two leaves named `"USA"`), the attachment logic at [packages/treetime/src/ancestral/fitch.rs#L53-L57](../../packages/treetime/src/ancestral/fitch.rs#L53-L57) uses `.find()` which returns the first match. The second node silently receives the same sequence or fails.

### Duplicate FASTA names

When a FASTA file contains multiple sequences with the same name, the first matching sequence is used. Subsequent sequences with that name are ignored without warning.

### Unnamed internal nodes

Newick internal nodes often have no name or arbitrary labels like `"Node_42"`. There is no mechanism to attach ancestral sequences from external sources to these nodes, as name matching requires names.

### Whitespace and encoding mismatches

Names that differ only in whitespace or encoding are treated as distinct:

- `"Sample 1"` vs `"Sample_1"` (space vs underscore)
- `"Sample  1"` (double space) vs `"Sample 1"` (single space)
- Trailing whitespace, invisible Unicode characters
- Different Unicode normalization forms (NFC vs NFD)

### Case sensitivity

The comparison at [packages/treetime/src/ancestral/fitch.rs#L55](../../packages/treetime/src/ancestral/fitch.rs#L55) is case-sensitive. `"USA"` and `"usa"` are treated as different names with no standard convention for which is correct.

## Current behavior

Attachment fails with an error when a leaf node name has no matching sequence:

```rust
.ok_or_else(|| make_report!("Leaf sequence not found: '{leaf_name}'"))?;
```

This error is actionable but provides no guidance on near-matches or potential causes (case, whitespace, duplicates).

## Impact

Users with mismatched names must manually edit either the tree or the FASTA file to make names match exactly. The error message does not help identify which names are close matches or suggest corrections.

## Possible improvements

1. Build a name index upfront and report all mismatches at once (not one at a time)
2. Detect and report duplicate names before attachment
3. Provide fuzzy matching suggestions for near-misses
4. Support explicit name mapping via config file (see [config-file-multi-partition proposal](../proposals/config-file-multi-partition.md))
5. Support case-insensitive matching as an option

## Related issues

- [Column auto-detection gaps in CSV readers](M-dates-column-auto-detection-gaps.md) - similar name matching issues for metadata columns
- [Multi-segment genome input not wired](N-io-multi-segment-genome-input.md) - related input architecture gap
- [Optimize command accepts only a single alignment](N-optimize-multi-alignment-input.md) - related input limitation

## Related

- [N-io-name-reconciliation-duplicated.md](N-io-name-reconciliation-duplicated.md) - reconciliation logic duplicated across subsystems
- Design: [kb/proposals/input-name-matching-validation.md](../proposals/input-name-matching-validation.md) - ecosystem survey, four design axes, architecture analysis
- [Configuration file format for multi-partition analysis](../proposals/config-file-multi-partition.md) - could include explicit name mappings
