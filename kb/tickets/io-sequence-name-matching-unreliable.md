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

## Tests

### Happy paths

- All names match exactly: attachment succeeds, no warnings
- One mismatch: error lists the unmatched name with fuzzy suggestions
- Multiple mismatches: error lists all unmatched names at once (not one-at-a-time)
- Case-only mismatch: `USA` in tree, `usa` in FASTA -> suggestion "did you mean 'usa'?"
- Underscore/space swap: `Sample_1` in tree, `Sample 1` in FASTA -> suggestion
- Extra whitespace: `name ` (trailing space) vs `name` -> suggestion
- Fuzzy suggestion ranking: closest match listed first

### Edge cases -- duplicates

- Duplicate leaf names in tree: `(USA:0.1,USA:0.2);` -> error listing duplicates before attachment
- Duplicate FASTA names: two sequences named `A` -> error or warning listing duplicates
- Duplicate in both tree and FASTA: report both sets independently
- Leaf name appears once in tree but twice in FASTA -> warning (ambiguous which sequence to use)

### Edge cases -- near-misses

- Double space vs single: `Sample  1` vs `Sample 1` -> suggestion
- Tab vs space: `Sample\t1` vs `Sample 1` -> suggestion
- Unicode normalization: NFC vs NFD of accented characters (e.g. `e` + combining acute vs `e-acute`) -> suggestion
- Invisible Unicode: zero-width space, BOM, soft hyphen in name -> suggestion or strip
- Substring match: `Sample_1` in tree, `Sample_1_duplicate` in FASTA -> no false positive (substring is not a match)
- Swapped names: tree has `A` and `B`, FASTA has `B` and `A` -> exact match (order doesn't matter)

### Edge cases -- boundary

- Tree with 1 leaf, FASTA with 1 sequence, names match -> success
- Tree with 1 leaf, FASTA with 0 sequences -> error
- Tree with 0 external sequences needed (all internal) -> success (nothing to attach)
- Empty leaf name in tree -> error (cannot match empty name)
- FASTA name that is only whitespace -> treated as empty, reported

### Pathological

- Tree with 10K+ leaves, all mismatched -> error message lists all, performance of fuzzy matching acceptable (not O(n^2) against FASTA names)
- FASTA with 10K+ sequences, 1 mismatch -> single mismatch found without scanning entire FASTA for each leaf
- Names that are 1000+ characters -> fuzzy matching still works, no truncation of suggestions
- Names containing only special characters: `---`, `***`, `[brackets]`

### Regression

- Existing tests with matching names still pass
- Error format change: tests checking exact error message text need updating (list this explicitly)

## Related issues

- Source: [M-io-sequence-name-matching-unreliable.md](../issues/M-io-sequence-name-matching-unreliable.md) -- delete after full resolution
- [Column auto-detection gaps in CSV readers](../issues/M-dates-column-auto-detection-gaps.md) - similar name matching issues for metadata columns
- [Multi-segment genome input not wired](../issues/N-io-multi-segment-genome-input.md) - related input architecture deficiency
- [Optimize command accepts only a single alignment](../issues/N-optimize-multi-alignment-input.md) - related input limitation

## Related proposals

- [Configuration file format for multi-partition analysis](../proposals/config-file-multi-partition.md) - could include explicit name mappings
