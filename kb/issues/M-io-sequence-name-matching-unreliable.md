# Sequence-to-node name matching is unreliable

The current input model reads tree topology from Newick and sequences from FASTA as separate files, then reconciles them by matching sequence names to leaf node names. This name-based matching is inherently fragile.

## Failure modes

### Unnamed internal nodes

Newick internal nodes often have no name or arbitrary labels like `"Node_42"`. There is no mechanism to attach ancestral sequences from external sources to these nodes, as name matching requires names.

### Whitespace and encoding mismatches

Names that differ only in whitespace or encoding are treated as distinct:

- `"Sample 1"` vs `"Sample_1"` (space vs underscore)
- `"Sample  1"` (double space) vs `"Sample 1"` (single space)
- Trailing whitespace, invisible Unicode characters
- Different Unicode normalization forms (NFC vs NFD)

### Case sensitivity

The name lookup of `fn pair_by_name()` [packages/treetime-graph/src/pair_by_name.rs#L27](../../packages/treetime-graph/src/pair_by_name.rs#L27), which every command uses to pair FASTA records with leaves, is case-sensitive. `"USA"` and `"usa"` are treated as different names with no standard convention for which is correct.

## Current behavior

`ancestral` and `homoplasy` fill a leaf without a matching record with unknown characters and warn for each such leaf; they stop when more than a third of the leaves have no record, unless `--ignore-missing-alns` is set (`fn complete_alignment_for_leaves()` [packages/treetime/src/ancestral/attach.rs#L18](../../packages/treetime/src/ancestral/attach.rs#L18)). `timetree`, `optimize` and `prune` stop at the first leaf without a record:

- Dense partitions: `Leaf sequence not found: '<name>'` [packages/treetime/src/partition/marginal/dense/partition.rs#L60](../../packages/treetime/src/partition/marginal/dense/partition.rs#L60)
- Sparse partitions: `Leaf sequence not found after alignment completion: '<name>'` [packages/treetime/src/partition/fitch/passes.rs#L68](../../packages/treetime/src/partition/fitch/passes.rs#L68)

These messages give no guidance on near-matches or potential causes (case, whitespace).

## Impact

Users with mismatched names must manually edit either the tree or the FASTA file to make names match exactly. The error message does not help identify which names are close matches or suggest corrections.

## Possible improvements

1. Build a name index upfront and report all mismatches at once (not one at a time)
2. Provide fuzzy matching suggestions for near-misses
3. Support explicit name mapping via config file (see [config-file-multi-partition proposal](../proposals/config-file-multi-partition.md))
4. Support case-insensitive matching as an option

## Related issues

- [M-dates-date-column-requires-exact-header-match.md](M-dates-date-column-requires-exact-header-match.md) - date column detection matches exact headers only
- [M-mugration-column-detection-no-positional-fallback.md](M-mugration-column-detection-no-positional-fallback.md) - mugration column detection errors instead of falling back to header position
- [Multi-segment genome input not wired](N-io-multi-segment-genome-input.md) - related input architecture gap
- [Optimize command accepts only a single alignment](N-optimize-multi-alignment-input.md) - related input limitation

## Related

- [N-io-missing-name-policy-differs-across-subsystems.md](N-io-missing-name-policy-differs-across-subsystems.md) - reconciliation logic duplicated across subsystems
- Design: [kb/proposals/input-name-matching-validation.md](../proposals/input-name-matching-validation.md) - ecosystem survey, four design axes, architecture analysis
- [Configuration file format for multi-partition analysis](../proposals/config-file-multi-partition.md) - could include explicit name mappings
