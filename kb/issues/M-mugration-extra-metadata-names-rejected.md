# Mugration rejects metadata entries not present in tree

Mugration's `validate_trait_names()` in [`packages/treetime/src/partition/marginal_discrete.rs#L251-L284`](../../packages/treetime/src/partition/marginal_discrete.rs#L251) hard-errors in both directions: tree leaves missing from metadata, and metadata names missing from tree. The second check causes failures on datasets where the metadata file contains samples that were pruned from the tree, which is common in practice.

## Cross-command inconsistency

v1 handles extra-metadata and missing-data scenarios differently across subsystems. The mugration validator is the strictest, while FASTA and dates follow v0's lenient approach.

| Subsystem                                                                                                       | Extra external names (not in tree) | Tree leaves without external data                    |
| --------------------------------------------------------------------------------------------------------------- | ---------------------------------- | ---------------------------------------------------- |
| FASTA ([`attach.rs#L21-L22`](../../packages/treetime/src/ancestral/attach.rs#L21))                              | Silently ignored                   | Warning + ambiguous fill; error if >1/3 missing      |
| Dates ([`date_constraints.rs#L92-L113`](../../packages/treetime/src/clock/date_constraints.rs#L92))             | Warning only                       | Marked as bad branch; error if insufficient coverage |
| Mugration ([`marginal_discrete.rs#L275-L281`](../../packages/treetime/src/partition/marginal_discrete.rs#L275)) | **Hard error**                     | **Hard error**                                       |

## v0 behavior

v0 mugration at [`wrappers.py#L777-L782`](../../packages/legacy/treetime/treetime/wrappers.py#L777) iterates tree terminals and builds pseudo-sequences by looking up each leaf name in the traits dict. Leaves without a trait entry get `missing_char`. Extra metadata entries are never accessed. Output: "Assigned discrete traits to 18 out of 18 taxa."

The tb/20 dataset has 20 metadata rows but only 18 tree leaves (G22694 and G22721 were pruned upstream). v0 succeeds; v1 errors.

## Correct behavior

Extra metadata entries should be a warning (or silent), matching FASTA and dates behavior. Metadata files commonly contain more samples than the tree.

Tree leaves missing from metadata should assign the missing-data character and warn, matching both v0 and the FASTA subsystem's approach. Whether to error when coverage is very low (as FASTA does at >1/3 missing) is a separate design question.

## Related issues

- [N-io-name-reconciliation-duplicated.md](N-io-name-reconciliation-duplicated.md) -- reconciliation logic duplicated across subsystems with inconsistent policies
- [M-io-sequence-name-matching-unreliable.md](M-io-sequence-name-matching-unreliable.md) -- broader name matching reliability
