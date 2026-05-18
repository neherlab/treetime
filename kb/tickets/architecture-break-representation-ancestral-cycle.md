# Break representation/ -> ancestral/ dependency

## Description

Two imports in `representation/` reach up into `ancestral/` (layer 3 importing from layer 4):

- `partition/fitch.rs:2:` -> `ancestral::fitch::compress_sequences`
- `partition/marginal_dense.rs:4:` -> `ancestral::fitch_indel::{resolve_indels_backward, resolve_indels_forward}`

## Fix

### resolve_indels_backward/forward

Pure functions on gap ranges and deletion maps. Input types are `Vec<(usize, usize)>`, `BTreeMap<(usize, usize), Deletion>`, `usize`. No dependency on ancestral-specific types. Move to `seq/indel.rs` (alongside existing `InDel` type) or `partition/indel_resolution.rs`.

### compress_sequences

Called by `PartitionFitch::compress()` convenience constructor. Two options:

- Remove the convenience constructor. Callers already can call `ancestral::fitch::compress_sequences` directly with a `PartitionFitch` reference.
- Move the constructor call to a free function in `ancestral/` that creates and returns a `PartitionFitch`. Inverts the dependency direction.

## Related issues

- Source: [M-core-remaining-architectural-debt-after-extraction](../issues/M-core-remaining-architectural-debt-after-extraction.md)
