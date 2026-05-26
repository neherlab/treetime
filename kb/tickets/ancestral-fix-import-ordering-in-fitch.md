# Fix import ordering in `fitch.rs`

`use std::collections::BTreeSet;` at `packages/treetime/src/ancestral/fitch.rs:2:` is placed between two `use crate::` imports. Move after the `crate::` block to restore standard import grouping, or run `rustfmt` with import grouping enabled.

## Related issues

Source: `.memory/15-pr-714/synthesis1.md` finding L1
