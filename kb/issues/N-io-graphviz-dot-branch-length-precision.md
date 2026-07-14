# Graphviz DOT writer uses Newick precision defaults for edge labels

The Graphviz DOT writer at `packages/treetime-io/src/graphviz.rs` calls `format_weight` with `NwkWriteOptions::default()`. DOT output is for visualization, not data exchange -- full-precision labels like `0.00123456789012345` on every edge make graphs unreadable.

Currently this produces 3-significant-digit labels (matching the Newick default). After the Newick precision default is fixed ([N-io-nwk-writer-3-sigfig-default-truncates-precision.md](N-io-nwk-writer-3-sigfig-default-truncates-precision.md)), the DOT writer would produce full-precision labels unless it sets its own explicit precision.

## Fix

Use explicit precision for visualization output (e.g. 6 significant digits):

```rust
format_weight(weight, &NwkWriteOptions {
  weight_significant_digits: Some(6),
  ..NwkWriteOptions::default()
})
```

## Locations

- `packages/treetime-io/src/graphviz.rs`

## Related

- Blocked by: [N-io-nwk-writer-3-sigfig-default-truncates-precision.md](N-io-nwk-writer-3-sigfig-default-truncates-precision.md)
- Design: [kb/proposals/newick-branch-length-precision.md](../proposals/newick-branch-length-precision.md) (Follow-up section)

No implementation ticket is ready until the display precision is approved independently from the exchange-format default.
