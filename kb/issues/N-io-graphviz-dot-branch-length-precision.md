# Graphviz DOT writer uses Newick precision defaults for edge labels

> [!IMPORTANT]
> **Decision required.** The display precision of DOT edge labels (for example 6 significant digits, as under "Fix") needs approval independently from the Newick exchange-format default in [M-io-nwk-writer-3-sigfig-default-truncates-precision.md](M-io-nwk-writer-3-sigfig-default-truncates-precision.md).

The Graphviz DOT writer calls `format_weight` with `NwkWriteOptions::default()` for the `xlabel` of each edge (`fn print_edges()`, [packages/treetime-io/src/graphviz.rs#L108](../../packages/treetime-io/src/graphviz.rs#L108)). `fn format_weight()` falls back to 3 significant digits when the options give none ([packages/treetime-io/src/nwk.rs#L212](../../packages/treetime-io/src/nwk.rs#L212)). DOT output is for visualization, not data exchange -- full-precision labels like `0.00123456789012345` on every edge make graphs unreadable.

Currently this produces 3-significant-digit labels (matching the Newick default). After the Newick precision default is fixed ([M-io-nwk-writer-3-sigfig-default-truncates-precision.md](M-io-nwk-writer-3-sigfig-default-truncates-precision.md)), the DOT writer would produce full-precision labels unless it sets its own explicit precision.

## Fix

Use explicit precision for visualization output (e.g. 6 significant digits):

```rust
format_weight(weight, &NwkWriteOptions {
  weight_significant_digits: Some(6),
  ..NwkWriteOptions::default()
})
```

## Locations

- `packages/treetime-io/src/graphviz.rs` `fn print_edges()`
- `packages/treetime-io/src/nwk.rs` `fn format_weight()`

## Related

- Blocked by: [M-io-nwk-writer-3-sigfig-default-truncates-precision.md](M-io-nwk-writer-3-sigfig-default-truncates-precision.md)
- Design: [kb/proposals/newick-branch-length-precision.md](../proposals/newick-branch-length-precision.md) (Follow-up section)

Implementation waits until the display precision is approved independently from the exchange-format default.

## Proposed implementation plan

The display precision is open (Decision required above). The plan is the same for any approved value `N`.

1. In `fn print_edges()` (`packages/treetime-io/src/graphviz.rs`), pass `NwkWriteOptions { weight_significant_digits: Some(N), ..NwkWriteOptions::default() }` to `format_weight`, so that the DOT label stops depending on the Newick default. Keep the `weight` attribute at full precision, because Graphviz uses it for layout only
2. Land this before or together with the Newick default change of [M-io-nwk-writer-3-sigfig-default-truncates-precision.md](M-io-nwk-writer-3-sigfig-default-truncates-precision.md), so that DOT labels never switch to full precision
3. Test in `packages/treetime-io/src/__tests__/test_graphviz.rs`: an edge weight `0.00123456789012345` gives `xlabel="<value rounded to N significant digits>"` in the whole expected DOT text (oracle: decimal rounding of the input by hand)
4. Validate: `just lint-rs`, `just test-rs -E 'package(treetime-io)'`
