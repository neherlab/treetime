# Newick writer defaults to 3 significant digits, truncating branch lengths

> [!IMPORTANT]
> **Decision required.** The default branch-length precision of the Newick exchange format is not approved. The tool survey and the options are in [kb/proposals/newick-branch-length-precision.md](../proposals/newick-branch-length-precision.md).

The Newick and Nexus writers format branch lengths through `fn newick_write_options()` at [packages/treetime-io/src/nwk.rs#L201-L210](../../packages/treetime-io/src/nwk.rs#L201-L210), and the Graphviz DOT writer through `fn format_weight()` at [packages/treetime-io/src/nwk.rs#L212-L218](../../packages/treetime-io/src/nwk.rs#L212-L218). Both default to 3 significant digits:

```rust
significant_digits: options.weight_significant_digits.or(Some(3)),
```

No production caller sets `weight_significant_digits`; only tests set it (to 17). So every Newick and Nexus tree output of every command writes branch lengths with 3 significant digits. v0 writes them through Biopython with 5 decimal places.

`0.123456` becomes `0.123`. This silently destroys precision on read-write-read cycles.

Major parsers write higher precision by default: Biopython uses 5 decimal places (`%1.5f`), ETE uses 6 significant digits (`%g`), DendroPy defaults to full precision via Python's `str(float)`, gotree writes full precision. See [kb/proposals/newick-branch-length-precision.md](../proposals/newick-branch-length-precision.md) for a full 10-tool source-verified survey.

## Impact

Branch lengths lose precision on every write cycle. Repeated read-write-read degrades values. Scientific reproducibility affected when trees are re-read after writing.

A public mugration report describes small branch lengths being written as zero in annotated Nexus output [[issue](https://github.com/neherlab/treetime/issues/271)]. The maintainer response attributes the observed output to branch collapsing or output precision and reports increasing the precision [[comment](https://github.com/neherlab/treetime/issues/271#issuecomment-2282346277)]. That report concerns the Python Nexus path; it is related evidence for preserving branch-length precision, not evidence that both implementations share a formatter or root cause.

## Locations

- `fn newick_write_options()` [packages/treetime-io/src/nwk.rs#L201-L210](../../packages/treetime-io/src/nwk.rs#L201-L210): used by the Newick writer and by the Nexus writer [packages/treetime-io/src/nex.rs#L34](../../packages/treetime-io/src/nex.rs#L34)
- `fn format_weight()` [packages/treetime-io/src/nwk.rs#L212-L218](../../packages/treetime-io/src/nwk.rs#L212-L218): used by the DOT writer with `NwkWriteOptions::default()` [packages/treetime-io/src/graphviz.rs#L108](../../packages/treetime-io/src/graphviz.rs#L108)
- `NwkWriteOptions.weight_significant_digits` defaults to `None`, and both sites replace `None` with `Some(3)`
- `fn float_to_digits()` in [packages/treetime-utils/src/fmt/float.rs#L41-L43](../../packages/treetime-utils/src/fmt/float.rs#L41-L43) sets the same 3-digit default again when no limit is given

Implementation waits until the exchange-format default is approved.

## Proposed implementation plan

The default precision is open (Decision required above). The plan holds for any approved default `D` (full round-trip precision, or a fixed count of significant or decimal digits).

1. One default in one place: set `D` as the `#[default]` of `NwkWriteOptions.weight_significant_digits` / `weight_decimal_digits` (`packages/treetime-io/src/nwk.rs`), and remove the `.or(Some(3))` fallbacks in `fn newick_write_options()` and `fn format_weight()`. For full precision, make the `None` path of `fn float_to_digits()` (`packages/treetime-utils/src/fmt/float.rs`) write the shortest round-trip representation instead of 3 digits, after checking its other callers
2. Land [N-io-graphviz-dot-branch-length-precision.md](N-io-graphviz-dot-branch-length-precision.md) first or in the same change, so that DOT labels keep a display precision
3. Tests in `packages/treetime-io/src/__tests__/`: write and read back a tree with branch lengths `0.123456789012345` and `1e-9`; with full precision the parsed values equal the input exactly, with a fixed `D` the written text equals the hand-rounded value (oracle: decimal rounding by hand, Biopython `%1.5f` for v0 parity if `D` is 5 decimals). Nexus writer: same case through `packages/treetime-io/src/nex.rs`
4. Validate: `just lint-rs`, `just test-rs`, and `just smoke`, which shows the expected change in every Newick and Nexus output; update the golden outputs only after checking that branch lengths are the only change
