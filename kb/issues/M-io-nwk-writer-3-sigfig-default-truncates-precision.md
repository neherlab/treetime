# Newick writer defaults to 3 significant digits, truncating branch lengths

> [!IMPORTANT]
> **Decision required.** The default branch-length precision of the Newick exchange format is not approved. The tool survey and the options are in [kb/proposals/newick-branch-length-precision.md](../proposals/newick-branch-length-precision.md).

The Newick and Nexus writers format branch lengths through `fn newick_write_options()` at [packages/treetime-io/src/nwk.rs#L222-L231](../../packages/treetime-io/src/nwk.rs#L222-L231), and the Graphviz DOT writer through `fn format_weight()` at [packages/treetime-io/src/nwk.rs#L233-L239](../../packages/treetime-io/src/nwk.rs#L233-L239). Both default to 3 significant digits:

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

- `fn newick_write_options()` [packages/treetime-io/src/nwk.rs#L222-L231](../../packages/treetime-io/src/nwk.rs#L222-L231): used by the Newick writer and by the Nexus writer [packages/treetime-io/src/nex.rs#L35](../../packages/treetime-io/src/nex.rs#L35)
- `fn format_weight()` [packages/treetime-io/src/nwk.rs#L233-L239](../../packages/treetime-io/src/nwk.rs#L233-L239): used by the DOT writer with `NwkWriteOptions::default()` [packages/treetime-io/src/graphviz.rs#L114](../../packages/treetime-io/src/graphviz.rs#L114)
- `NwkWriteOptions.weight_significant_digits` defaults to `None`, and both sites replace `None` with `Some(3)`
- `fn float_to_digits()` in [packages/treetime-utils/src/fmt/float.rs#L40-L42](../../packages/treetime-utils/src/fmt/float.rs#L40-L42) sets the same 3-digit default again when no limit is given

Implementation waits until the exchange-format default is approved.
