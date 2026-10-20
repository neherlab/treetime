# Newick writer defaults to 3 significant digits, truncating branch lengths

> [!IMPORTANT]
> **Decision required.** The default branch-length precision of the Newick exchange format is not approved. The tool survey and the options are in [kb/proposals/newick-branch-length-precision.md](../proposals/newick-branch-length-precision.md).

`fn format_weight()` at [packages/treetime-io/src/nwk.rs#L199-L205](../../packages/treetime-io/src/nwk.rs#L199-L205) defaults to 3 significant digits:

```rust
float_to_digits(
  weight,
  options.weight_significant_digits.or(Some(3)),
  options.weight_decimal_digits,
)
```

`0.123456` becomes `0.123`. This silently destroys precision on read-write-read cycles.

Major parsers write higher precision by default: Biopython uses 5 decimal places (`%1.5f`), ETE uses 6 significant digits (`%g`), DendroPy defaults to full precision via Python's `str(float)`, gotree writes full precision. See [kb/proposals/newick-branch-length-precision.md](../proposals/newick-branch-length-precision.md) for a full 10-tool source-verified survey.

## Impact

Branch lengths lose precision on every write cycle. Repeated read-write-read degrades values. Scientific reproducibility affected when trees are re-read after writing.

A public mugration report describes small branch lengths being written as zero in annotated Nexus output [[issue](https://github.com/neherlab/treetime/issues/271)]. The maintainer response attributes the observed output to branch collapsing or output precision and reports increasing the precision [[comment](https://github.com/neherlab/treetime/issues/271#issuecomment-2282346277)]. That report concerns the Python Nexus path; it is related evidence for preserving branch-length precision, not evidence that both implementations share a formatter or root cause.

## Locations

- [packages/treetime-io/src/nwk.rs#L199-L205](../../packages/treetime-io/src/nwk.rs#L199-L205)
- `NwkWriteOptions.weight_significant_digits` defaults to `None` but `format_weight` overrides with `Some(3)`
- `fn float_to_digits()` in [packages/treetime-utils/src/fmt/float.rs#L45-L47](../../packages/treetime-utils/src/fmt/float.rs#L45-L47) sets the same 3-digit default again when no limit is given

Implementation waits until the exchange-format default is approved.
