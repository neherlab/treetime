# Newick writer defaults to 3 significant digits, truncating branch lengths

`format_weight` at `packages/treetime-io/src/nwk.rs#L274` defaults to 3 significant digits:

```rust
float_to_digits(
  weight,
  options.weight_significant_digits.or(Some(3)),
  options.weight_decimal_digits,
)
```

`0.123456` becomes `0.123`. This silently destroys precision on read-write-read cycles.

Major parsers write higher precision by default: Biopython uses 5 decimal places (`%1.5f`), ETE uses 6 significant digits (`%g`), DendroPy uses 10 digits (`%.10e`), gotree writes full precision.

## Impact

Branch lengths lose precision on every write cycle. Repeated read-write-read degrades values. Scientific reproducibility affected when trees are re-read after writing.

## Locations

- `packages/treetime-io/src/nwk.rs#L268-L277`
- `NwkWriteOptions.weight_significant_digits` defaults to `None` but `format_weight` overrides with `Some(3)`
