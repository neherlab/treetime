# Newick writer does not quote names with special characters

`nwk_write_with` at `packages/treetime-io/src/nwk.rs#L242-L243` writes node names bare:

```rust
if let Some(name) = name {
  write!(writer, "{name}")?;
}
```

Names containing `(),:;[]' ` or whitespace produce invalid Newick that no parser can read back. The Newick standard requires quoting such names with single quotes (`'name with spaces'`), with internal quotes escaped as `''`.

All major parsers quote names on write: Biopython (`NewickIO.quote()`), ETE4 (`quote()`), DendroPy, gotree.

## Impact

Any node name containing special characters produces output that fails to roundtrip through the reader. Currently no bundled dataset triggers this because taxon names are simple identifiers, but external trees with spaces in names will break on re-read.

## Locations

- Writer: `packages/treetime-io/src/nwk.rs#L242-L243`
- Reader (handles quoting correctly): `packages/util-newick/src/parse.rs` `quoted_label` rule

## Related issues

- [kb/tickets/io-nwk-writer-style-dispatch.md](../tickets/io-nwk-writer-style-dispatch.md)
