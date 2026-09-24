# `parse_delimited_str` drops a trailing empty field but keeps a leading one

## Summary

`parse_delimited_str` in `packages/treetime-io/src/parse_delimited.rs` treats empty fields differently depending on their position. The unit tests in `packages/treetime-io/src/__tests__/test_parse_delimited.rs` fix this behavior in place without stating why:

- `"A,B,"` parses to `["A", "B"]`: the empty field after the trailing delimiter is dropped
- `",A,B"` parses to `["", "A", "B"]`: the empty field before the leading delimiter is kept
- `",,"` parses to `["", ""]`: two fields, although the input has three field positions

## Consequence

The `prune` command in the CLI, server, and NAPI parses its `--prune-nodes-list` value with this function. A list with a trailing delimiter adds only the named nodes to the prune set, while the same list with a leading delimiter also adds the empty name `""`. How `prune` treats a name that matches no node has not been checked.

## Open question

Choose one rule for empty fields and apply it at both ends:

- **Tolerate edge delimiters**: drop empty fields at the start and the end, keep interior empty fields
- **Reject empty fields**: return an error that names the position of the empty field
- **Keep all fields**: return every field, including empty ones at both ends, and let callers validate
