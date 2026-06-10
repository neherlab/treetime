# Add NWK writer annotation style dispatch

Add `NwkStyle` enum (Plain, Annotated, Nhx) and replace the hardcoded BEAST format string at `treetime-io/src/nwk.rs:255` with a dispatch on style. Default to Plain, fixing the round-trip bug where annotated output is rejected by the reader.

## Done (`util-newick` crate)

The `util-newick` crate (`packages/util-newick/`) implements the standalone writer:

- `NwkStyle` enum: `Plain`, `Beast`, `Nhx` with configurable float precision
- BEAST2 canonical annotation placement (node attrs before `:`, branch attrs after `:` before length)
- BEAST key/value escaping (boolean/numeric lookalikes quoted, embedded quotes escaped)
- NHX serialization with reserved-char rejection (fallible `newick_to_string` returning `Result`)
- Nexus streaming writer with Taxa block, Translate table, multi-tree support

## Remaining (integration into `treetime-io`)

- Add `NwkStyle` parameter to `nwk_write_with`, `nwk_write_file_with`, `nex_write_with`, `nex_write_file_with`
- Replace the format string at `nwk.rs:255` with dispatch to `util-newick` writer
- Default: `Plain` -- suppresses all annotations, fixing the round-trip incompatibility
- Update all call sites in command `run.rs` files to pass `NwkStyle::Plain`

## Original scope

- Define `NwkStyle` enum: `Plain` (no annotations), `Annotated` (BEAST `[&key="value"]`), `Nhx` (NHX `[&&NHX:key=value:...]`)
- Add `NwkStyle` parameter to `nwk_write_with`, `nwk_write_file_with`, `nex_write_with`, `nex_write_file_with`
- Replace the format string at `nwk.rs:255` (`format!("[&{key}=\"{val}\"]")`) with a match on `NwkStyle`
- Implement NHX serialization (colon-separated, `&&NHX` prefix)
- Default: `Plain` -- suppresses all annotations, fixing the round-trip incompatibility (issue #8 in tickets3)
- Update all call sites in command `run.rs` files to pass `NwkStyle::Plain` (preserving current effective behavior since the reader can't parse annotations anyway)

## Locations

- Writer: `packages/treetime-io/src/nwk.rs:251-259` (comment serialization)
- Comment merge: `packages/treetime-io/src/nwk.rs:239` (two-source merge into BTreeMap)
- Shared output: `packages/treetime-io/src/graph.rs:18-53` (`write_graph_files_with_options`)
- Command call sites: `ancestral/run.rs`, `timetree/run.rs`, `optimize/run.rs`, `mugration/run.rs`, `clock/run.rs`, `prune/run.rs`

## Tests

### Happy paths

- Plain: `BTreeMap` with entries -> no `[...]` in output
- Annotated single key: `{"mutations": "A55G"}` -> `[&mutations="A55G"]`
- Annotated multiple keys: `{"mutations": "A55G", "date": "2020.5"}` -> `[&mutations="A55G",date="2020.5"]`
- Nhx single key: `{"S": "human"}` -> `[&&NHX:S=human]`
- Nhx multiple keys: `{"S": "human", "D": "Y"}` -> `[&&NHX:S=human:D=Y]`
- Empty CommentProviders: no annotations regardless of style
- Both sources merged: `NodeToNwk::nwk_comments()` returns `{"date": "2020.5"}`, `CommentProviders` returns `{"mutations": "A55G"}` -> both appear in output
- Nexus output: embedded Newick string uses same style as standalone NWK

### Edge cases

- Value containing double quotes: must be escaped or handled in BEAST style
- Value containing commas: must not break BEAST key-value parsing (`[&label="A,B"]`)
- Value containing colons: must not break NHX key-value parsing
- Value containing `]`: must not prematurely close the comment
- Value containing `[`: must not open a nested comment
- Empty string value: `{"key": ""}` -> `[&key=""]` or `[&&NHX:key=]`
- Key containing dots: `{"posterior.prob": "0.95"}` -> valid in both dialects
- Unicode in values: `{"location": "Zurich"}` (umlaut, CJK, emoji)
- Single-entry BTreeMap
- Large BTreeMap (50+ entries)
- Whitespace in values: `{"label": "New York"}`

### Pathological

- Value containing newlines
- Value containing NUL bytes
- Empty key in BTreeMap
- Thousands of annotations on a single node (performance)

### Regression

- All existing NWK/Nexus output tests pass with Plain default
- Round-trip: write with Annotated, read back with current parser (fails until ticket #3 ships -- documenting expected failure)

## Related issues

- Source: resolved `M-io-newick-output-incompatible-with-reader` (reader no longer rejects annotated output; writer format issues remain)
- Source: [kb/issues/N-io-nexus-writer-no-annotation-advantage.md](../issues/N-io-nexus-writer-no-annotation-advantage.md)
- Source: [kb/issues/M-io-nwk-writer-annotation-format-nonstandard.md](../issues/M-io-nwk-writer-annotation-format-nonstandard.md)
- Source: [kb/issues/M-io-nwk-writer-annotations-after-branch-length.md](../issues/M-io-nwk-writer-annotations-after-branch-length.md)
- Source: [kb/issues/M-io-nwk-writer-no-quoting-special-chars.md](../issues/M-io-nwk-writer-no-quoting-special-chars.md)
