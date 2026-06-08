# Replace bio crate Newick parser with custom dialect-aware parser

Replace the `bio::io::newick::read` dependency with a custom Newick parser that handles all comment dialects: plain comments, BEAST-style `[&key=value]`, NHX `[&&NHX:key=value]`, Extended Newick `#H1` reticulate markers, and Rich Newick `[&R]`/`[&U]` rooting markers.

## Scope

- Remove `bio` crate dependency for Newick parsing (keep if used elsewhere)
- Implement a Newick parser that produces the same `Graph<N, E, D>` output as the current `nwk_read`
- Parse and preserve `[&...]` and `[&&NHX:...]` annotations into the `BTreeMap<String, String>` comments on `NodeFromNwk::from_nwk` (currently hardcoded to `btreemap! {}` at [packages/treetime-io/src/nwk.rs#L73](../../packages/treetime-io/src/nwk.rs#L73))
- Auto-detect dialect from content (not file extension): `[&&NHX` -> NHX mode, `[&` -> BEAST mode, `[...]` without `&` -> strip
- Handle all value types: scalars, quoted strings, `{...}` arrays, `#RRGGBB` colors
- Handle node vs branch annotation placement (BEAST2 canonical: `nodeMeta` before `:`, `lengthMeta` after `:` before length)
- Handle eNewick `#Type Integer` markers for hybrid nodes (future network support)

## Grammar specifications

All dialect grammars are documented in [kb/reports/newick-annotation-dialects.md](../reports/newick-annotation-dialects.md) with production rules, value type catalogs, escaping conventions, and source code references from 9 tool repositories.

## Implementation approach

BEAST2's ANTLR4 grammar at [NewickParser.g4](https://github.com/CompEvol/beast2/blob/9321c88df6e5e90da4db5c0fc4872576990c015c/src/beast/base/evolution/tree/treeparser/NewickParser.g4) is a formal grammar covering all BEAST-style features. It can serve as the reference for the parser structure. Consider `pest` or `nom` as Rust parser options.

## Related issues

- Source: [M-io-bio-crate-newick-rejects-comments.md](../issues/M-io-bio-crate-newick-rejects-comments.md)
- Source: [M-io-newick-output-incompatible-with-reader.md](../issues/M-io-newick-output-incompatible-with-reader.md)
