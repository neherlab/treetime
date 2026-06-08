# bio crate Newick parser rejects standard comments

The `bio` crate (v2.3.0) `newick::read` function fails on any `[...]` content in a Newick string. The Newick and NEXUS specifications define `[...]` as comments that parsers should strip and ignore, but `bio` treats `[` as an unexpected token.

This is stricter than the spec and prevents reading annotated trees from BEAST, MrBayes, FigTree, IQ-TREE, TreeTime v0, and any tool that writes `[&key=value]` or `[&&NHX:key=value]` annotations.

## Locations

- Reader entry point: `bio::io::newick::read` called from [packages/treetime-io/src/nwk.rs#L51](../../packages/treetime-io/src/nwk.rs#L51)
- Dependency: `bio = { version = "=2.3.0", features = ["phylogeny"] }` in workspace `Cargo.toml`
- Comment handling TODO at [packages/treetime-io/src/nwk.rs#L73](../../packages/treetime-io/src/nwk.rs#L73): `let comments = btreemap! {}; // TODO: parse nwk comments`

## Impact

- v1 cannot read annotated trees from any external tool
- v1 cannot read its own annotated output (see [M-io-newick-output-incompatible-with-reader.md](M-io-newick-output-incompatible-with-reader.md))
- v1 cannot read trees with plain comments `[this is a comment]` that some tools emit

## Resolution

Replace `bio::io::newick` with a custom parser that handles all Newick comment dialects. Grammar specifications for all dialects are documented in [kb/reports/newick-annotation-dialects.md](../reports/newick-annotation-dialects.md).
