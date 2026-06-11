# Nexus writer produces double semicolon in TREES block

`nex_write_with` at `packages/treetime-io/src/nex.rs` embeds the Newick string into a template `Tree tree1={nwk};`. The Newick writer (`nwk_write_with`) appends its own trailing `;`, so the output is `Tree tree1=(...)root;;` with two semicolons.

Most Nexus parsers accept this (PAUP\*, Mesquite, FigTree all tolerate `;;`). Strictly, the NEXUS spec treats `;` as a command terminator, not a separator, so a bare `;` after the tree statement is a valid empty command.

Fix: either strip the trailing `;` from the Newick string before embedding, or change `nwk_write_with` to accept a "no trailing semicolon" option for embedded use.

## Locations

- Template: `packages/treetime-io/src/nex.rs` `nex_write_with()` `Tree tree1={nwk};` line
- Trailing semicolon: `packages/treetime-io/src/nwk.rs` `nwk_write_with()` `write!(writer, ";")`
