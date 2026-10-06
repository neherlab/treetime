# Nexus parser reads trees only

> [!IMPORTANT]
> **Decision required.** The parser reads the `Trees` blocks of a Nexus file and skips every other block. Whether the tree-only scope is the complete Nexus contract of TreeTime is not approved.

`fn nexus_from_string()` in [packages/util-newick/src/nexus.rs](../../packages/util-newick/src/nexus.rs) parses the block and command structure of a Nexus file with a `pest` grammar that skips comments and quoted text. It reads `Translate` tables per `Trees` block and `Tree` and `UTree` commands. A block without its end, or a command without its `;`, is an error ([kb/decisions/multi-format-tree-io.md](../decisions/multi-format-tree-io.md)).

## Not read

- `Taxa` blocks: tree labels are not checked against `TaxLabels`
- `Characters`, `Data`, `Sets`, `Assumptions` and program blocks (FigTree, MrBayes, PAUP*): skipped command by command

## Decision axes

- **Tree-only scope**: keep the tree-only reader, or read the `Taxa` block and check the tree labels against it
- **Other blocks**: skip them, or read alignments from `Characters` and `Data` blocks as an alignment input
