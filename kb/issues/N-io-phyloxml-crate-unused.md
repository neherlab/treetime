# PhyloXML crate is not used by any command

`packages/util-phyloxml` has a PhyloXML reader and writer (`fn phyloxml_read()` and `fn phyloxml_write()` in [packages/util-phyloxml/src/lib.rs](../../packages/util-phyloxml/src/lib.rs)), but no package depends on it. The graph converter was removed with other unused read paths in commit `0460901b`. The pathogen genomics tools that TreeTime works with (augur, Nextclade, UShER, Taxonium, IQ-TREE 3) have no PhyloXML support in their sources.

> [!IMPORTANT]
> **Decision required.** The crate is a workspace member ([Cargo.toml#L44](../../Cargo.toml#L44)) without a consumer. The project rules require the removal of unused code, while [kb/decisions/multi-format-tree-io.md](../decisions/multi-format-tree-io.md) describes the PhyloXML parser as a library crate that no command uses yet. The options are:
>
> - Wire PhyloXML into the tree output and input paths
> - Remove the crate

## Requirements if wired

- Parse the XML Schema boolean forms `true`, `false`, `1` and `0` for the TreeTime properties `bad_branch` and `date_inferred`. An earlier reader treated every value other than `true` as false
- Reject invalid values of recognized properties, with the property `ref` in the error message

## Related

- [kb/proposals/io-format-coverage.md](../proposals/io-format-coverage.md): ranking of format work
- [kb/decisions/multi-format-tree-io.md](../decisions/multi-format-tree-io.md): PhyloXML background
