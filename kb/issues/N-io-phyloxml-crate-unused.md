# PhyloXML crate is not used by any command

`packages/util-phyloxml` has a PhyloXML reader and writer (`fn phyloxml_read()` and `fn phyloxml_write()` in [packages/util-phyloxml/src/lib.rs](../../packages/util-phyloxml/src/lib.rs)), but no package depends on it. The graph converter was removed with other unused read paths in commit `0460901b`. The pathogen genomics tools that TreeTime works with (augur, Nextclade, UShER, Taxonium, IQ-TREE 3) have no PhyloXML support in their sources.

> [!IMPORTANT]
> **Decision required.** The crate is a workspace member ([Cargo.toml#L47](../../Cargo.toml#L47)) without a consumer. The project rules require the removal of unused code, while [kb/decisions/multi-format-tree-io.md](../decisions/multi-format-tree-io.md) describes the PhyloXML parser as a library crate that no command uses yet. The options are:
>
> - Wire PhyloXML into the tree output and input paths
> - Remove the crate

## Requirements if wired

- Parse the XML Schema boolean forms `true`, `false`, `1` and `0` for the TreeTime properties `bad_branch` and `date_inferred`. An earlier reader treated every value other than `true` as false
- Reject invalid values of recognized properties, with the property `ref` in the error message

## Related

- [kb/proposals/io-format-coverage.md](../proposals/io-format-coverage.md): ranking of format work
- [kb/decisions/multi-format-tree-io.md](../decisions/multi-format-tree-io.md): PhyloXML background

## Proposed implementation plan

Wire or remove is open (Decision required above). One short plan per option.

- Option "remove the crate": delete `packages/util-phyloxml`, its member line in `Cargo.toml`, and its workspace dependencies that no other crate uses; edit the PhyloXML paragraph of [kb/decisions/multi-format-tree-io.md](../decisions/multi-format-tree-io.md) and the row of [kb/proposals/io-format-coverage.md](../proposals/io-format-coverage.md); delete this issue. Validate: `just lint-rs`, `just test-rs`
- Option "wire PhyloXML":
  1. Add a converter between `Phyloxml` and the graph, names, and branch lengths in `packages/treetime-io/src/` beside `nwk.rs`, and select it by the `.xml`/`.phyloxml` extension in the tree input and tree output paths (`packages/app-output/src/tree_output.rs`)
  2. Apply the two requirements above to the property reader
  3. Tests: read a PhyloXML literal with `bad_branch` set to each of `true`, `false`, `1`, `0` and assert the whole parsed result; an invalid value gives an error that names the property `ref` (oracle: XML Schema `xs:boolean`); write and read back a tree and compare with the Newick result of the same tree
  4. Validate: `just lint-rs`, `just test-rs`, `just generated-check`
