# Implement augur node data JSON serde types

Define typed Rust serde structs for producing augur-compatible node data JSON output. These types live in a new `packages/util-augur-node-data-json/` crate, following the `util-*` naming convention used by `util-phyloxml` and `util-usher-mat` for auxiliary format crates.

## Design

See [../reports/augur-node-data-json.md](../reports/augur-node-data-json.md) section "Rust type design" for the full design. Summary:

- Generic container `AugurNodeDataJson<M, N>` with `generated_by`, `nodes: BTreeMap<String, N>`, and `#[serde(flatten)] metadata: M`
- Per-command type aliases: `AugurNodeDataJsonAncestral`, `AugurNodeDataJsonRefine`, `AugurNodeDataJsonTraits`
- Per-command metadata structs: `AugurNodeDataJsonAncestralMeta`, `AugurNodeDataJsonRefineMeta`, `AugurNodeDataJsonTraitsMeta`
- Per-command node structs: `AugurNodeDataJsonAncestralNode`, `AugurNodeDataJsonRefineNode`, `AugurNodeDataJsonTraitsNode`
- All structs carry `#[serde(flatten)] pub other: BTreeMap<String, Value>` for forward compatibility
- `BTreeMap` key ordering gives sorted keys matching augur's `sort_keys=True`

AA fields (`aa_muts`, `aa_sequences`) should be present as `Option` fields with `#[serde(skip_serializing_if)]`, marked with a code comment linking to [../proposals/node-data-json-aa-reconstruction.md](../proposals/node-data-json-aa-reconstruction.md).

## Scope

Types and serialization only. No output wiring (that belongs in per-command tickets). Includes unit tests verifying serialization matches augur test fixture format.

## Reference

- Format specification: [../reports/augur-node-data-json.md](../reports/augur-node-data-json.md)
- Augur test fixtures for validation: [`tests/functional/ancestral/data/`](https://github.com/nextstrain/augur/tree/024292af6daf/tests/functional/ancestral/data), [`tests/functional/refine/data/`](https://github.com/nextstrain/augur/tree/024292af6daf/tests/functional/refine/data), [`tests/functional/traits/`](https://github.com/nextstrain/augur/tree/024292af6daf/tests/functional/traits)
