# Implement tree inference from alignment

`load_input_data()` in the timetree initialization path hits a `todo!()` when the `--tree` flag is omitted. The `--tree` CLI argument is declared as `Option<PathBuf>`, implying it is optional, but the code path for tree inference from alignment data is not implemented.

v0 infers a tree from the alignment using neighbor-joining or other methods when no tree is provided.

## Impact

Users who invoke `treetime timetree` without `--tree` (expecting tree inference from alignment) get a panic instead of a helpful error message or actual tree inference. This blocks a standard v0 workflow.

## Fix

Either implement tree inference from alignment (a large feature tracked in [unimplemented.md](../algo/unimplemented.md)) or replace the `todo!()` with a descriptive error: "Tree inference from alignment is not yet implemented. Provide a tree with --tree."

## Related issues

- Source: [H-timetree-tree-inference-unimplemented.md](../issues/H-timetree-tree-inference-unimplemented.md) -- delete after full resolution
- [unimplemented.md](../algo/unimplemented.md) -- tree inference algorithms listed as unimplemented
