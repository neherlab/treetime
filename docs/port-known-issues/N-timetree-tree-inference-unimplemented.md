# No tree inference from alignment

When no `--tree` is provided, the code panics with `todo!()`. v0 can infer a
tree from the alignment using neighbor-joining or other methods.
