# Dropping a deep Auspice tree overflows the stack

`json_read()` and `json_write()` in `packages/treetime-utils/src/io/json.rs` use deser-json, whose drivers keep the nesting on the heap, so reading and writing a deeply nested Auspice tree need no extra stack. Dropping the parsed `AuspiceTree` is recursive and has no such protection: `AuspiceTreeNode` has the compiler-generated drop and a derived `Clone`, both recursive through `children: Vec<AuspiceTreeNode>`, so each node drops its `children` vector on the current stack [packages/treetime-io/src/auspice_types.rs#L179-L190](../../packages/treetime-io/src/auspice_types.rs#L179-L190).

## Evidence

Caterpillar trees (each internal node has one leaf and one internal child), read with `json_read_file()` into `AuspiceTree` on a thread with a 2 MiB stack, the stack size of a `spawn_blocking` thread of the app server:

| Tree depth | Dev profile                    | Release profile                                 |
| ---------- | ------------------------------ | ----------------------------------------------- |
| 5,000      | parse and drop succeed         | parse and drop succeed                          |
| 10,000     | parse succeeds, drop overflows | parse and drop succeed                          |
| 20,000     | not measured                   | parse and drop succeed                          |
| 30,000     | not measured                   | parse succeeds, drop overflows                  |
| 50,000     | not measured                   | parse and compact write succeed, drop overflows |

The depths up to 30,000 were measured with the earlier serde_json reader, whose drop behavior is the same, because `Drop` of the tree types did not change. The process aborts with `thread has overflowed its stack`. The trees of the example datasets are much shallower: 157 for `sc2/4500`, 73 for `rsv/a/2000`, 33 for `mpox/clade-ii/500`.

## Other recursive operations

- `ResultTree::from_auspice()` and `result_colorings()` in the app back end walk the tree with the iterative `fn preorder()`, which keeps its stack on the heap [packages/app-backend/src/results/tree.rs#L111-L120](../../packages/app-backend/src/results/tree.rs#L111-L120), so they do not limit the depth

> [!IMPORTANT]
> **Investigation required.** The derived `Clone` of `AuspiceTreeNode` is recursive like the drop and is not measured. The depth of the deepest trees that users load (for example ladder-like trees of large outbreaks) is not known.
