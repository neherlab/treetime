# resolve_tip_keys error and ambiguity paths untested

`packages/treetime/src/optimize/pipeline.rs:255` `resolve_tip_keys` is only tested with valid leaf names. Missing coverage for: missing tip name error, empty tip list, duplicate names in the tree, internal-node name matches. Also duplicated with `packages/treetime/src/clock/reroot.rs` (`find_node_key_by_name`), tracked in `kb/issues/N-reroot-duplicated-tip-name-resolution.md`.
