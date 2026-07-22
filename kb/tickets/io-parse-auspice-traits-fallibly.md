# Parse Auspice traits fallibly

Reject malformed present Auspice traits instead of dropping the trait or individual confidence entries.

> [!NOTE]
> This ticket predates the tree-output refactor, which removed the `tree_ir` layer and now writes formats directly from the graph in `packages/treetime/src/commands/shared/tree_output.rs` (readers in `packages/treetime-io/src/auspice.rs` and `packages/treetime-io/src/phyloxml.rs`). Its parent issue's current status is unconfirmed; re-derive the steps and code locations against the current code before executing.

## Required changes

- Change `parse_trait()` [packages/treetime-io/src/tree_ir/auspice.rs#L324](../../packages/treetime-io/src/tree_ir/auspice.rs#L324) to return `Result<TreeIrTrait, Report>` for a present trait object.
- Require a string `value` and validate the types of present `confidence` and `entropy` fields.
- Collect every confidence entry fallibly; reject the complete trait on the first non-numeric or non-finite member.
- Reject non-numeric or non-finite entropy.
- Propagate errors from `auspice_node_to_graph_components()` [packages/treetime-io/src/tree_ir/auspice.rs#L218](../../packages/treetime-io/src/tree_ir/auspice.rs#L218) with the node name, trait attribute, field path, confidence key where applicable, and original value.

## Validation

- Complete traits, absent traits, and present traits with missing or malformed values.
- Valid, empty, partially malformed, and wrong-container confidence cases.
- Valid and malformed entropy cases.
- A whole-document test proving a present malformed trait fails instead of disappearing.
- Whole-value round trips for valid documents.
- Full lint and test suite.

## Related issues

- Source: [kb/issues/M-io-auspice-trait-confidence-silently-coerced.md](../issues/M-io-auspice-trait-confidence-silently-coerced.md)
