# Auspice trait fields are silently dropped when malformed

The Auspice reader treats malformed trait objects as absent and drops individual non-numeric confidence entries. A document can therefore return success with a partial probability map or without a trait that was present in the input.

## Evidence

- `fn parse_trait()` [packages/treetime-io/src/tree_ir/auspice.rs#L324](../../packages/treetime-io/src/tree_ir/auspice.rs#L324) returns `Option<TreeIrTrait>`. A missing or non-string `value` returns `None`, and `filter_map()` silently removes confidence members whose values are not JSON numbers.
- `fn AuspiceReader::auspice_node_to_graph_components()` [packages/treetime-io/src/tree_ir/auspice.rs#L218](../../packages/treetime-io/src/tree_ir/auspice.rs#L218) inserts a trait only when `parse_trait()` returns `Some`, so malformed present input becomes indistinguishable from absence.

## Required contract

- Absence of a trait attribute remains valid.
- Once a trait object is present, parsing returns `Result<TreeIrTrait, Report>` and requires a string `value`.
- If `confidence` is present, it must be an object and every member must be a numeric finite value; one malformed member rejects the trait and document.
- If `entropy` is present, it must be a numeric finite value.
- Errors identify the node, trait attribute, and failing field or confidence key.

## Options

- **Reject malformed present traits:** return an error with field context.
- **Preserve malformed raw traits:** retain them beside the typed representation. This keeps source data but permits an invalid scientific field to travel through typed TreeIR.

## Recommendation

Reject malformed present traits. Presence is an explicit request to parse typed scientific metadata; partial confidence maps are invalid states.

## Validation

- Complete traits with and without confidence and entropy.
- Missing, non-string, and valid string `value` cases.
- Confidence objects containing valid numbers, one non-number, and the empty map.
- Non-object confidence and non-number entropy.
- A whole-document test proving malformed present traits return an error rather than disappearing.

## Related issues

- [N-io-tree-ir-architecture-unapproved.md](N-io-tree-ir-architecture-unapproved.md)
