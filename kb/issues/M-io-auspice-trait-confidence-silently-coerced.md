# Auspice trait fields are silently dropped when malformed

The Auspice reader treats malformed trait objects as absent and drops individual non-numeric confidence entries. A document can therefore return success with a partial probability map or without a trait that was present in the input.

> [!NOTE]
> The tree-output refactor removed the former `tree_ir` Auspice reader this issue originally cited. The current reader lives in [`packages/treetime-io/src/auspice.rs`](../../packages/treetime-io/src/auspice.rs). Whether it rejects or still silently drops malformed trait and confidence entries is **not yet confirmed**; the required contract below stands regardless.

## Evidence

- The former reader returned `Option` from trait parsing: a missing or non-string `value` yielded `None`, and `filter_map()` silently removed confidence members whose values were not JSON numbers. A trait was inserted only when parsing returned `Some`, so malformed present input became indistinguishable from absence.

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

