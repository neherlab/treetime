# Auspice trait reader contract is unverified after the tree I/O refactor

The app back end reads Auspice JSON into the typed structs of [packages/treetime-io/src/auspice_types.rs](../../packages/treetime-io/src/auspice_types.rs). It reads only the Auspice output of a TreeTime run: `fn read_auspice()` and `fn RunOutputs::read_file()` load the file that the run wrote [packages/app-backend/src/results/outputs.rs#L20-L24](../../packages/app-backend/src/results/outputs.rs#L20-L24), [packages/app-backend/src/results/outputs.rs#L60-L65](../../packages/app-backend/src/results/outputs.rs#L60-L65). The consumers are `fn ResultTree::from_auspice()` and `fn result_colorings()` in [packages/app-backend/src/results/tree.rs](../../packages/app-backend/src/results/tree.rs) and the mugration summary in [packages/app-backend/src/results/mugration.rs](../../packages/app-backend/src/results/mugration.rs).

## Current behavior

- `fn AuspiceTreeNodeAttrs::attr()` deserializes a named node attribute and returns an error with the attribute key when the value does not fit [packages/treetime-io/src/auspice_types.rs#L254-L271](../../packages/treetime-io/src/auspice_types.rs#L254-L271). An absent attribute gives `None`
- `struct AuspiceTreeNodeAttr` holds `value` as a `String` [packages/treetime-io/src/auspice_types.rs#L274-L280](../../packages/treetime-io/src/auspice_types.rs#L274-L280). A numeric `value`, as Auspice uses for continuous attributes, fails with an error
- `fn AuspiceTreeNodeAttr::state_confidence()` reads `confidence` only when it is a map of state probabilities [packages/treetime-io/src/auspice_types.rs#L294-L301](../../packages/treetime-io/src/auspice_types.rs#L294-L301). A map with non-numeric entries fails with an error. A `confidence` of any other JSON type gives `Ok(None)` without a message, so the mugration summary omits the ancestor probabilities of that node [packages/app-backend/src/results/mugration.rs#L55-L57](../../packages/app-backend/src/results/mugration.rs#L55-L57)
- `entropy` is not read. It stays in the untyped `other` map

TreeTime writes `value` as a string and `confidence` as a map, so its own outputs do not reach the silent path.

## Remaining problem

The reader has no stated contract for Auspice JSON that TreeTime did not write: a non-map `confidence` is dropped without a message, and a numeric `value` is an error.

> [!IMPORTANT]
> **Decision required.** It is not decided whether the app reads Auspice JSON from other tools. Today it reads only the output of its own runs. The options are:
>
> - Keep the reader limited to TreeTime outputs: return an error for a `confidence` that is not a map, and document that the reader accepts only the shapes TreeTime writes
> - Accept foreign Auspice JSON: take the accepted shapes of `value`, `confidence` and `entropy` from the Auspice v2 schema, which is not checked yet, and type the reader for them

## Validation

- Whole-document tests of the reader with a non-map `confidence` and a numeric `value`, with the outcome that the decision selects
