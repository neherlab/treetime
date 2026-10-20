# GTR API casing is inconsistent

The same domain uses the all-caps type name `GTR` [`packages/treetime/src/gtr/gtr.rs#L12`](../../packages/treetime/src/gtr/gtr.rs#L12) alongside Rust-style `GtrOutput`, `GtrModelType`, and `GtrModelName` [`packages/treetime/src/gtr/get_gtr.rs#L79`](../../packages/treetime/src/gtr/get_gtr.rs#L79) [`packages/treetime/src/gtr/get_gtr.rs#L130`](../../packages/treetime/src/gtr/get_gtr.rs#L130) [`packages/treetime/src/gtr/get_gtr.rs#L154`](../../packages/treetime/src/gtr/get_gtr.rs#L154). Model methods likewise mix `expQt_with_rate` and `expQt` [`packages/treetime/src/gtr/gtr.rs#L147`](../../packages/treetime/src/gtr/gtr.rs#L147) [`packages/treetime/src/gtr/gtr.rs#L184`](../../packages/treetime/src/gtr/gtr.rs#L184) with snake-case methods such as `propagate_profile` [`packages/treetime/src/gtr/gtr.rs#L169`](../../packages/treetime/src/gtr/gtr.rs#L169). The model parameter structs mix all-caps acronyms (`JC69Params`, `HKY85Params`, `TN93Params`) with `Jtt92Params` [`packages/treetime/src/gtr/get_gtr.rs#L187`](../../packages/treetime/src/gtr/get_gtr.rs#L187) [`packages/treetime/src/gtr/get_gtr.rs#L266`](../../packages/treetime/src/gtr/get_gtr.rs#L266) [`packages/treetime/src/gtr/get_gtr.rs#L375`](../../packages/treetime/src/gtr/get_gtr.rs#L375) [`packages/treetime/src/gtr/get_gtr.rs#L425`](../../packages/treetime/src/gtr/get_gtr.rs#L425).

## Impact

The API does not provide one predictable transformation from domain notation to Rust identifiers. Search must account for `GTR`, `Gtr`, `Qt`, and single-capital method forms, and new code has no unambiguous precedent.

## Required convention

Use Rust casing for public identifiers: `Gtr...` for types and snake case for functions and methods. Preserve mathematical symbols such as $Q$ and $e^{Qt}$ in documentation and equations. Rename the complete public family together so the cleanup does not introduce another mixed state.

## Validation

- Semantic rename covers imports, trait methods, generated documentation, tests, and CLI/schema references.
- Search confirms the old identifier family is absent outside citations to v0 or external formulas.
- Numerical GTR tests prove the change is behavior-preserving.
