# GTR API casing is inconsistent

The same domain uses all-caps type names such as `GTR`, `GTRParams`, and `GTRSiteSpecific` alongside Rust-style `GtrModelName` and `GtrOutput`. Model methods likewise mix `expQt`, `expQt_raw`, and `Q` with snake-case methods such as `average_rate` [`packages/treetime/src/gtr/gtr.rs#L109`](../../packages/treetime/src/gtr/gtr.rs#L109) [`packages/treetime/src/gtr/get_gtr.rs#L26`](../../packages/treetime/src/gtr/get_gtr.rs#L26) [`packages/treetime/src/gtr/gtr_site_specific.rs#L49`](../../packages/treetime/src/gtr/gtr_site_specific.rs#L49).

## Impact

The API does not provide one predictable transformation from domain notation to Rust identifiers. Search must account for `GTR`, `Gtr`, `Qt`, and single-capital method forms, and new code has no unambiguous precedent.

## Required convention

Use Rust casing for public identifiers: `Gtr...` for types and snake case for functions and methods. Preserve mathematical symbols such as $Q$ and $e^{Qt}$ in documentation and equations. Rename the complete public family together so the cleanup does not introduce another mixed state.

## Validation

- Semantic rename covers imports, trait methods, generated documentation, tests, and CLI/schema references.
- Search confirms the old identifier family is absent outside citations to v0 or external formulas.
- Numerical GTR tests prove the change is behavior-preserving.
