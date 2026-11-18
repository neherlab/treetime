# GTR API casing is inconsistent

The same domain uses the all-caps type name `GTR` [`packages/treetime/src/gtr/gtr.rs#L12`](../../packages/treetime/src/gtr/gtr.rs#L12) alongside Rust-style `GtrOutput`, `GtrModelType`, and `GtrModelName` [`packages/treetime/src/gtr/get_gtr.rs#L140`](../../packages/treetime/src/gtr/get_gtr.rs#L140) [`packages/treetime/src/gtr/get_gtr.rs#L199`](../../packages/treetime/src/gtr/get_gtr.rs#L199) [`packages/treetime/src/gtr/get_gtr.rs#L297`](../../packages/treetime/src/gtr/get_gtr.rs#L297). Model methods likewise mix `expQt_with_rate` and `expQt` [`packages/treetime/src/gtr/gtr.rs#L121`](../../packages/treetime/src/gtr/gtr.rs#L121) [`packages/treetime/src/gtr/gtr.rs#L158`](../../packages/treetime/src/gtr/gtr.rs#L158) with snake-case methods such as `propagate_profile` [`packages/treetime/src/gtr/gtr.rs#L143`](../../packages/treetime/src/gtr/gtr.rs#L143). The model parameter structs mix all-caps acronyms (`JC69Params`, `HKY85Params`, `TN93Params`) with `Jtt92Params` [`packages/treetime/src/gtr/get_gtr.rs#L375`](../../packages/treetime/src/gtr/get_gtr.rs#L375) [`packages/treetime/src/gtr/get_gtr.rs#L454`](../../packages/treetime/src/gtr/get_gtr.rs#L454) [`packages/treetime/src/gtr/get_gtr.rs#L613`](../../packages/treetime/src/gtr/get_gtr.rs#L613) [`packages/treetime/src/gtr/get_gtr.rs#L563`](../../packages/treetime/src/gtr/get_gtr.rs#L563).

## Impact

The API does not provide one predictable transformation from domain notation to Rust identifiers. Search must account for `GTR`, `Gtr`, `Qt`, and single-capital method forms, and new code has no unambiguous precedent.

## Required convention

Use Rust casing for public identifiers: `Gtr...` for types and snake case for functions and methods. Preserve mathematical symbols such as $Q$ and $e^{Qt}$ in documentation and equations. Rename the complete public family together so the cleanup does not introduce another mixed state.

## Validation

- Semantic rename covers imports, trait methods, generated documentation, tests, and CLI/schema references.
- Search confirms the old identifier family is absent outside citations to v0 or external formulas.
- Numerical GTR tests prove the change is behavior-preserving.

## Proposed implementation plan

All cited identifiers exist at the cited lines. `expQt` and `expQt_with_rate` are `pub(crate)`, the other identifiers are `pub`. `struct GTR` appears in 82 Rust files; the `*Params` family in `packages/treetime/src/gtr/get_gtr.rs` is `JC69Params`, `K80Params`, `F81Params`, `HKY85Params`, `T92Params`, `Jtt92Params`, `TN93Params`.

1. Rename `struct GTR` to `Gtr` in `packages/treetime/src/gtr/gtr.rs` and every use (rust-analyzer rename or `rg -l '\bGTR\b' packages --type rust` plus `sed -i 's/\bGTR\b/Gtr/g'`, then review the diff for string literals and doc text that name the model, which keep `GTR`).
2. Rename `fn GTR.expQt()` to `exp_qt()` and `fn GTR.expQt_with_rate()` to `exp_qt_with_rate()`. The workspace sets `non_snake_case = "allow"` in `Cargo.toml`, so the compiler does not flag the old names; keep or tighten that setting in a separate change.
3. Rename the parameter structs to `Jc69Params`, `K80Params`, `F81Params`, `Hky85Params`, `T92Params`, `Jtt92Params`, `Tn93Params`. Keep the serialized names (`jc69`, `hky85`, `tn93`) unchanged: check the `#[deser(rename = ...)]` and `#[schemars(...)]` attributes, so that config files and `packages/schemas` stay the same.
4. Search for the old family in docs and generated files: `rg '\bGTR\b|expQt|JC69Params|HKY85Params|TN93Params' packages docs dev`. Keep `GTR` in prose, CLI help text, and citations of v0.

Tests: no new tests (behavior-preserving rename). Oracle: the existing GTR unit and golden master tests in `packages/treetime/src/gtr/__tests__/`.

Validation: `just lint-rs`, `just test-rs`, `just generated-check` (schema text must not change), and the search of step 4 returns only prose and v0 citations.
