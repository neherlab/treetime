# Clippy fails on the dylint libraries in `dev/lints`

`just lint-dev-rs` runs clippy with warnings denied in each dylint library of `dev/lints`, each in its own workspace with its pinned nightly ([justfile](../../justfile), recipe `lint-dev-rs`). `check` and `check-all` leave the recipe out. It fails on two of the three libraries:

- `dylint-custom`: 25 warnings in the library and 4 in the `pub-unused-report` binary
- `dylint-mordant`: 271 warnings in the library
- `dylint-trailofbits`: no warnings

Most of the warnings are `clippy::disallowed_types` for `std::collections::HashMap`, `std::collections::HashSet`, and `std::sync::Mutex`. Clippy reads the nearest `clippy.toml` in the crate directory or a parent directory, and the lint libraries have none, so the bans of the root [clippy.toml](../../clippy.toml) apply to them. Those bans keep TreeTime outputs deterministic and steer locks to `parking_lot`.

The other warnings in `dylint-custom` are ordinary clippy findings:

- `clippy::collapsible_if` in `src/lints/topological_ordering.rs` and `src/lints/unclear_exports.rs`
- `clippy::question_mark` in `src/lints/topological_ordering.rs`
- `clippy::arc_with_non_send_sync` in `src/lints/no_comments.rs`
- `clippy::manual_pattern_char_comparison` in `src/lints/prefer_error_macros.rs`

> [!IMPORTANT]
> **Decision required.** Should the root `clippy.toml` bans apply to the lint libraries? The libraries run inside the compiler on rustc types, and their output order does not reach TreeTime results. Options:
>
> - Keep the root bans for the lint libraries and replace the banned types in `dylint-custom` and `dylint-mordant`
> - Add a `clippy.toml` in `dev/lints` without the root bans, then fix the remaining ordinary findings
