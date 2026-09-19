# Rust crate template

Starting point for a new workspace crate. It carries the workspace-inherited
package fields, `[lints] workspace = true`, the `#[ctor]` test setup block that
runs `global_init` and pins Rayon to one thread per test binary, and an empty
`__tests__` module for test files.

To add a crate:

1. Copy the directory to `packages/<crate-name>`.
2. Replace `CRATE_NAME` in `Cargo.toml` with the crate name.
3. Add the crate path to `members` in the workspace `Cargo.toml`.
4. Add real dependencies under `[dependencies]`, each `{ workspace = true }`.

This directory is listed under workspace `exclude`, so cargo never builds the
template in place.
